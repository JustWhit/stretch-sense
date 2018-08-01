
import javax.xml.xpath.*

Folder = 'Y:\GitRepositories\stretch-sense\Data';
sFolder = '\RawSpirometry';

% % Xiphoid
XMLfile1 = char(fullfile(Folder, sFolder, 'Spiro_6_26_18.xml'));
XMLfile2 = char(fullfile(Folder, sFolder, 'Spiro_5_31_18.xml'));

%MACfilename = '/Users/justinschaffner/Desktop/GitRepositories/stretch-sense/Data/SenseAppData/CAP_2018-03-07542162368_U_R_SIDE.csv';
Filename1 = '\RawSensorData\Abdominal\CAP_2018-06-19_JUSTIN_SVC.csv';
Gfilename1 = '\RawSensorData\Abdominal\GT_2018-06-19_JUSTIN_SVC.csv';
Filename2 = '\RawSensorData\Xiphoid\CAP\CAP\CAP_2018-07-25_JUSTIN_SVC_TEST8.csv';
Gfilename2 = '\RawSensorData\Xiphoid\CAP\GT\GT_2018-07-25JUSTIN_SVC_TEST8.csv';



xmlDoc1 = xmlread(XMLfile1);
xmlDoc2 = xmlread(XMLfile2);

factory = XPathFactory.newInstance;
xpath = factory.newXPath;


FlowData = xpath.compile('//ChannelVolume/SamplingValues');


FlowNodes1 = FlowData.evaluate(xmlDoc1, XPathConstants.NODESET);
FlowNodes2 = FlowData.evaluate(xmlDoc2, XPathConstants.NODESET);

% % extract Spiro traces from two seperate Tests
SpiroTS3 = getSpiroTS(FlowNodes1,1);
SpiroTS2 = getSpiroTS(FlowNodes1,2);
SpiroTS1 = getSpiroTS(FlowNodes1,3);
SpiroTS6 = getSpiroTS(FlowNodes2,1);
SpiroTS5 = getSpiroTS(FlowNodes2,3);
SpiroTS4 = getSpiroTS(FlowNodes2,4);

%extracting capacitance traces from the Stretch Sense data file
capfile1=strcat(Folder,Filename1);
GTfile1 = strcat(Folder,Gfilename1);
capfile2=strcat(Folder,Filename2);
GTfile2 = strcat(Folder,Gfilename2);

% % Read data
CapTb1=readtable(capfile1);
GT1=readtable(GTfile1);
CapTb2=readtable(capfile2);
GT2=readtable(GTfile2);

CapTS1 = getCapTS(CapTb1);
CapTS2 = getCapTS(CapTb2);

% % extract ground truth from ground truth file
GTtime1 = GT1{:,2};
GTtime2 = GT2{:,2};

% % Tlabel = GT{:,1};
GTtime1 = setTimeStamps(GTtime1);
GTtime2 = setTimeStamps(GTtime2);

%extract individual effort traces from capacitance data using the ground
%truth as a guide and resize the cap and spiro samples to the same length for
%correlation
[Cap1,SD1] = getPeaksVals(CapTS1,SpiroTS1,GTtime1,1);
[Cap2,SD2] = getPeaksVals(CapTS1,SpiroTS2,GTtime1,2);
[Cap3,SD3] = getPeaksVals(CapTS1,SpiroTS3,GTtime1,3);
[Cap4,SD4] = getPeaksVals(CapTS2,SpiroTS4,GTtime2,1);
[Cap5,SD5] = getPeaksVals(CapTS2,SpiroTS5,GTtime2,2);
[Cap6,SD6] = getPeaksVals(CapTS2,SpiroTS6,GTtime2,4);

% % Build Data for Model
inputs = [Cap1;Cap2;Cap3;Cap4];
targets = [SD1;SD2;SD3;SD4];
test = [Cap5;Cap6];
actual = [SD5;SD6];

% % Train\Test Model
mdl = fitlm(inputs,targets);
[pred, predCI] = predict(mdl,test);

plotTrainingFit(inputs,targets,1);
plotActualVSPred(actual,pred);
diagnosticPlots(mdl);
residualPlots(mdl);
errorHistogram(pred,actual);

% % returns offset for aligning traces
function a = getOffset(A,Cp,Ctime)
    Ap = A.Data;
    Atime = A.Time;
    start=Atime(1);
    for n = 1:length(Atime)
        Atime(n)=Atime(n)-start;
    end
    
% % Find peaks in data to align traces
    [pks,locs,~,~] = findpeaks(Ap, Atime, 'MinPeakProminence',3);
    [Cpks,Clocs,~,~] = findpeaks(Cp, Ctime, 'MinPeakProminence',0.1);
%     [~,Ai]=max(pks);
%     [~,Ci]=max(Cpks);
    a = Clocs(1)-locs(1);
end

% % returns Offset Capacitance Trace
function a = getCapTrace(CapTS,GTtime,SpiroTS,Sample)
    Cap = getsampleusingtime(CapTS,GTtime(2*Sample),GTtime((2*Sample)+1));
    Spiro = SpiroTS.Data;
    Corr = getOffset(Cap,Spiro,SpiroTS.Time);

    Cap = getsampleusingtime(CapTS,GTtime(2*Sample)-Corr,GTtime((2*Sample)+1));
%     figure,hold on, plot(cumsum(diff(A.Data))), plot(C.^2), title(Sample);
    a = Cap;

end

% % % Pair Spiro and Cap traces
% function a = pairTraces(CapTS,SpiroTS,GTtime,Sample)
%     SD = SpiroTS.Data;
%     C = getCapTrace(CapTS,GTtime,SpiroTS,Sample);
%     length1 = min([numel(C) numel(SD)]);
%     C = C(1:length1);
%     SD = SD(1:length1);
%     a = [C,SD];
% end



% % Returns Peak Valley pairs for Cap and Spiro traces
function [a, b] = getPeaksVals(CapTS, SpiroTS, GTtime,Sample)
    %[C,SD] = pairTraces(CapTS,SpiroTS,GTtime,Sample);
    C = getCapTrace(CapTS,GTtime,SpiroTS,Sample);
    C.Time = setTimeStamps(C.Time);
    % % Find peaks in data
    [Cpks,Clocs,~,~] = findpeaks(C.Data, C.Time, 'MinPeakProminence',3,'MinPeakDistance',2);
    [Cvals,Cvlocs,~,~] = findpeaks(-C.Data,C.Time,'MinPeakProminence',3,'MinPeakDistance',2);
    
    [Spks,Slocs,~,~] = findpeaks(SpiroTS.Data, SpiroTS.Time, 'MinPeakProminence',0.1);
    [Svals,Svlocs,~,~] = findpeaks(-SpiroTS.Data,SpiroTS.Time,'MinPeakProminence',0.1);
    
%     plotPeaksVals(C,Cpks,Clocs,Cvals,Cvlocs,SpiroTS,Spks,Slocs,-Svals,Svlocs);
    length = min([numel(Spks) numel(Svals)])-1; % % subtract 1 to remove the last big Inhale
    Cpks = Cpks(1:length);
    Cvals = Cvals(1:length);
    Spks = Spks(1:length);
    Svals = Svals(1:length);
    plotPeaksVals(C,Cpks,Clocs(1:length),Cvals,Cvlocs(1:length),SpiroTS,Spks,Slocs(1:length),-Svals,Svlocs(1:length));
%     disp(Cvals);
    CapFeatures = abs(Cpks + Cvals);
    SpiroFeatures = abs(Spks + Svals);
    a = CapFeatures;
    b = SpiroFeatures;
end



function a = setTimeStamps(time)
    start = time(1);
    for n = 1:length(time)
        time(n)=time(n)-start;
    end
    a = time;
end


function a = getCapTS(CapT)

    time = CapT{:,2};
    cap = CapT{:,1};
    sRate = 1/100;
    %reset timestamps to start at zero
    time = setTimeStamps(time);

    %resample cap data to match sample rate from spirometer 
    tsvector = 0:sRate:time(length(time));

    ts=timeseries(cap,time);
    a = resample(ts,tsvector);

end


function a = getSpiroTS(FlowNodes, Offset)
    node = FlowNodes.item(FlowNodes.getLength-Offset);
    Spiro = strsplit(char(node.getFirstChild.getNodeValue));
    Spiro = str2double(Spiro);
    Spiro = transpose(Spiro);
    time = linspace(0,length(Spiro)/100,length(Spiro));
    time = transpose(time);
    a=timeseries(Spiro,time);
end

function plotTrainingFit(inputs,targets, polynomial)
    figure,
    plot(inputs,targets,'o');
    grid on; 
    xlabel('Capacitance(pF)');
    ylabel('Volume(L)');
    title('Linear Fit of Training Data Model 3');
    linearCo = polyfit(inputs,targets,polynomial);
    xFit = linspace(min(inputs),max(inputs),50); 
    yFit = polyval(linearCo, xFit); 
    hold on;
    plot(xFit, yFit, 'b--','MarkerSize',1,'LineWidth',1);
    legend('Training Set','Fit','Location','Northwest');

end

function plotActualVSPred(actual,pred)
    figure, 
    hold on, 
    plot(pred,'o'),
    plot(actual,'o'), 
    title('Capacitance with Predicted and Actual Volume'), 
    legend(['Predicted';'Actual   ']),
    xlabel('Time(s)'), 
    ylabel('Volume');

end

function diagnosticPlots(mdl)
    figure,
    subplot(2,1,1), plotDiagnostics(mdl);
    subplot(2,1,2), plotDiagnostics(mdl, 'cookd');
end

function errorHistogram(pred, actual)
    errors = gsubtract(pred,actual);
    figure, ploterrhist(errors);
end

function residualPlots(mdl)
    figure, 
    subplot(3,1,1), plotResiduals(mdl);
    subplot(3,1,2), plotResiduals(mdl, 'probability');
    subplot(3,1,3), plotResiduals(mdl, 'caseorder');
end

function plotPeaksVals(CapTS,Cpks,Clocs,Cvals,Cvlocs,SpiroTS, Spks,Slocs,Svals,Svlocs)
    figure,
    hold on,
    plot(CapTS.Time,CapTS.Data,Clocs,Cpks,'o'); 
    ylabel('Capacitance(pF)');
    xlabel('Time(s)'); 
    grid on; pFREQ = numel(Cpks)/max(Clocs); title(['Peak Detection Freq Estimate = ', num2str(pFREQ), ' Hz; RR est = ', num2str(pFREQ*60),'(1/min)']);
    plot(Cvlocs,-Cvals,'r*');
    plot(SpiroTS.Time,SpiroTS.Data,Slocs,Spks,'o');
    plot(Svlocs,-Svals,'r*');
end