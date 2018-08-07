


TFolder = 'Y:\GitRepositories\stretch-sense\Data\ModelData';
FileListT = dir(char(fullfile(TFolder,'*.csv')));
% % Xiphoid


% % Build Data for Model
Features = [];
for n=1:length(FileListT)
   Table = readtable(char(fullfile(TFolder,FileListT(n).name))); 
   disp(FileListT(n).name);
   Time = [Table{:,1}];
   Sensors = [Table{:,2:3}];
   SensorTS = timeseries(Sensors,Time);
   [C,S] = getPeaksVals(SensorTS);
   % % % Make changes her to adjust the Features Used by the Model
   
%    % % Difference between Peaks and Valleys (inhale) minus deep breath
%     CapF = abs(C(:,1) + C(:,3));
%     SpiroF = abs(S(:,1) + S(:,3));
%     Features = [Features; [CapF SpiroF]];

   % % Values between Peaks and Valleys (inhale)
    for i=1:length(C(:,1))
        Features = [Features; [SensorTS.Data(C(i,4):C(i,2),1) SensorTS.Data(S(i,4):S(i,2),2)]];
    end
    
    
    
end
Features = shuffleRow(Features);
CapFeatures = Features(:,1);
SpiroFeatures = Features(:,2);
TRsetSize = floor(.8 * length(CapFeatures));
inputs = CapFeatures(1:TRsetSize);
targets = SpiroFeatures(1:TRsetSize);
test = CapFeatures(TRsetSize+1:end);
actual = SpiroFeatures(TRsetSize+1:end);


% % Train\Test Model
mdl = fitlm(inputs,targets);
[pred, predCI] = predict(mdl,test);

plotTrainingFit(inputs,targets,1);
plotActualVSPred(actual,pred);
% diagnosticPlots(mdl);
% residualPlots(mdl);
% errorHistogram(pred,actual);
plotScatter(inputs,targets);



% % Returns teh percent difference between two values
function a = pDiff(v1, v2)
    a = (abs(v1-v2)/((v1+v2)/2))*100;
end

% % Shuffles the Data
function ret = shuffleRow(mat)

    [r c] = size(mat);
    shuffledRow = randperm(r);
    ret = mat(shuffledRow, :);
end

% % Returns Peak Valley pairs for Cap and Spiro traces
function [a, b] = getPeaksVals(SensorsTS)
    
    % % Find peaks in data
    [Cpks,Clocs,~,~] = findpeaks(SensorsTS.Data(:,1), SensorsTS.Time, 'MinPeakProminence',.75,'MinPeakDistance',1);
    [Cvals,Cvlocs,~,~] = findpeaks(-SensorsTS.Data(:,1),SensorsTS.Time,'MinPeakProminence',.75,'MinPeakDistance',1);
    
    [Spks,Slocs,~,~] = findpeaks(SensorsTS.Data(:,2), SensorsTS.Time, 'MinPeakProminence',0.1);
    [Svals,Svlocs,~,~] = findpeaks(-SensorsTS.Data(:,2),SensorsTS.Time,'MinPeakProminence',0.1);
    
%     plotPeaksVals(C,Cpks,Clocs,Cvals,Cvlocs,SpiroTS,Spks,Slocs,-Svals,Svlocs);
    length = min([numel(Spks) numel(Svals)])-1; % % subtract 1 to remove the last big Inhale
    Cpks = Cpks(1:length);
    Clocs = Clocs(1:length);
    Cvals = Cvals(1:length);
    Cvlocs = Cvlocs(1:length);
    Spks = Spks(1:length);
    Slocs = Slocs(1:length);
    Svals = Svals(1:length);
    Svlocs = Svlocs(1:length);
    plotPeaksVals(SensorsTS,Cpks,Clocs(1:length),Cvals,Cvlocs(1:length),Spks,Slocs(1:length),Svals,Svlocs(1:length));
%     disp(Cvals);
    a = [Cpks Clocs Cvals Cvlocs];
    b = [Spks Slocs Svals Svlocs];
end



function a = setTimeStamps(time)
    start = time(1);
    for n = 1:length(time)
        time(n)=time(n)-start;
    end
    a = time;
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
    title('Predicted and Actual Volume'), 
    legend({'Predicted', 'Actual'}),
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

function plotPeaksVals(SensorTS,Cpks,Clocs,Cvals,Cvlocs, Spks,Slocs,Svals,Svlocs)
    figure,
    hold on,
    plot(SensorTS.Time,SensorTS.Data(:,1),Clocs,Cpks,'o'); 
    ylabel('Capacitance(pF)');
    xlabel('Time(s)'); 
    grid on; pFREQ = numel(Cpks)/max(Clocs); title(['Peak Detection Freq Estimate = ', num2str(pFREQ), ' Hz; RR est = ', num2str(pFREQ*60),'(1/min)']);
    plot(Cvlocs,-Cvals,'r*');
    plot(SensorTS.Time,SensorTS.Data(:,2),Slocs,Spks,'o');
    plot(Svlocs,-Svals,'r*');
end

function plotScatter(Capacitance, Spirometer)
    figure;
    scatter(Capacitance,Spirometer);
    title('Scatter Plot of Training Data');
    xlabel('Capacitance(pF)');
    ylabel('Volume(L)');
end