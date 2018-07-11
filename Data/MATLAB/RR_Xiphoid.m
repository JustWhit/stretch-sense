
%% get Respiratory Rate from Signal
import javax.xml.xpath.*

Folder = 'Z:\GitRepositories\stretch-sense\Data';
wFolder = '\Spirometry';


XMLfile1 = char(fullfile(Folder, wFolder, 'Spiro_6_26_18.xml'));
XMLfile2 = char(fullfile(Folder, wFolder, 'Spiro_5_31_18.xml'));

Folder = 'Z:\GitRepositories\stretch-sense\Data\';
%MACfilename = '/Users/justinschaffner/Desktop/GitRepositories/stretch-sense/Data/SenseAppData/CAP_2018-03-07542162368_U_R_SIDE.csv';
Filename1 = 'SenseAppData\Xiphoid\NoVideo\CAP_2018-06-26_JUSTIN_SVC.csv';
Gfilename1 = 'SenseAppData\Xiphoid\NoVideo\GT_2018-06-26_JUSTIN_SVC.csv';
Filename2 = 'SenseAppData\Xiphoid\NoVideo\CAP_2018-05-31_JUSTIN_SVC.csv';
Gfilename2 = 'SenseAppData\Xiphoid\NoVideo\GT_2018-05-31_Justin_SVC.csv';

xmlDoc1 = xmlread(XMLfile1);
xmlDoc2 = xmlread(XMLfile2);

factory = XPathFactory.newInstance;
xpath = factory.newXPath;


FlowData = xpath.compile('//ChannelVolume/SamplingValues');


FlowNodes1 = FlowData.evaluate(xmlDoc1, XPathConstants.NODESET);
FlowNodes2 = FlowData.evaluate(xmlDoc2, XPathConstants.NODESET);

% extracting Spirometer traces from XML file
node3 = FlowNodes1.item(FlowNodes1.getLength-1);
C3 = strsplit(char(node3.getFirstChild.getNodeValue));
C3 = str2double(C3);
C3 = transpose(C3);
Fs = 100;
sRate = 1/Fs;
C3time = linspace(0,length(C3)/100,length(C3));
C3time = transpose(C3time);
tsC3=timeseries(C3,C3time);

node2 = FlowNodes1.item(FlowNodes1.getLength-2);
C2 = strsplit(char(node2.getFirstChild.getNodeValue));
C2 = str2double(C2);
C2 = transpose(C2);
C2time = linspace(0,length(C2)/100,length(C2));
C2time = transpose(C2time);
tsC2=timeseries(C2,C2time);

node1 = FlowNodes1.item(FlowNodes1.getLength-3);
C1 = strsplit(char(node1.getFirstChild.getNodeValue));
C1 = str2double(C1);
C1 = transpose(C1);
C1time = linspace(0,length(C1)/100,length(C1));
C1time = transpose(C1time);
tsC1=timeseries(C1, C1time);

% % Extracting from second XML file
node6 = FlowNodes2.item(FlowNodes2.getLength-1);
C6 = strsplit(char(node6.getFirstChild.getNodeValue));
C6 = str2double(C6);
C6 = transpose(C6);
C6time = linspace(0,length(C6)/100,length(C6));
C6time = transpose(C6time);
tsC6=timeseries(C6,C6time);

node5 = FlowNodes2.item(FlowNodes2.getLength-3);
C5 = strsplit(char(node5.getFirstChild.getNodeValue));
C5 = str2double(C5);
C5 = transpose(C5);
C5time = linspace(0,length(C5)/100,length(C5));
C5time = transpose(C5time);
tsC5=timeseries(C5,C5time);

node4 = FlowNodes2.item(FlowNodes2.getLength-4);
C4 = strsplit(char(node4.getFirstChild.getNodeValue));
C4 = str2double(C4);
C4 = transpose(C4);
C4time = linspace(0,length(C4)/100,length(C4));
C4time = transpose(C4time);
tsC4=timeseries(C4, C4time);


%extracting capacitance traces from the Stretch Sense data file
capfile1=strcat(Folder,Filename1);
GTfile1 = strcat(Folder,Gfilename1);
capfile2=strcat(Folder,Filename2);
GTfile2 = strcat(Folder,Gfilename2);

% Read data
T=readtable(capfile1);
GT=readtable(GTfile1);
T2=readtable(capfile2);
GT2=readtable(GTfile2);

time = T{:,2};
cap = T{:,1};
time2 = T2{:,2};
cap2 = T2{:,1};

%reset timestamps to start at zero
start=time(1);
for n = 1:length(time)
   time(n) = time(n)-start;
end
start2=time2(1);
for n = 1:length(time2)
   time2(n) = time2(n)-start2;
end

%extract ground truth from ground truth file
Ttime = GT{:,2};
Tlabel = GT{:,1};
Ttime2 = GT2{:,2};
Tlabel2 = GT2{:,1};

start=Ttime(1);
for n = 1:length(Ttime)
    Ttime(n)=Ttime(n)-start;
end
start2=Ttime2(1);
for n = 1:length(Ttime2)
    Ttime2(n)=Ttime2(n)-start2;
end

%resample cap data to match sample rate from spirometer
tsvector = 0:sRate:time(length(time));
tsvector2 = 0:sRate:time2(length(time2));

ts=timeseries(cap,time);
ts2=timeseries(cap2,time2);

ts = resample(ts,tsvector);
ts2 = resample(ts2,tsvector2);
    





%extract individual effort traces from capacitance data using the ground
%truth as a guide and resize the two samples to the same length for
%correlation

C1p = tsC1.Data;
A1p = getCapTrace(ts,Ttime,tsC1,1);
length1 = min([numel(A1p) numel(C1p)]);
A1p = A1p(1:length1);
C1p = C1p(1:length1);

A2p= getCapTrace(ts,Ttime,tsC2,2);
C2p = tsC2.Data;
length2 = min([numel(A2p) numel(C2p)]);
A2p = A2p(1:length2);
C2p = C2p(1:length2);

A3p = getCapTrace(ts,Ttime,tsC3,3);
C3p = tsC3.Data;
length3 = min([numel(A3p) numel(C3p)]);
A3p = A3p(1:length3);
C3p = C3p(1:length3);

A4p = getCapTrace(ts2,Ttime2,tsC4,1);
C4p = tsC4.Data;
length4 = min([numel(A4p) numel(C4p)]);
A4p = A4p(1:length4);
C4p = C4p(1:length4);

A5p = getCapTrace(ts2,Ttime2,tsC5,2);
C5p = tsC5.Data;
length5 = min([numel(A5p) numel(C5p)]);
A5p = A5p(1:length5);
C5p = C5p(1:length5);

A6p = getCapTrace(ts2,Ttime2,tsC6,4);
C6p = tsC6.Data;
length6 = min([numel(A6p) numel(C6p)]);
A6p = A6p(1:length6);
C6p = C6p(1:length6);




% % Moving Average Filter
% b = (1/2)*ones(1,2);
% a = 1;
% A1p = filter(b,a,A1p);
% figure; hold on;
% plot(A1p,'r');
% plot(y,'b');
% legend('unfiltered','filtered');
% A1p = A3p;
% %remove DC bias
% meanAmp = mean(A1p);
% A1p = A1p-meanAmp;

% % Hilbert Transform

% y = hilbert(A1p);
% inst_amp= abs(y);
% inst_phase = unwrap(angle(y));
% inst_freq = diff(inst_phase)/(2*pi)*Fs;
% figure; hold on; plot(inst_amp, 'r'); plot(inst_freq, 'b'); plot(cos(inst_phase),'g');
% figure; plot(inst_phase);
% figure;plot(y);
% figure; plot(inst_freq);



% A1p = [0.0; cumsum(diff(A1p))];
% A3p = smooth(A3p,5,'lowess');

% %Truncated 
% A1p = A1p(1:3500);
% C1p = C1p(1:3500);
% C1time = C1time(1:3500);
% A2p = A2p(1:2800);
% C2p = C2p(1:2800);
% C2time = C2time(1:2800);
% A3p = A3p(1:2700);
% C3p = C3p(1:2700);
% C3time = C3time(1:2700);

% % % Pwelch and Peak Plots
% [Pxx,F]=pwelch(A1p,50,[],[],Fs);
% mfreq = meanfreq(Pxx,F);
% figure;
% plot(F,Pxx); ylabel('PSD'); xlabel('Frequency(Hz)'); grid on; [~,loc] = max(Pxx); pwFREQ = F(loc); title(['PWelch Frequency estimate = ', num2str(pwFREQ),' Hz; MeanFreq est = ', num2str(mfreq)]);
% 
% % % Find peaks in data
%     [pks,locs,widths,proms] = findpeaks(A3p, C3time, 'MinPeakProminence',6);
% % % Find Valleys in Data
%     [vals,vlocs,vwidths,vproms] = findpeaks(-A3p,C3time,'MinPeakProminence',6);
%     figure,hold on,
%     plot(C3time,A3p,locs,pks,'o'); ylabel('Capacitance(pF)');xlabel('Time(s)'); grid on; pFREQ = numel(pks)/max(locs); title(['Peak Detection Freq Estimate = ', num2str(pFREQ), ' Hz; RR est = ', num2str(pFREQ*60),'(1/min)']);
%     plot(vlocs,-vals,'r*');


% % %Regression Analysis

% % Create a Fitting Network
% % C1p = C1p - C1p(1,:);
% A1p = A1p - A1p(1,:);
% % C2p = C2p - C2p(1,:);
% A2p = A2p - A2p(1,:);
% % C3p = C3p - C3p(1,:);
% A3p = A3p - A3p(1,:);
% inputs = [A1p;A2p];
% targets = [C1p;C2p];
% hiddenLayerSize = 10;
% net = fitnet(hiddenLayerSize);
% 
% % Set up Division of Data for Training, Validation, Testing
% net.divideParam.trainRatio = 70/100;
% net.divideParam.valRatio = 15/100;
% net.divideParam.testRatio = 15/100;
 
% % Train the Network
% [net,tr] = train(net,inputs,targets);
%  
% % Test the Network
% outputs = net(inputs);
% errors = gsubtract(outputs,targets);
% performance = perform(net,targets,outputs);
%  
% % View the Network
% view(net);



% % %Linear Regression Fit FULL TRACE
% % C1p = C1p - C1p(1,:);
% A1p = A1p - A1p(1,:);
% % C2p = C2p - C2p(1,:);
% A2p = A2p - A2p(1,:);
% 
% 
% % C3p = C3p - C3p(1,:);
% A3p = A3p - A3p(1,:);
% inputs = [A1p;A2p];
% targets = [C1p;C2p];
% mdl = fitlm(inputs,targets);
% [pred, predCI] = predict(mdl,A3p);
% figure, plotDiagnostics(mdl);
% figure, plotDiagnostics(mdl, 'cookd');
% figure, plotResiduals(mdl);
% figure, plotResiduals(mdl, 'probability');
% figure, plotResiduals(mdl, 'caseorder');
% errors = gsubtract(pred,C3p); 
% 
% % Plots
% % % Uncomment these lines to enable various plots.
% % figure, plotperform(tr)
% % figure, plottrainstate(tr)
% % figure, plotfit(net, inputs, targets)
%  figure, plotregression(C3p,pred)
%  figure, ploterrhist(errors)
% NFFT = 2^nextpow2(length(C3p)); % Next power of 2 from length of y
% Y = fft(inputs,NFFT)/length(C3p);
% f = Fs/2*linspace(0,1,NFFT/2+1);
% % Plot single-sided amplitude spectrum.
% figure, plot(f,2*abs(Y(1:NFFT/2+1))) 
% title('Single-Sided Amplitude Spectrum of y(t)')
% xlabel('Frequency (Hz)')
% ylabel('|Y(f)|')
% %I added the next lines to find the value
% %find maximum value, it should be the fundamental frequency (approximated)
% [C,I]= max(2*abs(Y(1:NFFT/2+1)));
% F = f(I);

% %Linear Regression Fit MAX MIN DIFFERENCE

% % Find peaks in data
    [Apks1,Alocs1,Awidths1,Aproms1] = findpeaks(A1p, C1time, 'MinPeakProminence',3);
    [Cpks1,Clocs1,Cwidths1,Cproms1] = findpeaks(C1p, C1time, 'MinPeakProminence',0.1);
    [Apks2,Alocs2,Awidths2,Aproms2] = findpeaks(A2p, C2time, 'MinPeakProminence',3);
    [Cpks2,Clocs2,Cwidths2,Cproms2] = findpeaks(C2p, C2time, 'MinPeakProminence',0.1);
    [Apks3,Alocs3,Awidths3,Aproms3] = findpeaks(A3p, C3time, 'MinPeakProminence',3);
    [Cpks3,Clocs3,Cwidths3,Cproms3] = findpeaks(C3p, C3time, 'MinPeakProminence',0.1);
    [Apks4,Alocs4,Awidths4,Aproms4] = findpeaks(A4p, C4time, 'MinPeakProminence',3);
    [Cpks4,Clocs4,Cwidths4,Cproms4] = findpeaks(C4p, C4time, 'MinPeakProminence',0.1);
    [Apks5,Alocs5,Awidths5,Aproms5] = findpeaks(A5p, C5time, 'MinPeakProminence',3);
    [Cpks5,Clocs5,Cwidths5,Cproms5] = findpeaks(C5p, C5time, 'MinPeakProminence',0.1);
    [Apks6,Alocs6,Awidths6,Aproms6] = findpeaks(A6p, C6time, 'MinPeakProminence',3);
    [Cpks6,Clocs6,Cwidths6,Cproms6] = findpeaks(C6p, C6time, 'MinPeakProminence',0.1);
% % Find Valleys in Data
    [Avals1,Avlocs1,Avwidths1,Avproms1] = findpeaks(-A1p,C1time,'MinPeakProminence',3);
    [Cvals1,Cvlocs1,Cvwidths1,Cvproms1] = findpeaks(-C1p,C1time,'MinPeakProminence',0.1);
    [Avals2,Avlocs2,Avwidths2,Avproms2] = findpeaks(-A2p,C2time,'MinPeakProminence',3);
    [Cvals2,Cvlocs2,Cvwidths2,Cvproms2] = findpeaks(-C2p,C2time,'MinPeakProminence',0.1);
    [Avals3,Avlocs3,Avwidths3,Avproms3] = findpeaks(-A3p,C3time,'MinPeakProminence',3);
    [Cvals3,Cvlocs3,Cvwidths3,Cvproms3] = findpeaks(-C3p,C3time,'MinPeakProminence',0.1);
    [Avals4,Avlocs4,Avwidths4,Avproms4] = findpeaks(-A4p,C4time,'MinPeakProminence',3);
    [Cvals4,Cvlocs4,Cvwidths4,Cvproms4] = findpeaks(-C4p,C4time,'MinPeakProminence',0.1);
    [Avals5,Avlocs5,Avwidths5,Avproms5] = findpeaks(-A5p,C5time,'MinPeakProminence',3);
    [Cvals5,Cvlocs5,Cvwidths5,Cvproms5] = findpeaks(-C5p,C5time,'MinPeakProminence',0.1);
    [Avals6,Avlocs6,Avwidths6,Avproms6] = findpeaks(-A6p,C6time,'MinPeakProminence',3);
    [Cvals6,Cvlocs6,Cvwidths6,Cvproms6] = findpeaks(-C6p,C6time,'MinPeakProminence',0.1);
    
    
%     figure,hold on,plot(C1time,A1p,Alocs1,Apks1,'o'); ylabel('Capacitance(pF)');xlabel('Time(s)'); grid on; pFREQ = numel(Apks1)/max(Alocs1); title(['Peak Detection Freq Estimate = ', num2str(pFREQ), ' Hz; RR est = ', num2str(pFREQ*60),'(1/min)']);
% %     plot(vlocs,-vals,'r*');
%     plot(C1time,C1p,Clocs1,Cpks1,'r*');


% % Pair Peaks to Valleys
lengthA = min([numel(Apks1) numel(Avals1)])-1;
Apks1 = Apks1(1:lengthA);
Avals1 = Avals1(1:lengthA);
Cpks1 = Cpks1(1:lengthA);
Cvals1 = Cvals1(1:lengthA);

lengthA = min([numel(Apks2) numel(Avals2)])-1;
Apks2 = Apks2(1:lengthA);
Avals2 = Avals2(1:lengthA);
Cpks2 = Cpks2(1:lengthA);
Cvals2 = Cvals2(1:lengthA);

lengthA = min([numel(Apks3) numel(Avals3)])-1;
Apks3 = Apks3(1:lengthA);
Avals3 = Avals3(1:lengthA);
Cpks3 = Cpks3(1:lengthA);
Cvals3 = Cvals3(1:lengthA);

lengthA = min([numel(Apks4) numel(Avals4)]);
Apks4 = Apks4(1:lengthA);
Avals4 = Avals4(1:lengthA);
Cpks4 = Cpks4(1:lengthA);
Cvals4 = Cvals4(1:lengthA);

lengthA = min([numel(Apks5) numel(Avals5)]);
Apks5 = Apks5(1:lengthA);
Avals5 = Avals5(1:lengthA);
Cpks5 = Cpks5(1:lengthA);
Cvals5 = Cvals5(1:lengthA);

lengthA = min([numel(Apks6) numel(Avals6)]);
Apks6 = Apks6(1:lengthA);
Avals6 = Avals6(1:lengthA);
Cpks6 = Cpks6(1:lengthA);
Cvals6 = Cvals6(1:lengthA);


% % Build Data for Model
inputs = [abs(Apks1 + Avals1);abs(Apks2 + Avals2);abs(Apks3 + Avals3);abs(Apks4 + Avals4)];
targets = [abs(Cpks1 + Cvals1); abs(Cpks2 + Cvals2);abs(Cpks3 + Cvals3);abs(Cpks4 + Cvals4)];
test = [abs(Apks5 + Avals5);abs(Apks6 + Avals6)];
actual = [abs(Cpks5 + Cvals5);abs(Cpks6 + Cvals6)];

% % Train\Test Model
mdl = fitlm(inputs,targets);
[pred, predCI] = predict(mdl,test);

% % Plots
figure, plotDiagnostics(mdl);
figure, plotDiagnostics(mdl, 'cookd');
figure, plotResiduals(mdl);
figure, plotResiduals(mdl, 'probability');
figure, plotResiduals(mdl, 'caseorder');
errors = gsubtract(pred,actual); 

 figure, plotregression(actual,pred)
 figure, ploterrhist(errors)
NFFT = 2^nextpow2(length(actual)); % Next power of 2 from length of y
Y = fft(inputs,NFFT)/length(actual);
f = Fs/2*linspace(0,1,NFFT/2+1);
% Plot single-sided amplitude spectrum.
figure, plot(f,2*abs(Y(1:NFFT/2+1))) 
title('Single-Sided Amplitude Spectrum of y(t)')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
%I added the next lines to find the value
%find maximum value, it should be the fundamental frequency (approximated)
[C,I]= max(2*abs(Y(1:NFFT/2+1)));
F = f(I);

figure, hold on, plot(test,'o'), plot(pred,'o'),plot(actual,'o'), title('Capacitance with Predicted and Actual Volume'), legend(['Capacitance(pF)';'Pvolume(L)     ';'Avolume(L)     ']),xlabel('Data'), ylabel('Value');




% % returns offset for aligning traces
function a = getOffset(A,Cp,Ctime)
    Ap = A.Data;
    Atime = A.Time;
    start=Atime(1);
    for n = 1:length(Atime)
        Atime(n)=Atime(n)-start;
    end
    
% % Find peaks in data to align traces on Biggest Peak (Deep Breath)
    [pks,locs,~,~] = findpeaks(Ap, Atime, 'MinPeakProminence',3);
    [Cpks,Clocs,~,~] = findpeaks(Cp, Ctime, 'MinPeakProminence',0.25);
%     [~,Ai]=max(pks);
%     [~,Ci]=max(Cpks);
    a = Clocs(1)-locs(1);
end

% % returns Offset Capacitance Trace
function a = getCapTrace(CapTS,GTtime,SpiroTS,Sample)
    A = getsampleusingtime(CapTS,GTtime(2*Sample),GTtime((2*Sample)+1));
    C = SpiroTS.Data;
    Corr = getOffset(A,C,SpiroTS.Time);

    A = getsampleusingtime(CapTS,GTtime(2*Sample)-Corr,GTtime((2*Sample)+1));
%     figure,hold on, plot(cumsum(diff(A.Data))), plot(C.^2), title(Sample);
    a = A.Data;

end