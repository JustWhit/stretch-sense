
%% get Respiratory Rate from Signal
import javax.xml.xpath.*

Folder = 'Z:\GitRepositories\stretch-sense\Data';
wFolder = '\Spirometry';


XMLfile = char(fullfile(Folder, wFolder, 'Spiro_6_19_18.xml'));

Folder = 'Z:\GitRepositories\stretch-sense\Data\';
%MACfilename = '/Users/justinschaffner/Desktop/GitRepositories/stretch-sense/Data/SenseAppData/CAP_2018-03-07542162368_U_R_SIDE.csv';
Filename = 'SenseAppData\Abdominal\CAP_2018-06-19_JUSTIN_SVC.csv';
Gfilename = 'SenseAppData\Abdominal\GT_2018-06-19_JUSTIN_SVC.csv';

xmlDoc = xmlread(XMLfile);

factory = XPathFactory.newInstance;
xpath = factory.newXPath;

FlowData = xpath.compile('//ChannelVolume/SamplingValues');

FlowNodes = FlowData.evaluate(xmlDoc, XPathConstants.NODESET);

% extracting Spirometer traces from XML file
node3 = FlowNodes.item(FlowNodes.getLength-1);
C3 = strsplit(char(node3.getFirstChild.getNodeValue));
C3 = str2double(C3);
C3 = transpose(C3);
Fs = 100;
sRate = 1/Fs;
C3time = linspace(0,length(C3)/100,length(C3));
C3time = transpose(C3time);
tsC3=timeseries(C3,C3time);

node2 = FlowNodes.item(FlowNodes.getLength-2);
C2 = strsplit(char(node2.getFirstChild.getNodeValue));
C2 = str2double(C2);
C2 = transpose(C2);
C2time = linspace(0,length(C2)/100,length(C2));
C2time = transpose(C2time);
tsC2=timeseries(C2,C2time);

node1 = FlowNodes.item(FlowNodes.getLength-4);
C1 = strsplit(char(node1.getFirstChild.getNodeValue));
C1 = str2double(C1);
C1 = transpose(C1);
C1time = linspace(0,length(C1)/100,length(C1));
C1time = transpose(C1time);
tsC1=timeseries(C1, C1time);


%extracting capacitance traces from the Stretch Sense data file
capfile=strcat(Folder,Filename);
GTfile = strcat(Folder,Gfilename);

% Read data
T=readtable(capfile);
GT=readtable(GTfile);

time = T{:,2};
cap = T{:,1};

%reset timestamps to start at zero
start=time(1);
for n = 1:length(time)
   time(n) = time(n)-start;
end

%extract ground truth from ground truth file
Ttime = GT{:,2};
Tlabel = GT{:,1};
% for n = 1:length(Ttime)
%    Ttime(n) = (100 * (Ttime(n) - floor(Ttime(n)))) + (floor(Ttime(n)) * 60);
% end
% start = cap(1);
% for i=1:length(cap)
%     cap(i)= cap(i)-start;
% end
start=Ttime(1);
for n = 1:length(Ttime)
    Ttime(n)=Ttime(n)-start;
end

%resample cap data to match sample rate from spirometer
Fs = 100;
sRate = 1/Fs;
tsvector = 0:sRate:time(length(time));

ts=timeseries(cap,time);

ts = resample(ts,tsvector);
    
%     for i=2:2:length(Tlabel)
%         A(i) = getsampleusingtime(ts, Ttime(i),Ttime(i+1));
%     end





%extract individual effort traces from capacitance data using the ground
%truth as a guide and resize the two samples to the same length for
%correlation


A1 = getsampleusingtime(ts,Ttime(2)-0.5,Ttime(3));

A1p = A1.Data;
C1p = tsC1.Data;
length1 = min([numel(A1p) numel(C1p)]);
A1p = A1p(1:length1);
C1p = C1p(1:length1);

A2=getsampleusingtime(ts,Ttime(4)-0.6,Ttime(5));
A2p = A2.Data;
C2p = tsC2.Data;
length2 = min([numel(A2p) numel(C2p)]);
A2p = A2p(1:length2);
C2p = C2p(1:length2);

A3 = getsampleusingtime(ts,Ttime(6)-0.67,Ttime(7));
A3p = A3.Data;
C3p = tsC3.Data;
length3 = min([numel(A3p) numel(C3p)]);
A3p = A3p(1:length3);
C3p = C3p(1:length3);

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

% %Truncated Abdominal
% A1p = A1p(1:3500);
% C1p = C1p(1:3500);
% C1time = C1time(1:3500);
% A2p = A2p(1:2800);
% C2p = C2p(1:2800);
% C2time = C2time(1:2800);
% A3p = A3p(1:2700);
% C3p = C3p(1:2700);
% C3time = C3time(1:2700);

[Pxx,F]=pwelch(A1p,50,[],[],Fs);
mfreq = meanfreq(Pxx,F);
figure;
plot(F,Pxx); ylabel('PSD'); xlabel('Frequency(Hz)'); grid on; [~,loc] = max(Pxx); pwFREQ = F(loc); title(['PWelch Frequency estimate = ', num2str(pwFREQ),' Hz; MeanFreq est = ', num2str(mfreq)]);

% % Find peaks in data
    [pks,locs,widths,proms] = findpeaks(A3p, C3time, 'MinPeakProminence',6);
% % Find Valleys in Data
    [vals,vlocs,vwidths,vproms] = findpeaks(-A3p,C3time,'MinPeakProminence',6);
    figure,hold on,
    plot(C3time,A3p,locs,pks,'o'); ylabel('Capacitance(pF)');xlabel('Time(s)'); grid on; pFREQ = numel(pks)/max(locs); title(['Peak Detection Freq Estimate = ', num2str(pFREQ), ' Hz; RR est = ', num2str(pFREQ*60),'(1/min)']);
    plot(vlocs,-vals,'r*');
% %Regression Analysis


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
    [Apks1,Alocs1,Awidths1,Aproms1] = findpeaks(A1p, C1time, 'MinPeakProminence',6);
    [Cpks1,Clocs1,Cwidths1,Cproms1] = findpeaks(C1p, C1time, 'MinPeakProminence',0.25);
    [Apks2,Alocs2,Awidths2,Aproms2] = findpeaks(A2p, C2time, 'MinPeakProminence',6);
    [Cpks2,Clocs2,Cwidths2,Cproms2] = findpeaks(C2p, C2time, 'MinPeakProminence',0.25);
    [Apks3,Alocs3,Awidths3,Aproms3] = findpeaks(A3p, C3time, 'MinPeakProminence',6);
    [Cpks3,Clocs3,Cwidths3,Cproms3] = findpeaks(C3p, C3time, 'MinPeakProminence',0.25);
% % Find Valleys in Data
    [Avals1,Avlocs1,Avwidths1,Avproms1] = findpeaks(-A1p,C1time,'MinPeakProminence',6);
    [Cvals1,Cvlocs1,Cvwidths1,Cvproms1] = findpeaks(-C1p,C1time,'MinPeakProminence',0.25);
    [Avals2,Avlocs2,Avwidths2,Avproms2] = findpeaks(-A2p,C2time,'MinPeakProminence',6);
    [Cvals2,Cvlocs2,Cvwidths2,Cvproms2] = findpeaks(-C2p,C2time,'MinPeakProminence',0.25);
    [Avals3,Avlocs3,Avwidths3,Avproms3] = findpeaks(-A3p,C3time,'MinPeakProminence',6);
    [Cvals3,Cvlocs3,Cvwidths3,Cvproms3] = findpeaks(-C3p,C3time,'MinPeakProminence',0.25);
    
    
    figure,hold on,
    plot(C1time,C1p,Clocs1,Cpks1,'o'); ylabel('Capacitance(pF)');xlabel('Time(s)'); grid on; pFREQ = numel(Cpks1)/max(Clocs1); title(['Peak Detection Freq Estimate = ', num2str(pFREQ), ' Hz; RR est = ', num2str(pFREQ*60),'(1/min)']);
    plot(Cvlocs1,-Cvals1,'r*');
    
% % Pair Peaks to Valleys
lengthA = min([numel(Apks1) numel(Avals1)]);
Apks1 = Apks1(1:lengthA);
Avals1 = Avals1(1:lengthA);
Cpks1 = Cpks1(1:lengthA);
Cvals1 = Cvals1(1:lengthA);

lengthA = min([numel(Apks2) numel(Avals2)]);
Apks2 = Apks2(1:lengthA);
Avals2 = Avals2(1:lengthA);
Cpks2 = Cpks2(1:lengthA);
Cvals2 = Cvals2(1:lengthA);

lengthA = min([numel(Apks3) numel(Avals3)]);
Apks3 = Apks3(1:lengthA);
Avals3 = Avals3(1:lengthA);
Cpks3 = Cpks3(1:lengthA);
Cvals3 = Cvals3(1:lengthA);

% % Build Data for Model
inputs = [abs(Apks1 + Avals1);abs(Apks2 + Avals2)];
targets = [abs(Cpks1 + Cvals1); abs(Cpks2 + Cvals2)];
test = abs(Apks3 + Avals3);
actual = abs(Cpks3 + Cvals3);

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
