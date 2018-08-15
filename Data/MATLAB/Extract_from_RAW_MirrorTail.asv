

Folder = 'Y:\GitRepositories\stretch-sense\Data';
sFolder = '\RawSpirometry';
dFolder = '\RawSensorData';
outFolder = '\ExSensorSpiroData';

% % Assign input and output files here
TestName = '08_06_18_JUSTIN_SVC_TEST9';
SpiroFiles = {'SVCVOLTest9T1_08_06_18.csv' 'SVCVOLTest9T2_08_06_18.csv' 'SVCVOLTest10T3_08_07_18.csv'};
capFile=char(fullfile(Folder,dFolder,'\Xiphoid\CAP_HEART_SPIRO\CAP_2018-08-07_TEST_10_Justin.csv'));
noteFile = char(fullfile(Folder,dFolder,'\Xiphoid\CAP_HEART_SPIRO\GT_2018-08-07_TEST_10_Justin.csv'));

SpiroTraces = {};
for n=1: length(SpiroFiles)
   temp = csvread(char(fullfile(Folder,sFolder,SpiroFiles{n})));
   SpiroTraces{n} = temp;
end


%Read Cap and Note files
capTbl = readtable(capFile);
noteTbl = readtable(noteFile);
noteTime = noteTbl{:,2};
noteLabel = noteTbl{:,1};

%Get Cap TimeSeries
[WSensorTS,noteTime] = getWholeSensorTS(capTbl, noteTime, noteLabel);


for i = 1:length(SpiroFiles)
    Sensors = extract(SpiroTraces, WSensorTS,noteTime,i);
    csvwrite(char(fullfile(Folder, outFolder,strcat(TestName,'_T', num2str(i),'.csv'))),Sensors); 
end


% % Main Function, Pairs Spirometer Data with Cap trace, saves result to a
% % CSV file. If Heart Rate or Audio Amplitude is included, make
% % adjustments where noted for each function
function a = extract(SpiroTraces, WSensorTS, NoteTime, Sample)
    SpiroTS = getSpiroTS(SpiroTraces,Sample);
    SensorTS = getSensorTraceTS(WSensorTS,NoteTime,SpiroTS,(Sample));
    % [Cap,Heart,AudioAmp] = getSensorTraces(WCapTS,NoteTime,Sample);
    length = min([numel(SensorTS.Data(:,1)) numel(SpiroTS.Data)]);
    Time = SpiroTS.Time(1:length);
%     Cap = SensorTS.Data(1:length);
    Cap = SensorTS.Data(1:length,1);
    Heart = SensorTS.Data(1:length,2);
%     AudioAmp = SensorTS.Data(1:length,1);
    Spiro = SpiroTS.Data(1:length);
%     a = [Time Cap Spiro]; 
    a =  [Time Cap Spiro Heart]; % AudioAmp];
    figure;hold on;plot(Time,Cap);plot(Time,Spiro.^2);title(num2str(Sample)); 
end

% % Finds the Exhale tail and flips it
function a = flipExhaleTail(CapTrace)
figure; hold on;
plot(CapTrace);
    [Cpks,Clocs,~,~] = findpeaks(CapTrace, 'MinPeakProminence',3,'MinPeakDistance',200);
    [Cvals,Cvlocs,~,~] = findpeaks(-CapTrace, 'MinPeakProminence',3,'MinPeakDistance',200);
    % find the peak index of the peak of inhale
    Mx = max(Cpks);
    assignin('base','Max',Mx);
    i = find(Cpks==Mx,1);
plot(Clocs(i),Cpks(i),'r*');
    % find the position in the trace of peak of inhale
    position = Clocs(i);
    % find the first valley after that peak in time
    firstindex = find(Cvlocs>position,1);
    firstVal = Cvlocs(firstindex);
    secondVal = Cvlocs(firstindex + 1);
    % get the segment from frist valley to second valley
    Segment = CapTrace(firstVal:secondVal,1);
    % mirror the segment
    Segment = mirror(Segment);
    % replace teh original segment with the flipped segment
    CapTrace(firstVal:secondVal,1) = Segment;
    a = CapTrace;
plot(CapTrace);
end


% % Mirrors a segment of data
function a = mirror(Segment)
% figure; hold on;
%     plot(Segment);
    Segment = transpose(Segment);
    assignin('base','SegmentA',Segment);
    
    % creat an x variable for the matrix operation
    x = linspace(1,numel(Segment),numel(Segment));
    % define variables for the slope of the line and the offset
    m = (Segment(end)-Segment(1))/(x(end)-x(1));
    b = Segment(end) - (m*x(end));
    % Create reflection matrix
    T = (1/(1+m.^2))*[[1-m.^2,2*m,-2*m*b];[2*m,m.^2-1,2*b];[0,0,1+m.^2]];
    assignin('base','T',T);
    assignin('base','x',x);
    assignin('base','m',m);
    assignin('base','b',b);
    Segment = [x ; Segment; (1 + zeros(1,numel(Segment)))];
    assignin('base','Segment',Segment);
    Segment = T*Segment;
    assignin('base','SegmentF',Segment);
    
    
    a = transpose(Segment(2,:));
%     plot(a);
end


% % % Pair Spiro and Cap traces
% function  = matchCapTrace(CapTS,SpiroTS,NoteTime,Sample)
%     SD = SpiroTS.Data;
%     C = getTraceCapTS(CapTS,NoteTime,SpiroTS,Sample);
%     length1 = min([numel(C) numel(SD)]);
%     C = C(1:length1);
%     SD = SD(1:length1);
%     
% end


function [time, noteTime] = setWholeTimeStamps(time, noteTime)
    start = min(time(1),noteTime(1));
    time = time - start;
    for n = 1:length(noteTime)
        noteTime(n)=noteTime(n)-start;
    end
end

function time = setTimeStamps(time)
    start = time(1);
    time = time - start;
end


function [ts, noteTime] = getWholeSensorTS(SenseTbl, noteTime, noteLabel)
    % % % % IF FILE HAS HEART RATE, ADJUST COLUMNS RIGHT BY 1, ADD HEART
    % % % % RATE VARIABLE
%     time = SenseTbl{:,2};
%     length(time)
%     cap = SenseTbl{:,1};

%     % For files with more sensors
    time = SenseTbl{:,3};
    cap = SenseTbl{:,2};
    heart = SenseTbl{:,1};
%     audioAmp = SenseTbl{:,1};
    Sensors = [cap heart];

    sRate = 1/100;
    %reset timestamps to start at zero
    [time, noteTime] = setWholeTimeStamps(time, noteTime);
    
    %resample cap data to match sample rate from spirometer
    tsvector = time(1):sRate:time(length(time));

%     ts=timeseries(cap,time);
    ts = timeseries(Sensors,time); % Use this instead when file has
% %     multiple sensors
    ts = resample(ts,tsvector);
    figure; hold on; plot(ts); title('Whole TS');
     linetype = {'g--','r--'}; %'g--','r--','c--','m--'};
    for n=1:length(noteTime)
          text(noteTime(n),min(ylim),noteLabel(n),'Rotation',90,'Color',[1,0,0]);
          vline(noteTime(n),linetype{2});
    end
    hold off;
end


% % returns Offset Capacitance Trace
function a = getSensorTraceTS(SenseTS,NoteTime,SpiroTS,Sample)
    Sensors = getsampleusingtime(SenseTS,NoteTime(2*Sample),NoteTime((2*Sample)+1));
    
    Sensors.Time = setTimeStamps(Sensors.Time);
    Spiro = SpiroTS.Data;
    Corr = getOffset(Sensors,Spiro,SpiroTS.Time);
%     Corr
%     (NoteTime(2*Sample)-Corr)
%     
    
    Sensors = getsampleusingtime(SenseTS,(NoteTime(2*Sample)-Corr),NoteTime((2*Sample)+1)+8);
    Sensors.Data(:,1) = flipExhaleTail(Sensors.Data(:,1));
    figure,hold on, plot(cumsum(diff(Sensors.Data(:,2)))), plot(Spiro.^2), title(Sample);
    Sensors.Time = setTimeStamps(Sensors.Time);
    a = Sensors;
     
end


% % returns offset for aligning traces
function a = getOffset(SensorTS,Spiro,Spirotime)
%     Cap = CapTS.Data;
    Cap = SensorTS.Data(:,1); % Use if CapTS has multiple sensors
    time = SensorTS.Time;
%     time = setTimeStamps(CapTS.Time);
    
    
% % Find peaks in data to align traces
    [pks,locs,~,~] = findpeaks(Cap, time, 'MinPeakProminence',3,'MinPeakDistance',1);
    [Spks,Slocs,~,~] = findpeaks(Spiro, Spirotime, 'MinPeakProminence',0.1);
    figure;plot(Spirotime,Spiro);
    figure; hold on; plot(SensorTS); plot(locs,pks,'r*');plot(Spirotime,Spiro);plot(Slocs,Spks,'r*');
%     [~,Ci]=max(pks);
%     [~,Si]=max(spks);
% Slocs(1)
% locs(1)
    a = Slocs(1)-locs(1);
end


function Spiro = getSpiroTS(SpiroTraces, Sample)
    S = SpiroTraces{Sample};
    S = transpose(S);
    time = linspace(0,length(S)/100,length(S));
    time = transpose(time);
    Spiro = timeseries(S, time);
    
end