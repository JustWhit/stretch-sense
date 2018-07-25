import javax.xml.xpath.*

Folder = 'Z:\GitRepositories\stretch-sense\Data';
sFolder = '\RawSpirometry';
dFolder = '\RawSensorData';
soutFolder = '\ExSpirometry';
outFolder = '\ExSensorData';

% % Assign test NAME HERE
TestName = '7_25_18_JUSTIN_SVC_TEST8';
XMLfile = char(fullfile(Folder, sFolder, 'Spiro_7_25_18.xml'));

capFile = '';
noteFile = ' ';
% % XML extraction
xmlDoc = xmlread(XMLfile);

factory = XPathFactory.newInstance;
xpath = factory.newXPath;

FlowData = xpath.compile('//ChannelVolume/SamplingValues');

FlowNodes = FlowData.evaluate(xmlDoc, XPathConstants.NODESET);

for i = 1:5
csvwrite(char(fullfile(Folder, strcat(TestName,'_T', num2str((5-i)+1),'.csv'))),getSpiro(FlowNodes,i));
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


% % returns Offset Capacitance Trace
function a = getCapTrace(CapTS,GTtime,SpiroTS,Sample)
    Cap = getsampleusingtime(CapTS,GTtime(2*Sample),GTtime((2*Sample)+1));
    Spiro = SpiroTS.Data;
    Corr = getOffset(Cap,Spiro,SpiroTS.Time);

    Cap = getsampleusingtime(CapTS,GTtime(2*Sample)-Corr,GTtime((2*Sample)+1));
%     figure,hold on, plot(cumsum(diff(A.Data))), plot(C.^2), title(Sample);
    a = Cap;

end


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

function Spiro = getSpiro(FlowNodes, Offset)
    node = FlowNodes.item(FlowNodes.getLength-Offset);
    S = strsplit(char(node.getFirstChild.getNodeValue));
    S = str2double(S);
    S = transpose(S);
    time = linspace(0,length(S)/100,length(S));
    time = transpose(time);
    Spiro = [S time];
    figure;plot(time,S);title(num2str((5-Offset)+1));
end