


import javax.xml.xpath.*

Folder = 'Z:\GitRepositories\stretch-sense\Data';
wFolder = '\Spirometry';

XMLfile = char(fullfile(Folder, wFolder, 'Spiro_XML_export.xml'));

xmlDoc = xmlread(XMLfile);

factory = XPathFactory.newInstance;
xpath = factory.newXPath;

FlowData = xpath.compile('//ChannelVolume/SamplingValues');

FlowNodes = FlowData.evaluate(xmlDoc, XPathConstants.NODESET);


node3 = FlowNodes.item(FlowNodes.getLength-3);
C3 = strsplit(char(node3.getFirstChild.getNodeValue));
C3 = str2double(C3);

node2 = FlowNodes.item(FlowNodes.getLength-4);
C2 = strsplit(char(node2.getFirstChild.getNodeValue));
C2 = str2double(C2);

node1 = FlowNodes.item(FlowNodes.getLength-5);
C1 = strsplit(char(node1.getFirstChild.getNodeValue));
C1 = str2double(C1);

plot(C1);


%%
%% Single Plot

Folder = 'Z:\GitRepositories\stretch-sense\Data\';
%MACfilename = '/Users/justinschaffner/Desktop/GitRepositories/stretch-sense/Data/SenseAppData/CAP_2018-03-07542162368_U_R_SIDE.csv';
Filename = 'SenseAppData\Xiphoid\CAP\CAP_2018-04-24_FVC_EX_Justin.csv'
Gfilename = 'SenseAppData\Xiphoid\GT\GT_2018-04-24_FVC_EX_Justin.csv'

%GTfile='/Users/justinschaffner/Desktop/GitRepositories/stretch-sense/Data/SenseAppData/GT_2018-03-07542161920_U_F_CENTER.csv';
%SGTfile='/Users/justinschaffner/Desktop/GitRepositories/stretch-sense/Data/SenseAppData/GT_2018-03-07542162176_U_R_CORNER.csv';
%capfiles = ["Xiphoid\CAP_2018-03-17_X_B_CENTER_EILEEN_31_5inch.csv","4thIntercostal\CAP_2018-03-17_4th_B_Center_EILEEN_34_5inch.csv"];
%GTfiles = ["Xiphoid\GT_2018-03-17_X_B_CENTER_EILEEN_31_5inch.csv","4thIntercostal\GT_2018-03-17_4th_B_Center_EILEEN_34_5inch.csv"];



capfile=strcat(Folder,Filename);
GTfile = strcat(Folder,Gfilename);

%Read data
T=readtable(capfile);
GT=readtable(GTfile);

time = T{:,2};
cap = T{:,1};

start=time(1);
for n = 1:length(time)
   time(n) = time(n)-start;
end

ts=timeseries(cap,time);

Ttime = GT{:,2};
Tlabel = GT{:,1};
%for n = 1:length(Ttime)
%    Ttime(n) = (100 * (Ttime(n) - floor(Ttime(n)))) + (floor(Ttime(n)) * 60);
% end
start = cap(1);
for i=1:length(cap)
    cap(i)= cap(i)-start;
end
start=Ttime(1);
for n = 1:length(Ttime)
    Ttime(n)=Ttime(n)-start;
end

    
%     for i=2:2:length(Tlabel)
%         A(i) = getsampleusingtime(ts, Ttime(i),Ttime(i+1));
%     end

figure; hold on;
A1 = getsampleusingtime(ts,Ttime(2),Ttime(3));

% A2=getsampleusingtime(ts,Ttime(7)-1,Ttime(8));
% plot(A2);



plot(A1);
% A3 = getsampleusingtime(ts,Ttime(8),Ttime(9));
% plot(A3);
    

    %  Plot data
%     plot(ts)
    
    % Find peaks in data
%     [pks,locs,widths,proms] = findpeaks(cap, time,'MinPeakDistance',2, 'SortStr','descend');
%     plot(time,cap,locs,pks,'o')

    %z=max(cap{:});
    %invertedcap= cellfun(@(x) z - x, cap, 'Un', false);
    %[val,vloc,vwidths,vproms] = findpeaks(invertedcap, time, 'MinPeakDistance',2,'SortStr','descend');
    % plot(time,cap,vlocs,val,'o')
%     if i==1
%         linetype = 'r--';
%     else
%         linetype = 'b:';
%     end
%     linetype = {'g--','r--'}; %'g--','r--','c--','m--'};
%     for n=1:length(Ttime)
% %         line([Ttime(n) Ttime(n)],ylim,linetype);
% %         text(Ttime(n),min(ylim),Tlabel(n),'Rotation',90);
%           vline(Ttime(n),linetype{2}, Tlabel(n));
%            
%     end
% 
% %     text(locs+10,pks,num2str((1:numel(pks))'))
% title('Male, 37yrs')
% xlabel('time(1/30 s)');
% ylabel('Capacitance (pF)');
% 
% for n = 1:4
%     plot(NaN,NaN,linetype{n});
% end
% 
% legend('Xiphoid', 'Move'); %'Start/End', 'Move', 'Deep Breath', 'Misc');


% xml = xml2struct(char(fullfile(Folder,wFolder,'Spiro_XML_export.xml')));

% XMLfile = char(fullfile(Folder, wFolder, 'Spiro_XML_export.xml'));
% 
% xmlDoc = xmlread(XMLfile);
% 
% docNode = xmlDoc.getDocumentElement;
% 
% entries = docNode.getChildNodes;
% 
% info = entries.item(0).getChildNodes;
% 
% nodes = info.getElementsByTagName('ChannelFlow');


%% CrossCorrelation

import javax.xml.xpath.*

Folder = 'Z:\GitRepositories\stretch-sense\Data';
wFolder = '\Spirometry';


XMLfile = char(fullfile(Folder, wFolder, 'Spiro_XML_export.xml'));

xmlDoc = xmlread(XMLfile);

factory = XPathFactory.newInstance;
xpath = factory.newXPath;

FlowData = xpath.compile('//ChannelVolume/SamplingValues');

FlowNodes = FlowData.evaluate(xmlDoc, XPathConstants.NODESET);

%extracting Spirometer Traces from the XML
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

node1 = FlowNodes.item(FlowNodes.getLength-3);
C1 = strsplit(char(node1.getFirstChild.getNodeValue));
C1 = str2double(C1);
C1 = transpose(C1);
C1time = linspace(0,length(C1)/100,length(C1));
C1time = transpose(C1time);
tsC1=timeseries(C1, C1time);



Folder = 'Z:\GitRepositories\stretch-sense\Data\';
%MACfilename = '/Users/justinschaffner/Desktop/GitRepositories/stretch-sense/Data/SenseAppData/CAP_2018-03-07542162368_U_R_SIDE.csv';
Filename = 'SenseAppData\Xiphoid\CAP\CAPJustin_T2_2018-04-27546540928.csv'
Gfilename = 'SenseAppData\Xiphoid\GT\GTJustin_T2_2018-04-27546540928.csv'

%GTfile='/Users/justinschaffner/Desktop/GitRepositories/stretch-sense/Data/SenseAppData/GT_2018-03-07542161920_U_F_CENTER.csv';
%SGTfile='/Users/justinschaffner/Desktop/GitRepositories/stretch-sense/Data/SenseAppData/GT_2018-03-07542162176_U_R_CORNER.csv';
%capfiles = ["Xiphoid\CAP_2018-03-17_X_B_CENTER_EILEEN_31_5inch.csv","4thIntercostal\CAP_2018-03-17_4th_B_Center_EILEEN_34_5inch.csv"];
%GTfiles = ["Xiphoid\GT_2018-03-17_X_B_CENTER_EILEEN_31_5inch.csv","4thIntercostal\GT_2018-03-17_4th_B_Center_EILEEN_34_5inch.csv"];


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

figure; %hold on;
A1 = getsampleusingtime(ts,Ttime(2)-1.45,Ttime(3)+1.7);

A1p = A1.Data;
C1p = tsC1.Data;
length = min([numel(A1p) numel(C1p)]);
A1p = A1p(1:length);
C1p = C1p(1:length);

A2=getsampleusingtime(ts,Ttime(4)-1.2,Ttime(5)+1.5);
A2p = A2.Data;
C2p = tsC2.Data;
length = min([numel(A2p) numel(C2p)]);
A2p = A2p(1:length);
C2p = C2p(1:length);

A3 = getsampleusingtime(ts,Ttime(6)-1,Ttime(7)+1.5);
A3p = A3.Data;
C3p = tsC3.Data;
length = min([numel(A3p) numel(C3p)]);
A3p = A3p(1:length);
C3p = C3p(1:length);

A = [A1p; A2p; A3p];
C = [C1p; C2p; C3p];

%[sA1,stsC1] = synchronize(A1,tsC1,'Uniform','Interval',Fs);
crosscorr(A,C);
 figure;
 scatter(A, C);
%plot(A1);
figure;hold on;

plot(A1);
plot(A2);
plot(A3);
figure;hold on;

plot(tsC1);
plot(tsC2);
plot(tsC3);
    

%     %  Plot data
%     plot(ts)
%     
%     % Find peaks in data
% %     [pks,locs,widths,proms] = findpeaks(cap, time,'MinPeakDistance',2, 'SortStr','descend');
% %     plot(time,cap,locs,pks,'o')
% 
%     %z=max(cap{:});
%     %invertedcap= cellfun(@(x) z - x, cap, 'Un', false);
%     %[val,vloc,vwidths,vproms] = findpeaks(invertedcap, time, 'MinPeakDistance',2,'SortStr','descend');
%     % plot(time,cap,vlocs,val,'o')
% %     if i==1
% %         linetype = 'r--';
% %     else
% %         linetype = 'b:';
% %     end
%     linetype = {'g--','r--'}; %'g--','r--','c--','m--'};
%     for n=1:length(Ttime)
% %         line([Ttime(n) Ttime(n)],ylim,linetype);
% %         text(Ttime(n),min(ylim),Tlabel(n),'Rotation',90);
%           vline(Ttime(n),linetype{2}, Tlabel(n));
%            
%     end
% 
% %     text(locs+10,pks,num2str((1:numel(pks))'))
% title('Male, 37yrs')
% xlabel('time(1/30 s)');
% ylabel('Capacitance (pF)');
% % 
% % for n = 1:4
% %     plot(NaN,NaN,linetype{n});
% % end
% % 
% % legend('Xiphoid', 'Move'); %'Start/End', 'Move', 'Deep Breath', 'Misc');
% 
% 
% % xml = xml2struct(char(fullfile(Folder,wFolder,'Spiro_XML_export.xml')));
% 
% % XMLfile = char(fullfile(Folder, wFolder, 'Spiro_XML_export.xml'));
% % 
% % xmlDoc = xmlread(XMLfile);
% % 
% % docNode = xmlDoc.getDocumentElement;
% % 
% % entries = docNode.getChildNodes;
% % 
% % info = entries.item(0).getChildNodes;
% % 
% % nodes = info.getElementsByTagName('ChannelFlow');