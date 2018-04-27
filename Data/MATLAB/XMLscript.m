
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
    ts=timeseries(cap,time);
    
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