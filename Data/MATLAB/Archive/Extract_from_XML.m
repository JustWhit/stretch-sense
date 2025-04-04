import javax.xml.xpath.*

Folder = 'Z:\GitRepositories\stretch-sense\Data';
sFolder = '\RawSpirometry';
dFolder = '\SensorData';

XMLfile = char(fullfile(Folder, sFolder, 'Spiro_7_25_18.xml'));
% % Assign test NAME HERE
TestName = 'Spiro_7_25_18_JUSTIN_SVC_TEST8';
xmlDoc = xmlread(XMLfile);

factory = XPathFactory.newInstance;
xpath = factory.newXPath;

FlowData = xpath.compile('//ChannelVolume/SamplingValues');

FlowNodes = FlowData.evaluate(xmlDoc, XPathConstants.NODESET);

for i = 1:5
csvwrite(char(fullfile(Folder, strcat(TestName,'_T', num2str((5-i)+1),'.csv'))),getSpiro(FlowNodes,i));
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