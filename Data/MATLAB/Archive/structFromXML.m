import javax.xml.xpath.*

Folder = 'Z:\GitRepositories\stretch-sense\Data';
wFolder = '\Spirometry';

XMLfile = char(fullfile(Folder, wFolder, 'Spiro_XML_export.xml'));

xmlDoc = xmlread(XMLfile);

factory = XPathFactory.newInstance;
xpath = factory.newXPath;

FlowData = xpath.compile('//Patient');

FlowNodes = FlowData.evaluate(xmlDoc, XPathConstants.NODESET);

for i=1:length(FlowNodes)
   Patients(i)=xml2struct(FlowNodes(i)); 
end




% for i=1:length(Patients)
%    figure;
%    plot(Patients(i).Tests
% end