


TFolder = 'Y:\Justin\GitRepositories\stretch-sense\Data\ModelData';
FileListT = dir(char(fullfile(TFolder,'*.csv')));
% % Xiphoid


% % Build Data for Model
Features = [];
% figure; hold on;
for n=1:length(FileListT)
   Table = readtable(char(fullfile(TFolder,FileListT(n).name))); 
   disp(FileListT(n).name);
   Time = [Table{:,1}];
   Sensors = [Table{:,2:4}];
   
%    Sensors(:,1) = offsetCapacitance(Sensors(:,1));
   
   SensorTS = timeseries(Sensors,Time);
   SensorTS.Data(:,1) = sgolayfilt(SensorTS.Data(:,1),3,25);

  
% % % \\\\\\\Make changes here to adjust the Features Used by the Model\\\
   
%    [C,S] = getPeaksVals(SensorTS);
% 
% % % Difference between Peaks and Valleys (inhale) minus deep breath
%     CapF = [abs(C(:,1) + C(:,3)) abs(C(:,2) - C(:,4))];
%     SpiroF = abs(S(:,1) + S(:,3));
%     Features = [Features; [CapF SpiroF]];

% % % Values between Peaks and Valleys (inhale)
% figure; hold on;
% plot(SensorTS.Data(:,1));
%     SensorTS.Data(:,1) = sgolayfilt(SensorTS.Data(:,1),3,51);
% plot(SensorTS.Data(:,1));
%     SensorTS.Data(:,1) = offsetData(SensorTS.Data(:,1));
%     SensorTS.Data(:,2) = offsetData(SensorTS.Data(:,2));
% plot(SensorTS.Data(:,2));
%     for i=1:numel(C(:,1))
%         
%         [CapF, SpiroF] = resizeDataSet(SensorTS.Data(C(i,4):C(i,2),1), SensorTS.Data(S(i,4):S(i,2),2) );
%         Features = [Features; [CapF SpiroF]];
%     end

% % All Values\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    
    CapF = SensorTS.Data(:,1);
% % \\\\\\\\\ CHANGE THE VALUE HERE TO SWITCH FROM VOLUME TO FLOW\\\\\\
    SpiroF = SensorTS.Data(:,3);
    
%     CapF = truncateTrace(CapF);
% % % Add Time Variable for CapF calculation\\\\\\\\\\\\\\\\\\\\\\\\
%     TimeF = SensorTS.Time;
%     TimeF = TimeF(1:numel(CapF));
   
%     figure;hold on; plot(CapF);plot(SpiroF.^2);
    CapF = offsetData(CapF);
%     SpiroF=offsetData(SpiroF);
    Length = numel(CapF);
    SpiroF = SpiroF(1:Length);
% %\\\\\\\\\\\ COMMENT OUT NEXT LINE IF USING FLOW\\\\\\\\\\     
%     SpiroF = offsetData(SpiroF);

%     CapF = sgolayfilt(CapF,3,11);
    Features = [Features; [CapF SpiroF]];
    
end
% Features = shuffleRow(Features);
% c1 = Features(:,1)< 0;
% c2 = Features(:,2)<0;
% Cneg = all(c1 & c2,2);
% Features(Cneg,:) = [];
% c3 = Features(:,1)>0;
% c4 = Features(:,2)>0;
% Cpos = all(c3 & c4,2);
% Features(Cpos,:)=[];
% c5 = Features(:,1)==0;
% c6 = Features(:,2)==0;
% Czero = all(c5 & c6);
% Features(Czero,:)=[];



TRsetSize = floor(.8 * length(Features(:,1)));
inputs = Features(1:TRsetSize,1);
targets = Features(1:TRsetSize,2);
test = Features(TRsetSize+1:end,1);
actual = Features(TRsetSize+1:end,2);

% [Param1,Param2] = ransac(transpose(Features),2,100,20,0.1);

% % Train\Test Model
% mdl = fitlm(inputs,targets)
f = fit(inputs,targets,'fourier8')
figure;plot(f,inputs,targets);
% [pred, predCI] = predict(mdl,test);
predf = fourier8Model2(test);
% pred = sgolayfilt(pred,3,51);
plotTrainingFit(inputs,targets,1);
% zerofree = find(pred);
plotActualVSPred(actual,predf);
%  diagnosticPlots(f);
% residualPlots(mdl);
% errorHistogram(pred,actual);
% plotScatter(inputs,targets);

function test = gauss8Model(test)
    a1 =      -0.405;
    b1 =      -0.653;
    c1 =   0.0004993;
    a2 =       2.838;
    b2 =     -0.8753;
    c2 =      0.2795;
    a3 =      -4.827;
    b3 =       7.488;
    c3 =       5.911;
    a4 =        1.12;
    b4 =     -0.4689;
    c4 =      0.1827;
    a5 =     0.08044;
    b5 =      -0.535;
    c5 =       0.188;
    a6 =      0.1087;
    b6 =     -0.4317;
    c6 =     0.00652;
    a7 =           0;
    b7 =      -1.802;
    c7 =      0.1367;
    a8 =        1.49; %  (0.3472, 2.632)
    b8 =     -0.1571; %  (-0.3009, -0.01335)
    c8 =      0.2414; %  (0.1521, 0.3307)
    for n=1:numel(test)
        x = test(n);
        test(n) = a1*exp(-((x-b1)/c1)^2) + a2*exp(-((x-b2)/c2)^2) + a3*exp(-((x-b3)/c3)^2) + a4*exp(-((x-b4)/c4)^2) + a5*exp(-((x-b5)/c5)^2) + a6*exp(-((x-b6)/c6)^2) + a7*exp(-((x-b7)/c7)^2) + a8*exp(-((x-b8)/c8)^2);
    end
end

function a = fourier8Model2(test)
    a0 =     0.06146; %  (0.01234, 0.1106)
       a1 =    -0.08757; %  (-0.1757, 0.0005277)
       b1 =     -0.8026 ; % (-0.831, -0.7742)
       a2 =     0.04069; %  (-0.0405, 0.1219)
       b2 =      0.1523; %  (0.1071, 0.1974)
       a3 =    -0.06337; %  (-0.1362, 0.009449)
       b3 =     -0.2392; %  (-0.2968, -0.1815)
       a4 =      0.0322; %  (-0.03075, 0.09515)
       b4 =     0.07387; %  (-0.00219, 0.1499)
       a5 =    -0.01994; %  (-0.08752, 0.04763)
       b5 =     -0.1429; %  (-0.2528, -0.0329)
       a6 =     0.04706; %  (-0.01123, 0.1053)
       b6 =      0.0124; %  (-0.1328, 0.1576)
       a7 =    -0.02861; %  (-0.06497, 0.007753)
       b7 =   -0.001317; %  (-0.1219, 0.1193)
       a8 =     0.01578; %  (-0.000806, 0.03237)
       b8 =    -0.05113; %  (-0.1101, 0.007809)
       w =       6.064; %  (5.534, 6.594)
       
       for n =1: numel(test)
           x = test(n);
          test(n) = a0 + a1*cos(x*w) + b1*sin(x*w) + a2*cos(2*x*w) + b2*sin(2*x*w) + a3*cos(3*x*w) + b3*sin(3*x*w) + a4*cos(4*x*w) + b4*sin(4*x*w) + a5*cos(5*x*w) + b5*sin(5*x*w) + a6*cos(6*x*w) + b6*sin(6*x*w) + a7*cos(7*x*w) + b7*sin(7*x*w) + a8*cos(8*x*w) + b8*sin(8*x*w); 
       end
       a = test;
       
end

function a = fourier8Model(test)
    a0 = 0.08946; %  (-0.2539, 0.4328)
       a1 =  -0.1482; %  (-0.8027, 0.5063)
       b1 =  -0.7914; %  (-0.8315, -0.7514)
       a2 =   0.112; %  (-0.4755, 0.6996)
       b2 =  0.1499; %  (0.03283, 0.2669)
       a3 =  -0.106; %  (-0.5897, 0.3777)
       b3 = -0.2166; %  (-0.4001, -0.03303)
       a4 =  0.0489; %  (-0.3478, 0.4456)
       b4 =   0.0359; %  (-0.2077, 0.2795)
       a5 = -0.03665; %  (-0.4009, 0.3276)
       b5 = -0.07187; %  (-0.3247, 0.1809)
       a6 =  0.05729; %  (-0.2114, 0.326)
       b6 = -0.02341; %  (-0.2506, 0.2037)
       a7 =  -0.03193; %  (-0.1802, 0.1163)
       b7 = -0.003389; %  (-0.1353, 0.1285)
       a8 = 0.02355; %  (-0.03114, 0.07824)
       b8 =    -0.02062; %  (-0.09149, 0.05025)
       w =       5.667; %  (4.238, 7.095)
       
       for n =1: numel(test)
           x = test(n);
          test(n) = a0 + a1*cos(x*w) + b1*sin(x*w) + a2*cos(2*x*w) + b2*sin(2*x*w) + a3*cos(3*x*w) + b3*sin(3*x*w) + a4*cos(4*x*w) + b4*sin(4*x*w) + a5*cos(5*x*w) + b5*sin(5*x*w) + a6*cos(6*x*w) + b6*sin(6*x*w) + a7*cos(7*x*w) + b7*sin(7*x*w) + a8*cos(8*x*w) + b8*sin(8*x*w); 
       end
       a = test;
end


% % Removes Portions of each trace
function Trace = truncateTrace(Trace)
% % % Removes Full inhale/exhale
% figure;hold on;
    [pks,locs,~,~] = findpeaks(Trace, 'MinPeakProminence',0.75,'MinPeakDistance',100);
    [vals,vlocs,~,~] = findpeaks(-Trace, 'MinPeakProminence',0.75,'MinPeakDistance',100);
    assignin('base','TraceIN',Trace);
% plot(Trace);
    maxi = max(pks);
    index = find(pks == maxi);
    locusMax = locs(index);
% plot(locusMax,maxi,'r*');
    fIndex = find(vlocs < locusMax, 1, 'last');
    final = vlocs(fIndex);
% plot(final,-vals(fIndex),'b*');
    Trace = Trace(1:final);
% plot(Trace);
    assignin('base','TraceOUT',Trace);

end

% % RANSAC function taken from RANSAC wikipedia page
function [bestParameter1,bestParameter2] = ransac(data,num,iter,threshDist,inlierRatio)
 % data: a 2xn dataset with #n data points
 % num: the minimum number of points. For line fitting problem, num=2
 % iter: the number of iterations
 % threshDist: the threshold of the distances between points and the fitting line
 % inlierRatio: the threshold of the number of inliers 
 
 % % Plot the data points
 figure;plot(data(1,:),data(2,:),'o');hold on;
 number = size(data,2); % Total number of points
 bestInNum = 0; % Best fitting line with largest number of inliers
 bestParameter1=0;bestParameter2=0; % parameters for best fitting line
 for i=1:iter
 % % Randomly select 2 points
     idx = randperm(number,num); sample = data(:,idx);   
 % % Compute the distances between all points with the fitting line 
     kLine = sample(:,2)-sample(:,1);% two points relative distance
     kLineNorm = kLine/norm(kLine);
     normVector = [-kLineNorm(2),kLineNorm(1)];%Ax+By+C=0 A=-kLineNorm(2),B=kLineNorm(1)
     distance = normVector*(data - repmat(sample(:,1),1,number));
 % % Compute the inliers with distances smaller than the threshold
     inlierIdx = find(abs(distance)<=threshDist);
     inlierNum = length(inlierIdx);
 % % Update the number of inliers and fitting model if better model is found     
     if inlierNum>=round(inlierRatio*number) && inlierNum>bestInNum
         bestInNum = inlierNum;
         parameter1 = (sample(2,2)-sample(2,1))/(sample(1,2)-sample(1,1));
         parameter2 = sample(2,1)-parameter1*sample(1,1);
         bestParameter1=parameter1; bestParameter2=parameter2;
     end
 end
 
 % % Plot the best fitting line
 xAxis = -number/2:number/2; 
 yAxis = bestParameter1*xAxis + bestParameter2;
 plot(xAxis,yAxis,'r-','LineWidth',2);
end


% % shifts capacitance to baseline
function a = offsetData(Data, Time)
% % %     Offset by Median Valley
%     [vals,locs,~,~] = findpeaks(-Data, 'MinPeakProminence',0.75,'MinPeakDistance',100);
%     base = median(-vals);
%     a=Data - (base - 20);

% % %     Offset by Max
%     offset = max(Data) - 100;
%     a = Data - offset;

% % %     Offset by Minimum
%     offset = min(Data);
%     a = Data - offset;

% % % Normalize by feature scaling
%     minimum = min(Data);
%     maximum = max(Data);
%     Data = (Data - minimum);
%     a = Data/(maximum - minimum);
    
% % % Normalize by Standardization
%     Mean = mean(Data);
%     Std = std(Data);
%     Data = Data - Mean;
%     a = Data / Std;

% % % % offset by first value
%     offset = Data(1);
%     Data = Data - offset;
%     a = Data;

% %  differences
    a = diff(Data);
    

end

% % Returns the two sets of data points, stretched to fit the longest set.
function [a, b] = resizeDataSet( setA, setB)
   
    assignin('base','setA',setA);
    assignin('base','setB',setB);
    
    if(numel(setA) > numel(setB))
        a = setA;
        n = numel(setB);
        b = transpose(interp1(1:n, setB, linspace(1,n,numel(setA)), 'nearest'));
    else
        b = setB;
        n = numel(setA);
        a = transpose(interp1(1:n, setA, linspace(1,n,numel(setB)),'nearest'));
    end
        
end

% % Returns teh percent difference between two values
function a = pDiff(v1, v2)
    a = (abs(v1-v2)/((v1+v2)/2))*100;
end

% % Shuffles the Data
function ret = shuffleRow(mat)

    [r, ~] = size(mat);
    shuffledRow = randperm(r);
    ret = mat(shuffledRow, :);
end

% % Returns Peak Valley pairs for Cap and Spiro traces
function [a, b] = getPeaksVals(SensorsTS)
    
    % % Find peaks in data
    [Cpks,Clocs,~,~] = findpeaks(SensorsTS.Data(:,1), 'MinPeakProminence',0.1); %,'MinPeakDistance',100);
    [Cvals,Cvlocs,~,~] = findpeaks(-SensorsTS.Data(:,1), 'MinPeakProminence',0.1); %,'MinPeakDistance',100);
    
    [Spks,Slocs,~,~] = findpeaks(SensorsTS.Data(:,2), 'MinPeakProminence',0.1);
    [Svals,Svlocs,~,~] = findpeaks(-SensorsTS.Data(:,2), 'MinPeakProminence',0.1);
    
% plotPeaksVals(SensorsTS,Cpks,Clocs,Cvals,Cvlocs,Spks,Slocs,Svals,Svlocs);
    length = min([numel(Spks) numel(Svals)])-1; % % subtract 1 to remove the last big Inhale
    Cpks = Cpks(1:length);
    Clocs = Clocs(1:length);
    Cvals = Cvals(1:length);
    Cvlocs = Cvlocs(1:length);
    Spks = Spks(1:length);
    Slocs = Slocs(1:length);
    Svals = Svals(1:length);
    Svlocs = Svlocs(1:length);
%     plotPeaksVals(SensorsTS,Cpks,Clocs,Cvals,Cvlocs,Spks,Slocs,Svals,Svlocs);
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
    title('Linear Fit of Training Data');
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
    title('Predicted and Actual Flow'), 
    legend({'Predicted', 'Actual'}),
    xlabel('Samples'), 
    ylabel('Flow(L/s)');

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
    plot(SensorTS.Data(:,1));
    plot(Clocs,Cpks,'o'); 
    ylabel('Capacitance(pF)');
    xlabel('Time(s)'); 
    grid on; pFREQ = numel(Cpks)/max(Clocs); title(['Peak Detection Freq Estimate = ', num2str(pFREQ), ' Hz; RR est = ', num2str(pFREQ*60),'(1/min)']);
    plot(Cvlocs,-Cvals,'r*');
    plot(SensorTS.Data(:,2));
    plot(Slocs,Spks,'o');
    plot(Svlocs,-Svals,'r*');
end

function plotScatter(Capacitance, Spirometer)
    figure;
    scatter(Capacitance,Spirometer);
    title('Scatter Plot of Training Data');
    xlabel('Capacitance(pF)');
    ylabel('Volume(L)');
end