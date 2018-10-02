
TFolder = 'Y:\Justin\GitRepositories\stretch-sense\Data\TestData';
FileListT = dir(char(fullfile(TFolder,'*.csv')));

for n=1:length(FileListT)
    Table = readtable(char(fullfile(TFolder,FileListT(n).name))); 
    disp(FileListT(n).name);
    Time = [Table{:,1}];
    Sensors = [Table{:,2:4}];
    Time = alterDuplicates(Time);


    SensorTS = timeseries(Sensors,Time);
%     SensorTS.Data(:,1) = sgolayfilt(SensorTS.Data(:,1),3,25);
    SensorTS.Data(:,1) = denoise(SensorTS.Data(:,1));
    CapF = SensorTS.Data(:,1);
% % \\\\\\\\\ CHANGE THE VALUE HERE TO SWITCH FROM VOLUME TO FLOW\\\\\\
    SpiroF = SensorTS.Data(:,3);
    VolF = SensorTS.Data(:,2);
    CapF = truncateTrace(CapF);
% % % Add Time Variable for CapF calculation\\\\\\\\\\\\\\\\\\\\\\\\
%     TimeF = SensorTS.Time;
%     TimeF = TimeF(1:numel(CapF));

%     figure;hold on; plot(CapF);plot(SpiroF.^2);
    [CapF, durations] = offsetData(CapF,SensorTS.Time);
    
    Length = numel(CapF);
    SpiroF = SpiroF(2:Length+1);
figure;hold on;plot(CapF);plot(SpiroF);
    VolF = VolF(2:Length+1);
    TimeF = SensorTS.Time;
    TimeF = TimeF(2:Length+1);
    SpiroTS = timeseries(SpiroF,TimeF);
    VOLTS = timeseries(VolF,TimeF);
% %\\\\\\\\\\\ COMMENT OUT NEXT LINE IF USING FLOW\\\\\\\\\\     
%     SpiroF = offsetData(SpiroF);

% % \\\Fourier Models

f = fourier8Model2(CapF);
f = sgolayfilt(f,3,25);
fTS = timeseries(f,TimeF);
x0 = xintercepts(f,TimeF);
estVT = getVT(fTS,durations,x0);
actVT = getVT(SpiroTS,durations,x0);
VOLVT = getVTfromVOL(VOLTS,x0);
RR = getRR(x0);
% MAC = smooth(f,75);


plotActualVSPred(SpiroTS,fTS, ':fouri:', n,x0, estVT,actVT,VOLVT,RR);


% % % Gaussian Model
% gauss = gauss8Model(CapF);
% gaussTS = timeseries(gauss, TimeF);
% x0 = xintercepts(gauss,TimeF);
% estVT = getVT(gaussTS,durations,x0);
% actVT = getVT(SpiroTS,durations,x0);
% VOLVT = getVTfromVOL(VOLTS,x0);
% RR = getRR(x0);
% plotActualVSPred(SpiroTS,gaussTS, ':gauss:',n,x0, estVT,actVT,VOLVT,RR);


% % % Sine Model
% sin = sin8Model(CapF);
% sinTS = timeseries(sin, TimeF);
% x0 = xintercepts(sin,TimeF);
% estVT = getVT(sinTS,durations,x0);
% actVT = getVT(SpiroTS,durations,x0);
% VOLVT = getVTfromVOL(VOLTS,x0);
% RR = getRR(x0);
% plotActualVSPred(SpiroTS,sinTS, ':Sin:',n,x0, estVT,actVT,VOLVT,RR);


end


% % all data, (neg neg, pos pos) included
function test = sin8Model(test)
    a1 =      0.3345; %  (-304.6, 305.3)
    b1 =       4.103; %  (-931.4, 939.6)
    c1 =       -1.46; %  (-2007, 2004)
    a2 =      0.1276; %  (-12.93, 13.18)
    b2 =       10.09; %  (-2480, 2500)
    c2 =      -2.455; %  (-224.2, 219.3)
    a3 =      0.1066; %  (-84.98, 85.19)
    b3 =       12.97; %  (-1131, 1156)
    c3 =      0.5137; %  (-50.36, 51.39)
    a4 =       1.528; %  (-116, 119.1)
    b4 =       1.883; %  (-957.9, 961.6)
    c4 =       3.098; %  (-39.04, 45.23)
    a5 =      0.3918; %  (-228.1, 228.9)
    b5 =       6.288; %  (-984.8, 997.4)
    c5 =       1.845; %  (-602.4, 606.1)
    a6 =      0.1257; %  (-6.32, 6.571)
    b6 =       16.62; %  (-31.25, 64.5)
    c6 =      -2.752; %  (-4.075, -1.43)
    a7 =    0.004874; %  (-0.001799, 0.01155)
    b7 =         795; %  (787.3, 802.7)
    c7 =      -1.864; %  (-3.22, -0.5076)
    a8 =    0.002297; %  (-0.004351, 0.008944)
    b8 =       960.6; %  (944.3, 976.9)
    c8 =     -0.6297; %  (-3.52, 2.26)

    for n=1:numel(test)
        x = test(n);
        test(n) = a1*sin(b1*x+c1) + a2*sin(b2*x+c2) + a3*sin(b3*x+c3) + a4*sin(b4*x+c4) + a5*sin(b5*x+c5) + a6*sin(b6*x+c6) + a7*sin(b7*x+c7) + a8*sin(b8*x+c8);
    end
end


% % all data, (neg neg, pos pos) included
function test = gauss8ModelMirrored(test)
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
        if x>0
            test(n) = a1*exp(-((x-b1)/c1)^2) + a2*exp(-((x-b2)/c2)^2) + a3*exp(-((x-b3)/c3)^2) + a4*exp(-((x-b4)/c4)^2) + a5*exp(-((x-b5)/c5)^2) + a6*exp(-((x-b6)/c6)^2) + a7*exp(-((x-b7)/c7)^2) + a8*exp(-((x-b8)/c8)^2);
        else
            x = -x;
            test(n) = -(a1*exp(-((x-b1)/c1)^2) + a2*exp(-((x-b2)/c2)^2) + a3*exp(-((x-b3)/c3)^2) + a4*exp(-((x-b4)/c4)^2) + a5*exp(-((x-b5)/c5)^2) + a6*exp(-((x-b6)/c6)^2) + a7*exp(-((x-b7)/c7)^2) + a8*exp(-((x-b8)/c8)^2));
        end
    end


end

% % all data, (neg neg, pos pos) included
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

% % Outliers (neg neg, pos pos) excluded
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

% % all data, (neg neg, pos pos) included
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
%     assignin('base','TraceIN',Trace);
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
%     assignin('base','TraceOUT',Trace);

end

% % shifts capacitance to baseline
function [CapF, Durations] = offsetData(Data, Time)
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
    CapF = diff(Data);
    Durations = diff(Time(1:numel(Data),:));
    
% % %  Difference over Time
%     assignin('base','DiffData',diff(Data));
%     assignin('base','DiffTime',diff(Time));
%     a = diff(Data)./diff(Time);
end

function plotActualVSPred(actual,pred, Model, iteration,x0,estVT,actVT,VOLVT,RR)
    y0 = zeros(numel(x0));
    RMSE = (sqrt(mean((actual.Data-pred.Data).^2)))/(max(pred.Data)-min(pred.Data));
    figure, 
    hold on, 
    plot(pred,'o'),
    plot(actual,'o'), 
    plot(x0,y0,'r*','markers',12);
%     plot(MAC,'r--');
    title(['Flow' Model num2str(iteration) ':rRMSE:' num2str(RMSE) ':estVT:' num2str(estVT) ':actVT:' num2str(actVT) ':VOLVT:' num2str(VOLVT) ':RR:' num2str(RR)]), 
    legend({'Predicted', 'Actual','X-Intercepts'}),
    xlabel('Samples'), 
    ylabel('Flow(L/s)');

end

function a = xintercepts(Pred,TimeF)
    a = [];
    for n=1:numel(Pred)-1
       
       if (Pred(n)>0 && Pred(n+1) < 0) || (Pred(n)<0 && Pred(n+1)>0)
           m = (Pred(n)-Pred(n+1))/(TimeF(n)-TimeF(n+1));
           b = Pred(n)-(TimeF(n)*m);
           x = -b/m;
           a = [a;x];
       elseif Pred(n)==0
           a = [a;TimeF(n)];    
       end
    end
end

% % returns the Estimated Tidal Volume of a trace 
function a  = getVT(PredTS, Duration, x0)
   data = PredTS.Data;
   assignin('base','data',data);
   PredTS.Data = data.*Duration;
    Len = numel(x0);
%     rem = mod(Len,3);
%     Len = Len-rem;
    if Len < 7
        a = 0;
        return;
    end
    VT = [];
    for n=1:2:5
% % %         Volume using CumSum
%         temp = getsampleusingtime(PredTS,x0(n),x0(n+1));
%         inhale = (cumsum(abs(temp.Data(:))));
%         temp = getsampleusingtime(PredTS, x0(n+1),x0(n+2));
%         exhale = (cumsum(abs(temp.Data(:))));

% % %    Volume using cumtrapz
%         temp = getsampleusingtime(PredTS,x0(n),x0(n+1));
%         inhale = cumtrapz(abs(temp.Data(:)));
%         temp = getsampleusingtime(PredTS, x0(n+1),x0(n+2));
%         exhale = cumtrapz(abs(temp.Data(:)));

% %    Volume using Trapz
        temp = getsampleusingtime(PredTS,x0(n),x0(n+1));
        inhale = trapz(abs(temp.Data(:)));
        temp = getsampleusingtime(PredTS, x0(n+1),x0(n+2));
        exhale = trapz(abs(temp.Data(:)));

        VT = [VT;inhale;exhale];
    end
    assignin('base','VTarray',VT);
    a = mean(VT);
end

% % get VT from Volume data
function a = getVTfromVOL(VOLTS, x0)
  
     Len = numel(x0);
%     rem = mod(Len,3);
%     Len = Len-rem;
    if Len < 7
        a = 0;
        return;
    end
    VT = [];
    for n=1:2:5
% % %         Volume using cumsum
%         temp = getsampleusingtime(VOLTS,x0(n),x0(n+1));
%         B1 = (cumsum(abs(diff(temp.Data(:)))));
%         temp = getsampleusingtime(VOLTS, x0(n+1),x0(n+2));
%         B2 = (cumsum(abs(diff(temp.Data(:)))));
        
% %         Volume using trapz
        temp = getsampleusingtime(VOLTS,x0(n),x0(n+1));
        B1 = (trapz(abs(diff(temp.Data(:)))));
        temp = getsampleusingtime(VOLTS, x0(n+1),x0(n+2));
        B2 = (trapz(abs(diff(temp.Data(:)))));
        
        VT = [VT;B1;B2];
        
    end
    a = mean(VT);
end

% % alters duplicate time stamps
function time = alterDuplicates(time)
    [uniqueA i j] = unique(time,'first');
    indexToDupes = find(not(ismember(1:numel(time),i)));
    time(indexToDupes) = time(indexToDupes)+0.00001;
end

% % Returns Respiratory Rate from x-intercepts
function a = getRR(x0)
    Len = numel(x0);
    if Len < 7
        a = 0;
        return;
    end
    duration = x0(7)-x0(1);
    a = (3/duration)*60;
end

% % Denoising the data
function a = denoise(x)
    level = 5;
    wname = 'sym4';
    tptr  = 'sqtwolog';
    sorh  = 's';
    npc_app = 1;
    npc_fin = 1;
%     npc_app = 'kais';
%     npc_fin = 'kais';
    a = wmulden(x, level, wname, npc_app, npc_fin, tptr, sorh);
    
end