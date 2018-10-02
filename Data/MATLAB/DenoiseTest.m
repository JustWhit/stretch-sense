

Folder = 'Y:\Justin\GitRepositories\stretch-sense\Data\ModelData';
FileList = dir(char(fullfile(Folder,'*.csv')));

Traces = {};

for n=1:length(FileList)
    figure;hold on;
    Table = readtable(char(fullfile(Folder,FileList(n).name))); 
    disp(FileList(n).name);
    Time = [Table{:,1}];
    Sensors = [Table{:,2:4}];
    Time = alterDuplicates(Time);
    
    
    SensorTS = timeseries(Sensors,Time);
    [CapF, ~] = offsetData(Sensors(:,1),Time);
    
%     CapF = truncateTrace(Sensors(:,1));
%     [c,l] = wavedec(CapF,3,'db2');
%     approx = appcoef(c,l,'db2');
%     [cd1,cd2,cd3] = detcoef(c,l,[1 2 3]);
%     subplot(5,1,1)
%     plot(approx)
%     title('Approximation Coefficients')
%     subplot(5,1,2)
%     plot(cd3)
%     title('Level 3 Detail Coefficients')
%     subplot(5,1,3)
%     plot(cd2)
%     title('Level 2 Detail Coefficients')
%     subplot(5,1,4)
%     plot(cd1)
%     title('Level 1 Detail Coefficients')
%     subplot(5,1,5)
%     plot(CapF)
%     title('Capacitance(pF)');
    plot(CapF);
%     SensorTS.Data(:,1) = sgolayfilt(SensorTS.Data(:,1),3,25);
    CapF = denoise(Sensors(:,1));
    [CapF, durations] = offsetData(CapF,Time);
%     CapF = denoise(CapF);
    plot(CapF);
%     CapF = SensorTS.Data(:,1);
% % \\\\\\\\\ CHANGE THE VALUE HERE TO SWITCH FROM VOLUME TO FLOW\\\\\\
    FlowF = SensorTS.Data(:,3);
    VolF = SensorTS.Data(:,2);
    CapF = truncateTrace(CapF);
% % % Add Time Variable for CapF calculation\\\\\\\\\\\\\\\\\\\\\\\\
%     TimeF = SensorTS.Time;
%     TimeF = TimeF(1:numel(CapF));

%     figure;hold on; plot(CapF);plot(SpiroF.^2);
%     [CapF, durations] = offsetData(CapF,SensorTS.Time);

    Length = numel(CapF);
%     FlowF = FlowF(2:Length);
% figure;hold on;plot(CapF);plot(FlowF);
%     VolF = VolF(2:Length);
    TimeF = SensorTS.Time;
%     TimeF = TimeF(2:Length);
    plot(FlowF);
    title(FileList(n).name, 'Interpreter', 'none');
    xlabel('Time(s)');
    ylabel('Value');
    legend({'\Delta Capacitance(pF)','Denoised \Delta Capacitance(pF)','Flow(L/s)'});
    
    Traces(n,:) = {CapF,FlowF,VolF,TimeF,durations};
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

function plotActualVSPred(actual,pred, Name,x0,estVT,actVT,RR)
    y0 = zeros(numel(x0));
    RMSE = (sqrt(mean((actual.Data-pred.Data).^2)))/(max(pred.Data)-min(pred.Data));
    figure, 
    hold on, 
    plot(pred,'o'),
    plot(actual,'o'), 
    plot(x0,y0,'r*','markers',12);
    title([Name ':rRMSE:' num2str(RMSE) ':estVT:' num2str(estVT) ':actVT:' num2str(actVT) ':RR:' num2str(RR)], 'Interpreter', 'none'), 
    legend({'Predicted', 'Actual','X-Intercepts'}),
    xlabel('Samples'), 
    ylabel('Flow(L/s)');

end

function x0 = xintercepts(Pred,TimeF)
    x0 = [];
    for n=1:numel(Pred)-1
       
       if (Pred(n)>0 && Pred(n+1) < 0) || (Pred(n)<0 && Pred(n+1)>0)
           m = (Pred(n)-Pred(n+1))/(TimeF(n)-TimeF(n+1));
           b = Pred(n)-(TimeF(n)*m);
           x = -b/m;
           x0 = [x0;x];
       elseif Pred(n)==0
           x0 = [x0;TimeF(n)];    
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

% % returns the Estimated Tidal Volume of a trace based on inhales while 
% % evaluating inhales for badness
function a  = getVT2(FlowTS, Duration, x0)
   data = FlowTS.Data;
   assignin('base','data',data);
   FlowTS.Data = data.*Duration;
    Len = numel(x0);
    rem = mod(Len,2);
    if Len < 7
        a = 0;
        return;
    end
    if rem > 0
       Len = Len-1; 
    end
    VT = [];
    prev = 0;
    pdur = 0;
    for n=1:2:Len-2
        inhale = 0;
        duration=0;
% %     check if enough inhales have been found to calc VT
        if numel(VT)>3
            break;
        end
% %    Volume using Trapz
% %     find inhales
        temp = getsampleusingtime(FlowTS,x0(n),x0(n+1));
        duration = x0(n+1)-x0(n);
        if temp.Data(2)<0
            inhale = trapz(abs(temp.Data(:)));
        else
            temp = getsampleusingtime(FlowTS, x0(n+1),x0(n+2));
            inhale = trapz(abs(temp.Data(:)));
            duration = x0(n+2)-x0(n+1);
        end
% %     determine badness of inhale using 200mL as a threshold
        if inhale<(0.5*mean(VT))
            continue;
        elseif duration < (pdur*0.2)
            continue;
        elseif prev==0
            VT = [VT;inhale];
            prev = inhale;
            pdur = duration;
            continue;
        elseif abs(prev-inhale)<0.2
            VT = [VT;inhale];
            prev = inhale;
            pdur = duration;
            continue;
        elseif abs(mean(VT)-inhale)<0.2
            VT = [VT;inhale];
            prev = inhale;
            pdur = duration;
        end
            

        
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