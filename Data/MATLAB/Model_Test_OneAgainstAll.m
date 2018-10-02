
Folder = 'Y:\Justin\GitRepositories\stretch-sense\Data\ModelData';
FileList = dir(char(fullfile(Folder,'*.csv')));

Traces = {};

for n=1:length(FileList)
    
    Table = readtable(char(fullfile(Folder,FileList(n).name))); 
    disp(FileList(n).name);
    Time = [Table{:,1}];
    Sensors = [Table{:,2:4}];
    Time = alterDuplicates(Time);

%     CapF = sgolayfilt(Sensors(:,1),3,25); 
%     CapF = denoise(Sensors(:,1));
    CapF = Sensors(:,1);

    FlowF = Sensors(:,3);    
    VolF = Sensors(:,2);
    
% figure;hold on;plot(Sensors(:,1));plot(y_dft);plot(CapF);
    CapF = truncateTrace(CapF);
    M = mean(CapF);
    CapF=CapF-M;
    CapF = fourierFilterThat(CapF,Time(1:numel(CapF)));
    CapF = transpose(CapF);
    [CapF, durations] = offsetData(CapF,Time(1:numel(CapF)));

% figure;plot(Time(1:length(CapF)),CapF); title(FileList(n).name, 'Interpreter','none');xlabel('Time(s)');ylabel('Capacitance(pF)');    
    Length = numel(CapF);
    
    FlowF = FlowF(2:Length+1);
% figure;hold on;plot(Sensors(:,1));plot(f);plot(CapF);
    VolF = VolF(2:Length+1);
    
    TimeF = Time(2:Length+1);


    
    Traces(n,:) = {CapF,FlowF,VolF,TimeF,durations};
end

[nr,~] = size(Traces);

for n=1:nr
    

    CapF = Traces{n,1};
    FlowTS = timeseries(Traces{n,2},Traces{n,4});
    VOLTS = timeseries(Traces{n,3},Traces{n,4});
    
% %     Build Training Set from Data not used in Test n
    Features = [];
    for i=1:nr
       if i~=n
           Features = [Features; [Traces{i,1} Traces{i,2}]];
       end
    end
    
% %     Remove pos pos and neg neg variables
    c1 = Features(:,1)< 0;
    c2 = Features(:,2)<0;
    Cneg = all(c1 & c2,2);
    Features(Cneg,:) = [];
    c3 = Features(:,1)>0;
    c4 = Features(:,2)>0;
    Cpos = all(c3 & c4,2);
    Features(Cpos,:)=[];
    c5 = Features(:,1)==0;
    c6 = Features(:,2)==0;
    Czero = all(c5 & c6);
    Features(Czero,:)=[];

% % \\\Fourier Models
% % Train model
inputs = Features(:,1);
targets = Features(:,2);
y_dft = fit(inputs,targets,'fourier8');

% % Get Predictions
pred = y_dft(Traces{n,1});

ci = predint(y_dft,Traces{n,1});

% % Smooth Predictions
% pred = sgolayfilt(pred,3,25);
PredTS = timeseries(pred,Traces{n,4});
% % get x intercepts
x0 = xintercepts(pred,Traces{n,4});
estVT = getVT(PredTS,Traces{n,5},x0);
actVT = getVT(FlowTS,Traces{n,5},x0);
%VOLVT = getVTfromVOL(VOLTS,x0);
RR = getRR(x0);

plotActualVSPred(FlowTS,PredTS, FileList(n).name,x0, estVT,actVT,RR);





end


% % % Fourier Transform and filter signal
function f=fourierFilterThat(y,time)
    Fs = numel(y)/(time(end)-time(1));
    % Compute DFT of x
    Y = fft(y); 
    L = length(y);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    
    f = Fs*(0:(L/2))/L;
    [pks,locs,width,~]=findpeaks(P1,f);
    
    i = pickAPeak(pks,locs);
    Freq = locs(i);
    bandW = 2*width(i);
    RR = Freq*60;
    
    f = FouFilter(y',time(end)-time(1),Freq,bandW,1,0);

end

% % % Select a suitable peak frequency
function i = pickAPeak(pks,locs)
     [~,i] = max(pks);
    if(locs(i)<0.2)
        [~,idx]=sort(pks,'descend');
        i=idx(2);
    end

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
    
% % % Diff with Wavelet Transform
%      CapF = diff(Data);
%     [c,l] = wavedec(CapF,3,'db2');
%     CapF = appcoef(c,l,'db2');
%     Durations = diff(Time(1:numel(Data),:));
    
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
    Rem = mod(Len,2);
    if Rem == 0
        Len = Len -1;
    end
    if Len < 7
        a = 0;
        return;
    end
    VT = [];
    for n=1:2:(Len-2)
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
        assignin('base','temp',temp);
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

function saveToCSV(PreFix,FileName,Trace)
    Folder = 'Y:\Justin\GitRepositories\stretch-sense\Data\ExSensorSpiroData\Truncated';
    csvwrite(char(fullfile(Folder,strcat(PreFix, FileName))),Trace);
end




function ry=FouFilter(y,samplingtime,centerfrequency,frequencywidth,shape,mode)
% Fourier filter function for time-series signal vector y; 'samplingtime'
% is the total duration of sampled signal in sec, millisec, or microsec;
% 'centerfrequency' and 'frequencywidth' are the center frequency and width 
% of the filter in Hz, KHz, or MHz, respectively; 'Shape' determines the 
% sharpness of the cut-off. If shape = 1, the filter is Gaussian; as  
% shape increases the filter shape becomes more and more rectangular. 
% Set mode = 0 for band-pass filter, mode = 1 for band-reject (notch) filter.  
% FouFilter returns the filtered signal.
%
% Example: Sine wave in noisy background.
% First half is just noise; sine wave starts halfway through.
% xx=[0:.001:2*pi]';
% signal=sin(20*xx);
% noise=randn(size(xx));
% x=1:2*length(xx)';
% y=[noise;signal+noise]; % sine wave is added halfway through.
% SignalToNoiseRatio=std(signal)/std(noise)
% FilteredSignal=foufilter(y',1,20,100,5,0);
% subplot(2,1,1)
% plot(x,y);
% title('First half is just noise; sine wave starts halfway through')
% subplot(2,1,2)
% plot(x,FilteredSignal);
% title('Signal filtered with FouFilter.m')
%
%  T. C. O'Haver (toh@umd.edu),  version 1.5, May, 2007

center=centerfrequency*samplingtime; %  center harmonic (fourier component)
width=frequencywidth*samplingtime; %  width of filter (in harmonics)

fy=fft(y); % Fourier transform of signal
lft1=[1:(length(fy)/2)];
lft2=[(length(fy)/2+1):length(fy)];
% Compute filter shape
ffilter1=ngaussian(lft1,center+1,width,shape);
ffilter2=ngaussian(lft2,length(fy)-center+1,width,shape);
ffilter=[ffilter1,ffilter2];
if mode==1, ffilter=1-ffilter; end
if length(fy)>length(ffilter), ffilter=[ffilter ffilter(1)];end
ffy=fy.*ffilter;  % Multiply filter by Fourier transform of signal
ry=real(ifft(ffy)); % Recover filter signal from Fourier transform

end
function g = ngaussian(x,pos,wid,n)
%  ngaussian(x,pos,wid) = peak centered on x=pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  Shape is Gaussian when n=1, becomes more rectangular as n increases.
% Example: ngaussian([1 2 3],1,2,1) gives result [1.0000    0.5000    0.0625]
g = exp(-((x-pos)./(0.6005615.*wid)) .^(2*round(n)));
end
