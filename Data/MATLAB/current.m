

A=[[bike.rmsAx] [bike.rmsAy] [bike.rmsAz] [bike.rmsGx] [bike.rmsGy] [bike.rmsGz] [run.rmsAx] [run.rmsAy] [run.rmsAz] [run.rmsGx] [run.rmsGy] [run.rmsGz] [swim.rmsAx] [swim.rmsAy] [swim.rmsAz] [swim.rmsGx] [swim.rmsGy] [swim.rmsGz] [transition.rmsAx] [transition.rmsAy] [transition.rmsAz] [transition.rmsGx] [transition.rmsGy] [transition.rmsGz]]
labels=["bikeAx", "bikeAy", "bikeAz", "bikeGx", "bikeGy", "bikeGz", "runAx", "runAy", "runAz", "runGx", "runGy", "runGz", "swimAx","swimAy","swimAz", "swimGx" "swimGy", "swimGz", "transAx", "transAy", "transAz", "transGx", "transGy", "transGz"]
boxplot(A, labels)
set(gca,'FontSize',10,'XTickLabelRotation',90)
xlabel('Feature by axis and activity')
ylabel('Amplitude')
title('Root Mean Square')


%%
A=[[bike.meanAx] [bike.meanAy] [bike.meanAz] [bike.meanGx] [bike.meanGy] [bike.meanGz] [run.meanAx] [run.meanAy] [run.meanAz] [run.meanGx] [run.meanGy] [run.meanGz] [swim.meanAx] [swim.meanAy] [swim.meanAz] [swim.meanGx] [swim.meanGy] [swim.meanGz] [transition.meanAx] [transition.meanAy] [transition.meanAz] [transition.meanGx] [transition.meanGy] [transition.meanGz]]
labels=["bikeAx", "bikeAy", "bikeAz", "bikeGx", "bikeGy", "bikeGz", "runAx", "runAy", "runAz", "runGx", "runGy", "runGz", "swimAx","swimAy","swimAz", "swimGx" "swimGy", "swimGz", "transAx", "transAy", "transAz", "transGx", "transGy", "transGz"]
boxplot(A, labels)
set(gca,'FontSize',10,'XTickLabelRotation',90)
xlabel('Feature by axis and activity')
ylabel('Amplitude')
title('Mean')


%%
A=[[bike.stddevAx] [bike.stddevAy] [bike.stddevAz] [bike.stddevGx] [bike.stddevGy] [bike.stddevGz] [run.stddevAx] [run.stddevAy] [run.stddevAz] [run.stddevGx] [run.stddevGy] [run.stddevGz] [swim.stddevAx] [swim.stddevAy] [swim.stddevAz] [swim.stddevGx] [swim.stddevGy] [swim.stddevGz] [transition.stddevAx] [transition.stddevAy] [transition.stddevAz] [transition.stddevGx] [transition.stddevGy] [transition.stddevGz]]
labels=["bikeAx", "bikeAy", "bikeAz", "bikeGx", "bikeGy", "bikeGz", "runAx", "runAy", "runAz", "runGx", "runGy", "runGz", "swimAx","swimAy","swimAz", "swimGx" "swimGy", "swimGz", "transAx", "transAy", "transAz", "transGx", "transGy", "transGz"]
boxplot(A, labels)
set(gca,'FontSize',10,'XTickLabelRotation',90)
xlabel('Feature by axis and activity')
ylabel('Amplitude')
title('Standard Deviation')

%%
A=[[bike.mediAx] [bike.mediAy] [bike.mediAz] [bike.mediGx] [bike.mediGy] [bike.mediGz] [run.mediAx] [run.mediAy] [run.mediAz] [run.mediGx] [run.mediGy] [run.mediGz] [swim.mediAx] [swim.mediAy] [swim.mediAz] [swim.mediGx] [swim.mediGy] [swim.mediGz] [transition.mediAx] [transition.mediAy] [transition.mediAz] [transition.mediGx] [transition.mediGy] [transition.mediGz]]
labels=["bikeAx", "bikeAy", "bikeAz", "bikeGx", "bikeGy", "bikeGz", "runAx", "runAy", "runAz", "runGx", "runGy", "runGz", "swimAx","swimAy","swimAz", "swimGx" "swimGy", "swimGz", "transAx", "transAy", "transAz", "transGx", "transGy", "transGz"]
boxplot(A, labels)
set(gca,'FontSize',10,'XTickLabelRotation',90)
xlabel('Feature by axis and activity')
ylabel('Amplitude')
title('Median')


%%
A=[[bike.avgdisAx] [bike.avgdisAy] [bike.avgdisAz] [bike.avgdisGx] [bike.avgdisGy] [bike.avgdisGz] [run.avgdisAx] [run.avgdisAy] [run.avgdisAz] [run.avgdisGx] [run.avgdisGy] [run.avgdisGz] [swim.avgdisAx] [swim.avgdisAy] [swim.avgdisAz] [swim.avgdisGx] [swim.avgdisGy] [swim.avgdisGz] [transition.avgdisAx] [transition.avgdisAy] [transition.avgdisAz] [transition.avgdisGx] [transition.avgdisGy] [transition.avgdisGz]]
labels=["bikeAx", "bikeAy", "bikeAz", "bikeGx", "bikeGy", "bikeGz", "runAx", "runAy", "runAz", "runGx", "runGy", "runGz", "swimAx","swimAy","swimAz", "swimGx" "swimGy", "swimGz", "transAx", "transAy", "transAz", "transGx", "transGy", "transGz"]
boxplot(A, labels)
set(gca,'FontSize',10,'XTickLabelRotation',90)
xlabel('Feature by axis and activity')
ylabel('Amplitude')
title('Avg Dis b/ Peaks')


%%
A=[[bike.stddisAx] [bike.stddisAy] [bike.stddisAz] [bike.stddisGx] [bike.stddisGy] [bike.stddisGz] [run.stddisAx] [run.stddisAy] [run.stddisAz] [run.stddisGx] [run.stddisGy] [run.stddisGz] [swim.stddisAx] [swim.stddisAy] [swim.stddisAz] [swim.stddisGx] [swim.stddisGy] [swim.stddisGz] [transition.stddisAx] [transition.stddisAy] [transition.stddisAz] [transition.stddisGx] [transition.stddisGy] [transition.stddisGz]]
labels=["bikeAx", "bikeAy", "bikeAz", "bikeGx", "bikeGy", "bikeGz", "runAx", "runAy", "runAz", "runGx", "runGy", "runGz", "swimAx","swimAy","swimAz", "swimGx" "swimGy", "swimGz", "transAx", "transAy", "transAz", "transGx", "transGy", "transGz"]
boxplot(A, labels)
set(gca,'FontSize',10,'XTickLabelRotation',90)
xlabel('Feature by axis and activity')
ylabel('Amplitude')
title('Std Dev of Dis b/ Peaks')


%%
A=[[bike.avgampAx] [bike.avgampAy] [bike.avgampAz] [bike.avgampGx] [bike.avgampGy] [bike.avgampGz] [run.avgampAx] [run.avgampAy] [run.avgampAz] [run.avgampGx] [run.avgampGy] [run.avgampGz] [swim.avgampAx] [swim.avgampAy] [swim.avgampAz] [swim.avgampGx] [swim.avgampGy] [swim.avgampGz] [transition.avgampAx] [transition.avgampAy] [transition.avgampAz] [transition.avgampGx] [transition.avgampGy] [transition.avgampGz]]
labels=["bikeAx", "bikeAy", "bikeAz", "bikeGx", "bikeGy", "bikeGz", "runAx", "runAy", "runAz", "runGx", "runGy", "runGz", "swimAx","swimAy","swimAz", "swimGx" "swimGy", "swimGz", "transAx", "transAy", "transAz", "transGx", "transGy", "transGz"]
boxplot(A, labels)
set(gca,'FontSize',10,'XTickLabelRotation',90)
xlabel('Feature by axis and activity')
ylabel('Amplitude')
title('Avg Amp of Peaks')


%%
A=[[bike.stdampAx] [bike.stdampAy] [bike.stdampAz] [bike.stdampGx] [bike.stdampGy] [bike.stdampGz] [run.stdampAx] [run.stdampAy] [run.stdampAz] [run.stdampGx] [run.stdampGy] [run.stdampGz] [swim.stdampAx] [swim.stdampAy] [swim.stdampAz] [swim.stdampGx] [swim.stdampGy] [swim.stdampGz] [transition.stdampAx] [transition.stdampAy] [transition.stdampAz] [transition.stdampGx] [transition.stdampGy] [transition.stdampGz]]
labels=["bikeAx", "bikeAy", "bikeAz", "bikeGx", "bikeGy", "bikeGz", "runAx", "runAy", "runAz", "runGx", "runGy", "runGz", "swimAx","swimAy","swimAz", "swimGx" "swimGy", "swimGz", "transAx", "transAy", "transAz", "transGx", "transGy", "transGz"]
boxplot(A, labels)
set(gca,'FontSize',10,'XTickLabelRotation',90)
xlabel('Feature by axis and activity')
ylabel('Amplitude')
title('Std Dev of Amp of Peaks')


%%
A=[[bike.skewAx] [bike.skewAy] [bike.skewAz] [bike.skewGx] [bike.skewGy] [bike.skewGz] [run.skewAx] [run.skewAy] [run.skewAz] [run.skewGx] [run.skewGy] [run.skewGz] [swim.skewAx] [swim.skewAy] [swim.skewAz] [swim.skewGx] [swim.skewGy] [swim.skewGz] [transition.skewAx] [transition.skewAy] [transition.skewAz] [transition.skewGx] [transition.skewGy] [transition.skewGz]]
labels=["bikeAx", "bikeAy", "bikeAz", "bikeGx", "bikeGy", "bikeGz", "runAx", "runAy", "runAz", "runGx", "runGy", "runGz", "swimAx","swimAy","swimAz", "swimGx" "swimGy", "swimGz", "transAx", "transAy", "transAz", "transGx", "transGy", "transGz"]
boxplot(A, labels)
set(gca,'FontSize',10,'XTickLabelRotation',90)
xlabel('Feature by axis and activity')
ylabel('Amplitude')
title('Skewness')

%%
A=[[bike.kurtAx] [bike.kurtAy] [bike.kurtAz] [bike.kurtGx] [bike.kurtGy] [bike.kurtGz] [run.kurtAx] [run.kurtAy] [run.kurtAz] [run.kurtGx] [run.kurtGy] [run.kurtGz] [swim.kurtAx] [swim.kurtAy] [swim.kurtAz] [swim.kurtGx] [swim.kurtGy] [swim.kurtGz] [transition.kurtAx] [transition.kurtAy] [transition.kurtAz] [transition.kurtGx] [transition.kurtGy] [transition.kurtGz]]
labels=["bikeAx", "bikeAy", "bikeAz", "bikeGx", "bikeGy", "bikeGz", "runAx", "runAy", "runAz", "runGx", "runGy", "runGz", "swimAx","swimAy","swimAz", "swimGx" "swimGy", "swimGz", "transAx", "transAy", "transAz", "transGx", "transGy", "transGz"]
boxplot(A, labels)
set(gca,'FontSize',10,'XTickLabelRotation',90)
xlabel('Feature by axis and activity')
ylabel('Amplitude')
title('Kurtosis')






%% Double Plot for Comparison
Folder = 'Z:\GitRepositories\stretch-sense\Data\SenseAppData\';
%filename = '/Users/justinschaffner/Desktop/GitRepositories/stretch-sense/Data/SenseAppData/CAP_2018-03-07542162368_U_R_SIDE.csv';

%GTfile='/Users/justinschaffner/Desktop/GitRepositories/stretch-sense/Data/SenseAppData/GT_2018-03-07542161920_U_F_CENTER.csv';
%SGTfile='/Users/justinschaffner/Desktop/GitRepositories/stretch-sense/Data/SenseAppData/GT_2018-03-07542162176_U_R_CORNER.csv';
capfiles = ["Xiphoid\CAP_2018-03-17_X_B_CENTER_EILEEN_31_5inch.csv","4thIntercostal\CAP_2018-03-17_4th_B_Center_EILEEN_34_5inch.csv"];
GTfiles = ["Xiphoid\GT_2018-03-17_X_B_CENTER_EILEEN_31_5inch.csv","4thIntercostal\GT_2018-03-17_4th_B_Center_EILEEN_34_5inch.csv"];
figure; hold on;

for i=1:length(capfiles)
    capfile=strcat(Folder,capfiles(i));
    GTfile = strcat(Folder,GTfiles(i));
    GT=readtable(GTfile);
    %Read data
    T=readtable(capfile);
    
    Ttime = GT{:,2};
    Tlabel = GT{:,1};
    time = T{:,2};
    cap = T{:,1};
    start=time(1);
    for n = 1:length(time)
       time(n) = time(n)-start;
    end
    %for n = 1:length(Ttime)
    %    Ttime(n) = (100 * (Ttime(n) - floor(Ttime(n)))) + (floor(Ttime(n)) * 60);
    % end

    start=Ttime(1);
    for n = 1:length(Ttime)
        Ttime(n)=Ttime(n)-start;
    end
    ts=timeseries(cap,time);

    %  Plot data
    plot(ts)

    % Find peaks in data
%     [pks,locs,widths,proms] = findpeaks(cap, time,'MinPeakDistance',2, 'SortStr','descend');
%     plot(time,cap,locs,pks,'o')

    %z=max(cap{:});
    %invertedcap= cellfun(@(x) z - x, cap, 'Un', false);
    %[val,vloc,vwidths,vproms] = findpeaks(invertedcap, time, 'MinPeakDistance',2,'SortStr','descend');
    % plot(time,cap,vlocs,val,'o')
    if i==1
        linetype = 'r--';
    else
        linetype = 'b:';
    end
    for n=1:length(Ttime)
%         line([Ttime(n) Ttime(n)],ylim,linetype);
%         text(Ttime(n),min(ylim),Tlabel(n),'Rotation',90);
         vline(Ttime(n),linetype,Tlabel(n));
    end
end
%     text(locs+10,pks,num2str((1:numel(pks))'))
title('Xiphoid vs 4th Intercostal::Back Center')
xlabel('timestamp');
ylabel('Capacitance (pF)');



legend('Xiphoid', '4th Intercostal');




%% Dynamic Time Warping 
Folder = 'Z:\GitRepositories\stretch-sense\Data\SenseAppData\';

capfiles = ["Xiphoid\CAP_2018-03-17_X_B_CENTER_EILEEN_31_5inch.csv","4thIntercostal\CAP_2018-03-17_4th_B_Center_EILEEN_34_5inch.csv"];
%GTfiles = ["Xiphoid\GT_2018-03-17_X_B_CENTER_EILEEN_31_5inch.csv","4thIntercostal\GT_2018-03-17_4th_B_Center_EILEEN_34_5inch.csv"];



fileOne=strcat(Folder,capfiles(1));
fileTwo=strcat(Folder,capfiles(2));
%Read data
One=readtable(fileOne);
Two=readtable(fileTwo);

tOne = One{:,2};
cOne = One{:,1};
start=tOne(1);
for n = 1:length(tOne)
   tOne(n) = tOne(n)-start;
end

tTwo = Two{:,2};
cTwo = Two{:,1};
start=tTwo(1);
for n = 1:length(tTwo)
   tTwo(n) = tTwo(n)-start;
end


tsOne = timeseries(cOne,tOne);
tsTwo = timeseries(cTwo,tTwo);

[d,i1,i2] = dtw(tsOne,tsTwo);

tsOneW = tsOne(i1);
tsTwoW = tsTwo(i2);

figure; hold on;
plot(tsOneW);
plot(tsTwoW);


title('Xiphoid vs 4th Intercostal::Back Center')
xlabel('timestamp');
ylabel('Capacitance (pF)');



legend('Xiphoid', '4th Intercostal');




%% Single Plot

Folder = 'Z:\GitRepositories\stretch-sense\Data\';
%MACfilename = '/Users/justinschaffner/Desktop/GitRepositories/stretch-sense/Data/SenseAppData/CAP_2018-03-07542162368_U_R_SIDE.csv';
Filename = 'RawSensorData\Xiphoid\CAP\CAP\CAP_2018-03-25_Britt_FC_Xiphoid_32_5inch.csv';
Gfilename = 'RawSensorData\Xiphoid\CAP\GT\GT_2018-03-25_Britt_FC_Xiphoid_32_5inch.csv';

%GTfile='/Users/justinschaffner/Desktop/GitRepositories/stretch-sense/Data/SenseAppData/GT_2018-03-07542161920_U_F_CENTER.csv';
%SGTfile='/Users/justinschaffner/Desktop/GitRepositories/stretch-sense/Data/SenseAppData/GT_2018-03-07542162176_U_R_CORNER.csv';
%capfiles = ["Xiphoid\CAP_2018-03-17_X_B_CENTER_EILEEN_31_5inch.csv","4thIntercostal\CAP_2018-03-17_4th_B_Center_EILEEN_34_5inch.csv"];
%GTfiles = ["Xiphoid\GT_2018-03-17_X_B_CENTER_EILEEN_31_5inch.csv","4thIntercostal\GT_2018-03-17_4th_B_Center_EILEEN_34_5inch.csv"];



    capfile=strcat(Folder,Filename);
    GTfile = strcat(Folder,Gfilename);
    
    %Read data
    T=readtable(capfile);
    GT=readtable(GTfile);
    
% %     % with heart rate
%     time = T{:,3};
%     cap = T{:,2};
%     heart = T{:,1};
    
%     % without heart rate
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

    start=Ttime(1);
    for n = 1:length(Ttime)
        Ttime(n)=Ttime(n)-start;
    end
    ts=timeseries(cap,time);
%     Hts = timeseries(heart,time);

    %  Plot data
%     subplot(2,1,1);
    plot(ts);
   
    
    
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
    linetype = {'g--','r--'}; %'g--','r--','c--','m--'};
    for n=1:length(Ttime)
%         line([Ttime(n) Ttime(n)],ylim,linetype);
          text(Ttime(n),min(ylim),Tlabel(n),'Rotation',90,'Color',[1,0,0]);
          vline(Ttime(n),linetype{2});
           
    end

%     text(locs+10,pks,num2str((1:numel(pks))'))
title('Effect of Talking on Capacitance Readings')
xlabel('Time(s)');
ylabel('Capacitance (pF)');
legend('Xiphoid'); %'Start/End', 'Move', 'Deep Breath', 'Misc');


% for n = 1:4
%     plot(NaN,NaN,linetype{n});
% end


 
%   subplot(2,1,2);
%  plot(Hts);
%  title('Effect of Talking on Heart Rate')
% xlabel('Time(s)');
% ylabel('Heart Rate(beats/min)');
% legend('Heart Rate'); %'Start/End', 'Move', 'Deep Breath', 'Misc');
%  subplot(2,1,2);
%  
% linetype = {'g--','r--'}; %'g--','r--','c--','m--'};
%     for n=1:length(Ttime)
% %         line([Ttime(n) Ttime(n)],ylim,linetype);
%           text(Ttime(n),min(ylim),Tlabel(n),'Rotation',90);
%           vline(Ttime(n),linetype{2});
%            
%     end

%% Stretch Test
path = '/Users/justinschaffner/Desktop/GitRepositories/stretch-sense/Data/StretchTesting.csv';

test=readtable(path);
Idist = test{:,2};
Icap = test{:,3};
Ddist = test{:,4};
Dcap = test{:,5};

IncrD =[0.0,0.0,0.0,0.0];
IncrC = [0.0,0.0,0.0,0.0];
Increrr = [0.0,0.0,0.0,0.0];
DecrD = [0.0,0.0,0.0,0.0];
DecrC = [0.0,0.0,0.0,0.0];
Decrerr = [0.0,0.0,0.0,0.0];
len = length(Idist);
errlen= len/4;
for j = 1:4
    idist = 0.0;
    icap = 0.0;
    ddist = 0.0;
    dcap = 0.0;
    ivar = zeros(errlen,1);
    dvar = zeros(errlen,1);
    count = 1;
    for i = j:4:len
        idist = idist + Idist(i);
        icap = icap + Icap(i);
        ivar(count) = Icap(i);
        ddist = ddist + Ddist(i);
        dcap = dcap + Dcap(i);
        dvar(count) = Dcap(i);
        count = count + 1;
    end    
    idist = idist / (length(Idist) / 4);
    IncrD(j)=idist;
    icap = icap / (length(Icap) / 4);
    IncrC(j) = icap;
    Increrr(j)= var(ivar);
    ddist = ddist / (length(Ddist) / 4);
    DecrD(j)=ddist;
    dcap = dcap / (length(Dcap) / 4);
    DecrC(j)=dcap;
    Decrerr(j)= var(dvar);
    
end

figure; hold on;
errorbar(IncrD,IncrC, Increrr);
errorbar(DecrD,DecrC, Decrerr);
title('Capacitance vs Distance')
xlabel('Distance(cm)');
ylabel('Capacitance (pF)');
legend('Increasing', 'Decreasing');







%% Convert Cap to Distance
Folder = 'Z:\GitRepositories\stretch-sense\Data\';
%MACfilename = '/Users/justinschaffner/Desktop/GitRepositories/stretch-sense/Data/SenseAppData/CAP_2018-03-07542162368_U_R_SIDE.csv';
Filename = "SenseAppData\Xiphoid\CAP_2018-03-25_Jan_Breathing_FC_Xiphoid_30inch.csv";
Gfilename = "SenseAppData\Xiphoid\GT_2018-03-25_Jan_Breathing_FC_Xiphoid_30inch.csv";
Sfilename = 'StretchTesting.csv';
%GTfile='/Users/justinschaffner/Desktop/GitRepositories/stretch-sense/Data/SenseAppData/GT_2018-03-07542161920_U_F_CENTER.csv';
%SGTfile='/Users/justinschaffner/Desktop/GitRepositories/stretch-sense/Data/SenseAppData/GT_2018-03-07542162176_U_R_CORNER.csv';
%capfiles = ["Xiphoid\CAP_2018-03-17_X_B_CENTER_EILEEN_31_5inch.csv","4thIntercostal\CAP_2018-03-17_4th_B_Center_EILEEN_34_5inch.csv"];
%GTfiles = ["Xiphoid\GT_2018-03-17_X_B_CENTER_EILEEN_31_5inch.csv","4thIntercostal\GT_2018-03-17_4th_B_Center_EILEEN_34_5inch.csv"];
figure; hold on;


    capfile=strcat(Folder,Filename);
    GTfile = strcat(Folder,Gfilename);
    Sfile = strcat(Folder,Sfilename);
    
    %Read data
    Tstretch=readtable(Sfile);
    Tcap=readtable(capfile);
    Tgt=readtable(GTfile);
    
    time = Tcap{:,2};
    cap = Tcap{:,1};
    start=time(1);
    for n = 1:length(time)
       time(n) = time(n)-start;
    end
%     for n = 1:length(Ttime)
%        Ttime(n) = (100 * (Ttime(n) - floor(Ttime(n)))) + (floor(Ttime(n)) * 60);
%     end

%     start=Ttime(1);
%     for n = 1:length(Ttime)
%         Ttime(n)=Ttime(n)-start;
%     end
    
    Ttime = Tgt{:,2};
    Tlabel = Tgt{:,1};
    
    Idist = Tstretch{:,2};
    Icap = Tstretch{:,3};
    Ddist = Tstretch{:,4};
    Dcap = Tstretch{:,5};
    ydist = [Idist; Ddist];
    xcap = [Icap; Dcap];
    Tmdl = table(xcap,ydist);
    %linear regression model
    mdl = fitlm(Tmdl);
    ypred = predict(mdl, cap);
    
    for n = 1:length(ypred)
        ypred(n)=ypred(n)-10;
    end
 
    
ts=timeseries(ypred,time);
    %  Plot data
    plot(ts)

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
    linetype = {'g--','r--','c--','m--'};
    for n=1:length(Ttime)
%         line([Ttime(n) Ttime(n)],ylim,linetype);
%         text(Ttime(n),min(ylim),Tlabel(n),'Rotation',90);
         vline(Ttime(n),linetype{2},Tlabel(n));
    end

%     text(locs+10,pks,num2str((1:numel(pks))'))
title('Female, 40yrs')
xlabel('Timestamp(s)');
ylabel('Distance(cm)');



for n = 1:4
    plot(NaN,NaN,linetype{n});
end

legend('Xiphoid', 'Start/End', 'Move', 'Deep Breath', 'Misc');




