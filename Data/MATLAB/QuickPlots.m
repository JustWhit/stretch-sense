current = CAPWalking;
GT = GTWalking;
GTLabel = GT.Time;
GTTime = GT.Label;

start = current.Time(1);
current.Time = current.Time - start;
GTTime = GTTime-start;
[pred,predCI] = predict(mdl, current.pF);
figure;
hold on;
plot(current.Time,current.pF);title('Capacitance');xlabel('Time(s)');ylabel('Capacitance(pF)');
for n=1:length(GTTime)
%           text(GTTime(n),min(ylim),GTLabel(n));
          vline(GTTime(n),'r--');
end
figure;
subplot(2,1,1);
hold on;
plot(current.Time,pred);title('Distance');xlabel('Time(s)');ylabel('Distance(cm)');
for n=1:length(GTTime)
          %text(GTTime(n),min(ylim),GTLabel(n),'Rotation',90,'Color',[1 0 0]);
          vline(GTTime(n),'r--');
end
subplot(2,1,2);
plot(current.Time,current.HR);title('Heart Rate');xlabel('Time(s)');ylabel('Heart Rate(bpm)');
for n=1:length(GTTime)
          %text(GTTime(n),min(ylim),GTLabel(n),'Rotation',90,'Color',[1 0 0]);
          vline(GTTime(n),'r--');
end

%% repeatability

path = 'Y:\GitRepositories\stretch-sense\Data\Archive\StretchTesting.csv';

test=readtable(path);

Dist = [test{:,2}; test{:,4}];
Cap = [test{:,3}; test{:,5}];

Set = unique(Dist);
CV = [];
for n=1 : length(Set)
    M = [];
    for i=1:length(Dist)
        if Dist(i) == Set(n)
           M = [M;Cap(i)]; 
        end
        
    end
    temp = 100*(std(M)/mean(M));
    CV = [CV;[Set(n) temp]];
end

disp(CV);
disp('MEAN:');
disp(mean(CV));


%% Test 8 subplots

start = TimeTest4T4(1);
TimeTest4T4 = TimeTest4T4 - start;
TEST4CAPTS = timeseries(CapTest4T4,TimeTest4T4);
TEST4VOLTS = timeseries(VolTest4T4, TimeTest4T4);
TEST4CAPTS = getsampleusingtime(TEST4CAPTS,0,29);
TEST4VOLTS = getsampleusingtime(TEST4VOLTS, 0,29);

figure;
subplot(2,1,1);
plot(TEST4CAPTS);xlabel('Time(s)');ylabel('Capacitance(pF)');title('');
subplot(2,1,2);
plot(TEST4VOLTS); xlabel('Time(s)');ylabel('Volume(L)');title('');
disp(corrcoef(TEST4CAPTS.Data,TEST4VOLTS.Data));







