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