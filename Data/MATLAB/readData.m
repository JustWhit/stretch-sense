filename = 'Deep_Breaths.csv';

% Read data
T=readtable(filename);

num = T{:,3};
cap = T{:,1};

%  Plot data
figure;
plot(num,cap)
title('Deep Breaths')
xlabel('Sample number');
ylabel('Capacitance (pF)');

% Find peaks in data
[pks,locs,widths,proms] = findpeaks(cap, num,'MinPeakDistance',50, 'SortStr','descend');
plot(num,cap,locs,pks,'o')
hold on
text(locs+10,pks,num2str((1:numel(pks))'))