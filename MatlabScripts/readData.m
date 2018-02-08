filename = 'Deep_Breaths.csv';

T=readtable(filename);

num = T{:,3};
cap = T{:,1};

figure;
plot(num,cap)
title('Deep Breaths')
xlabel('Sample number');
ylabel('Capacitance (pF)');