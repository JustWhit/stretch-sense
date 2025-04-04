%%
path = 'Y:\GitRepositories\stretch-sense\Data\Archive\StretchTesting.csv';

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

X = [Icap; Dcap];
Y = [Idist; Ddist];
mdl = fitlm(X,Y);
% start = test_CapvsDist.Time(1);
% test_CapvsDist.Time = test_CapvsDist.Time-start;
% Test = timeseries(test_CapvsDist.pF,test_CapvsDist.Time);
[pred, predCI] = predict(mdl,Test.Data);
Pred8 = timeseries(pred,Test.Time);
figure;hold on;plot(Test);title('Test Data');xlabel('Time(s)');ylabel('Capacitance(pF)');
figure;hold on;plot(Pred8);plot(Video_GT.Time,Video_GT.Dist,'--');plot(Video_GT.Time,Video_GT.Dist,'r*');title('Model Output');xlabel('Time(s)');ylabel('Distance(cm)');