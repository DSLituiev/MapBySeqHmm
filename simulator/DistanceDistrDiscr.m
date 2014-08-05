clear all
close all

N = 4000;
M = 200;
p = M/N;
repeat = 100;

 SamplingPoints = sort(randi(N,[M,repeat]));

% SamplingPoints = (N*sort(rand(M,repeat)));

Dist = floor(abs(diff(SamplingPoints)));

maxhist = 70;
[pEmp, x] = hist(Dist, 0:(maxhist+1) );
pEmp = mean(pEmp,2);
pTheor(:,1) = ((M-1) / N) *(1- (M-2)/N ).^x; %<== precise
pTheor(:,2) = (p) *( 1- p ).^x;              %<== approximate geometric distribution

sum(pTheor)

figure
bar(x, pEmp )
hold all
plot(x, M*pTheor(:,1),'r' )
plot(x, M*pTheor(:,2),'g' )

xlim([1, maxhist])
ylim([0, max(pEmp(2:end-1))])
