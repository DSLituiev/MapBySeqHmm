%%% Gillespie's Algorithm, Main Script, 1D
close all; clear all; %clc
addpath('C:\Users\Dmitri\Documents\MATLAB\Information_and_Statistics\H_M')

%== Parameters of the simulation:
Rec_max = 600; % maximal number of reactions per one path (also limited by max time, see loop)
T_max = 2;
TimePoints = 700;
SamplingPoints = 300;
rTotal = 1e4;
t0 = 200; 
r_expected = rTotal./SamplingPoints;
NFalsePositives = 3;
Nplants = 50;

%== the microscopic reaction rate
c = 1; 
%== transition matrices
S = [ 1, 0; 0, 1 ];    % 1->0
P = [ 0, 1; 1, 0 ];
Nmatr = int32(P - S);
%== Q-matrix
hmatr = Qmatrix(Nplants);
%== run Gillespie:
[Zmatr tmatr] = GillespieLoop(Rec_max, Nplants-NFalsePositives, Nmatr, T_max);

%% sampling
%== generate sampling points
SamplingT = T_max*sort(rand(SamplingPoints,1));
% == study the distribution of the intervals between the sampling points
Dist = diff(SamplingT);
dt = T_max/SamplingPoints/4;
[pEmp, t] = hist(Dist, (.5*dt):dt:(T_max/100));

tau = T_max/SamplingPoints;
pTheor =  SamplingPoints*(exp( - (t-.5*dt)./tau) - exp( - (t + .5*dt)./tau) );
%== plot the distribution
figure
bar(t', pEmp')
hold all
plot(t', pTheor', 'r-')
xlim([0,T_max/100])
%== sample
Z = bsxfun(@le, tmatr, SamplingT');
Zs = zeros(size(SamplingT));
% surf(double(Z'), 'linestyle', 'none'); view(0,90)

for  ii = 1:length(SamplingT)
    ind = find(Z(:,ii), 1, 'last');
    if ~isempty(ind)
        Zs(ii) = Zmatr(ind , 1 );
    else
        Zs(ii) = 0;
    end
end

%% emission
r = random('Poisson', r_expected, [SamplingPoints,1] );
q = bsxfun(@(x,y)random('binom',x,y), r, Zs/Nplants/2);
%== plot
figure 
stairs(tmatr, Zmatr(:,1), 'b')
hold on
plot(SamplingT, Zs, 'r')
plot(SamplingT, 2*Nplants*q./r, 'go','markersize',3,'MarkerFaceColor','g')
xlim([0, max(tmatr)])
ylim([0, Nplants+10])
xlabel('linkage, morgans')
