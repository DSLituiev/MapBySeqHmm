%%% Gillespie's Algorithm, Main Script, 1D

%=== Add FALSE POSITIVES option

close all; clear all; %clc

addpath('C:\Users\Dmitri\Documents\MATLAB\Information_and_Statistics\H_M')
SamplingPoints = 100; 
r_expected = 7;
% MaxT = 2000;
%== maximal number of reactions per one path (also limited by max time, see loop)
%== ! St_max must be even !!!
St_max = 16;

Nplants = 49; % if odd, generates stationary process

%== Generate the paths:
tmatr = cumsum(log(1./rand(St_max, Nplants,'single')));
statevector = false(St_max,1);
statevector(1:2:end) =  1;

[~, lastminind] = min(tmatr,[],2);
%== the cutoff time
T_max = tmatr( end , lastminind(end) );
[tvect, tind] = sort(tmatr(:));

%== plot the distribution of the waiting times
figure
hist(abs(diff(tvect)), (T_max/SamplingPoints/Nplants):(2*T_max/SamplingPoints/Nplants):(200*T_max/SamplingPoints/Nplants) )
xlim([0, (180*T_max/SamplingPoints/Nplants)])
xlabel('recombination distance')
ylabel('frequency')
%% assemble a stationary path
%== initialize the 'mutation matrix' 
%= indicating which chromosome stays after the meiosis
MutMatrCoeff = ones(size(tmatr));
%== assign the even columns (i.e. each second plant) '-1'
MutMatrCoeff(:, 2:2:end) = -1;
%== generate the state-derivtive matrix 
%= based on the order of the events:
%= events that appear on even 
DxMatr = zeros(size(tmatr));
DxMatr(tind) = (~mod(tind,2) + -mod(tind,2));
DxMatr = DxMatr.*MutMatrCoeff;
%== sum up the states changes in the order they occur:
x = sum(MutMatrCoeff(1,:)+1)/2 + cumsum( DxMatr(tind) );
%== plot the stationary distribution
figure
stairs(tvect, x)
hold on
plot([ 0, T_max ], Nplants./2*[1 1], 'r:' )
xlim([ 0, T_max ])
ylim([0, 10*round(0.1*Nplants)])

%% assemble a path under selection
%== select the points where at time t0 the state is 1
t0 = T_max/2;
t0ind = find(tvect >t0, 1,'first');
t0 = tvect(t0ind);
indt0 = zeros(Nplants,1);
for ii = 1:Nplants
    indt0(ii) = find( tmatr(:,ii) >t0, 1,'first');
   %  tmatr( indt0(ii),ii)
end


MutMatrCoeff = -bsxfun(@plus, ones(size(tmatr)), -2*mod(indt0,2)' );
DxMatr(tind) = (~mod(tind,2) + -mod(tind,2));
DxMatr = DxMatr.*MutMatrCoeff;

%== sum up the states changes in the order they occur:
cumdx = cumsum( DxMatr(tind) );
x = Nplants-cumdx(t0ind) + cumdx;
%== plot the distribution under selection
figure
stairs(tvect, x)
hold on
plot([ 0, T_max ], Nplants./2*[1 1], 'r:' )
xlim([ 0, T_max ])
ylim([0, 10*round(0.1*Nplants)])
%% sampling
%== generate sampling points
SamplingT = sort([T_max*rand(SamplingPoints-1,1); t0]);
% == study the distribution of the intervals between the sampling points
Dist = diff(SamplingT);
dt = T_max/SamplingPoints/10;
[pEmp, t] = hist(Dist, (.5*dt):dt:(T_max/10));

tau = T_max/SamplingPoints;
pTheor =  SamplingPoints*(exp( - (t-.5*dt)./tau) - exp( - (t + .5*dt)./tau) );
%== plot the distribution
figure
bar(t', pEmp')
hold all
plot(t', pTheor', 'r-')
 xlim([0,T_max/12])
%== sample
Z = bsxfun(@le, tvect, SamplingT');
Zs = zeros(size(SamplingT));
% surf(double(Z'), 'linestyle', 'none'); view(0,90)

for  ii = 1:length(SamplingT)
    ind = find(Z(:,ii), 1, 'last');
    if ~isempty(ind)
        Zs(ii) = x(ind  );
    else
        Zs(ii) = ceil(Nplants./2);
    end
end

%% emission
r = random('Poisson', r_expected, [SamplingPoints,1] );
q = bsxfun(@(x,y)random('binom',x,y), r, Zs/Nplants/2);
%== plot
figure 
stairs(tvect, x, 'b')
hold on
plot(SamplingT, Zs, 'rs-','markersize',3,'MarkerFaceColor','r')
plot(SamplingT, 2*Nplants*q./r, 'go','markersize',3,'MarkerFaceColor','g')
xlim([0, T_max])
ylim([0, 10*round(.15*Nplants)])
xlabel('linkage, morgans')


