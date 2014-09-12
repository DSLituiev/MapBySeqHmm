%==== A simulator script for the mapping process

close all; clear all; clc; %<= magic spells

%% add useful functions
addpath('D:\MATLABuserfunctions\mtimesx');
addpath('D:\MATLABuserfunctions\MinMaxSelection');
addpath('D:\MATLABuserfunctions');
%= add the parent path
addpath(cd(cd('..')))

addpath('../emission')

%% == Initialize parameters:
par.Nplants = 50;
% par.SNPperPlant = 10;
par.r_expected = 12; 
par.FalsePositives = 0;
par.SamplingPoints = 20;
par.cSNPpos = .5;
par.Selection = true;
par.T_max = 1; %<= maximal length of the path
par.Rec_max = par.T_max * 10; %<= maximal number of recombinations per one path

%% == generate a path 
%= tv : vector of mutation locations
%= xv : vector of the number of mutant plants
%= t0 : location of the mutation under selection
[tChr, fChr, t0, ~, T_max] = generatePath(par);


%% == plot the distribution of the waiting times
tau = plotWaitingTimes(tChr);
estWT = tau*par.Nplants;

fprintf('number of recombinant plants:\t%u\n', par.Nplants )
fprintf('number of recombination points:\t%u\n', numel(tChr))
fprintf('chromosome length:\t%g\n', T_max)
%% simulate sequencing
[tSampling, q, r, fAtSamplingPoints] = sequencePath(fChr, tChr, par);

%% == plot
if par.Selection
    fname = 'chromosome under selection';
else
    fname = 'stationary chromosome';
end

figure('name', fname) 
stairs(tChr, fChr, 'b', 'linewidth', 2.5)
hold on
plot(tSampling, fAtSamplingPoints, 'rs-','markersize',8,'MarkerFaceColor','r')
plot(tSampling, 2*par.Nplants*q./r, 'go','markersize',8,'MarkerFaceColor','g')
xlim([0, T_max])
ylim([0, 10*round(.15 * par.Nplants)])
xlabel('linkage, morgans')
legend({'k: chromosome chunks'; 'k: SNP genotyping'; 'f: mutant read frequency x 2\itN'})


%% Forward-backward
tic


AR = readDataVect(ones(1, numel(tSampling)), tSampling ,q,r);
ChrMap(1).nt = [0 1];
ChrMap(1).cM = [0 100];

AR.chrMap = ChrMap;

AR.pop = par.Nplants;
AR.Alpha = 17;
AR.emissionHandle = @(qq, rr, ff)emissionk0(qq, rr, AR.pop);
AR.calcEmission;
AR.run();

logPobsDy = AR.xPsel;
fprintf('F-B algorithm took\t%2.2g\ts\n', toc)
%% == plot
[~, maxind] = max(logPobsDy);

figure
stem(tSampling, logPobsDy,  'go', 'BaseValue', floor(min(logPobsDy))-2,...
    'markersize',3,'MarkerFaceColor','g')
hold on
plot(tSampling, logPobsDy, 'y-')
plot(tSampling(maxind), max(logPobsDy), 'bv')
plot(t0, max(logPobsDy), 'r^')

xlim([0, par.T_max])
ylim([ floor(min(logPobsDy(~isinf(logPobsDy)))), ceil(max(logPobsDy))] )


fprintf('t0 predicted:\t%g\n', tSampling(maxind))
fprintf('t0 set      :\t%g\n', t0)
