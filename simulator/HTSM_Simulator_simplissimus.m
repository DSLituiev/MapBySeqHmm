%==== A simulator script for the mapping process

close all; clear all; clc; %<= magic spells

%% add useful functions
addpath('D:\MATLABuserfunctions\mtimesx');
addpath('D:\MATLABuserfunctions\MinMaxSelection');
addpath('D:\MATLABuserfunctions');
%= add the parent path
addpath(cd(cd('..')))

%% == Initialize parameters:
par.Nplants = 50;
% par.SNPperPlant = 10;
par.r_expected = 60; 
par.FalsePositives = 0;

par.T_max = 1; %<= maximal length of the path
par.SamplingPoints = 100 ;
par.M = par.SamplingPoints;
par.cSNPpos = 0.8;
par.Selection = true;
par.Rec_max = par.T_max * 10; %<= maximal number of recombinations per one path

study.N = par.Nplants;
study.kvect = (0:study.N)';

%% == generate a path 
%= tv : vector of mutation locations
%= xv : vector of the number of mutant plants
%= t0 : location of the mutation under selection

%% simulate sequencing
flag = 1;
while flag    
    [tChr, fChr, t0, ~, T_max] = generatePath(par);
    flag = any(size(fChr) == 0 );
    if ~flag
        [tSampling, q, r, fAtSamplingPoints] = sequencePath(fChr, tChr, par);
    end    
end
%% == plot the distribution of the waiting times
% tau = plotWaitingTimes(tChr);
% estWT = tau*par.Nplants;
% 
% fprintf('number of recombinant plants:\t%u\n', par.Nplants )
% fprintf('number of recombination points:\t%u\n', numel(tChr))
% fprintf('chromosome length:\t%g\n', T_max)


%% == plot
if par.Selection
    fname = 'chromosome under selection';
else
    fname = 'stationary chromosome';
end

figure('name', fname) 
stairs([0; tChr; par.T_max], [fChr(1); fChr; fChr(end)], 'b', 'linewidth', 2.5)
hold on
plot(tSampling, fAtSamplingPoints, 'rs-','markersize',3,'MarkerFaceColor','r')
plot(tSampling, 2*par.Nplants*q./r, 'go','markersize',3,'MarkerFaceColor','g')
xlim([0, T_max])
ylim([0, 10*round(.15 * par.Nplants)])
xlabel('linkage, morgans')
legend({'k: chromosome chunks'; 'k: SNP genotyping'; 'f: mutant read frequency x 2\itN'})

%% emission
fl.CutReadNumberBinom = 300;
fnames = {'q','r'};

for ii = 1:length(fnames)
    limReads.max.(fnames{ii}) = max( eval(fnames{ii}) ,[],1);
    limReads.min.(fnames{ii}) = min( eval(fnames{ii}) ,[],1);
end

if fl.CutReadNumberBinom
    limReads.max.q = min( limReads.max.q, uint32(fl.CutReadNumberBinom) );
    limReads.max.r = min( limReads.max.r, uint32(fl.CutReadNumberBinom) );
end

EmMatrix = initializeEmissionMatrix(limReads, study);

%% Forward-backward
stationaryFlag = 0;
tic
logPobsDy = propagateProbAlongChr(tSampling, q, r,  EmMatrix, limReads, study.N, stationaryFlag)

logPobsDy = logPobsDy - calcMarginal(logPobsDy);

% [scaledAcum, scaledBcum, scaleA, scaleB] = calcCumMatrices(par.M , par.Nplants,E,T);
% logPobs = calcLogPobs(M, N, CM, T, stationaryFlag, varargin);
diff(logPobsDy)


% runFBsmall_backup20131205(tSampling, q, r,  par.Nplants);
fprintf('F-B algorithm took\t%2.2g\ts\n', toc)
% return
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
ylim([ floor(min(logPobsDy)), ceil(max(logPobsDy))] )


fprintf('t0 predicted:\t%g\n', tSampling(maxind))
fprintf('t0 set      :\t%g\n', t0)
