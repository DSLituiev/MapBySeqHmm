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
LINKAGE_LOSENING = 1;
xCausativeSNP = 0.2;

par = simPooledSeq(xCausativeSNP, 50, 12, 100);
par.Selection = false;
%% == generate a path 
par.generatePath();
%% == plot the distribution of the waiting times
tau = plotWaitingTimes(par.tChr);

fprintf('number of recombinant plants:\t%u\n', par.pop.N)
fprintf('number of recombination points:\t%u\n', numel(par.tChr))
fprintf('chromosome length:\t%g\n', par.T_max)
%% simulate sequencing
[tSampling, q, r, fAtSamplingPoints] = sequencePath(par);
%% plot hidden data and sampling/sequencing
f = plotSeqSim( par );
%% Forward-backward
tic
[xLogOdds, xkPout] = par.runHMM();
fprintf('F-B algorithm took\t%2.2g\ts\n', toc)
%% == plot
f = plotHmmInferenceRes( par );