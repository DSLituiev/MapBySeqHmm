%==== A simulator script for the mapping process
close all; clear all; clc; %<= magic spells

%% add useful functions
addpath('D:\MATLABuserfunctions\mtimesx');
addpath('D:\MATLABuserfunctions\MinMaxSelection');
addpath('D:\MATLABuserfunctions');
%= add the parent path
addpath(cd(cd('..')))
addpath('../emission')
FIGURES_PATH = 'sim_figures';
%% == Initialize parameters:
LINKAGE_LOSENING = 1;
xCausativeSNP = 0.2;

n_iter = 500;
n_crossovers_pop = 50;
chr_length = .5;
n_plants = 25;
r_expected = 12;
batchsim = batchSimPoolSeq(xCausativeSNP, n_plants, r_expected, n_crossovers_pop, chr_length);
batchsim.Selection = false;
%% == generate a path 
batchsim.iterate(n_iter, 'grid');
%% == plot the distribution of the waiting times
data_suffix = sprintf('-iter_%u-crossovers_%u-cM_%u-N_%u-r%u', n_iter, n_crossovers_pop,...
    round(100*chr_length), n_plants, r_expected);



f = figure('name', 'logOdds' );
plot( batchsim.n_t, batchsim.xnLogOdds, 'color', [1,1,1]*0.8 )

%% 
figure('name', 'P selection')
surf( bsxfun(@minus, batchsim.xnPsel, calcMarginal(batchsim.xnPsel,1) ), 'linestyle', 'none')
view(0, 90)

figure('name', 'logOdds')
surf( batchsim.xnLogOdds , 'linestyle', 'none')
view(0, 90)

xnPselNormalized = bsxfun(@minus, batchsim.xnPsel, calcMarginal(batchsim.xnPsel,1) );
[f, medP] = plotBatchProbT(batchsim.t, xnPselNormalized , 'P selection' );
hold all
a = 9.5;
yL = exp( - a*batchsim.t);
yR = exp( - a*(batchsim.T_max-batchsim.t));
y = a * ( yL + yR - 1 );
plot(batchsim.t, y , 'kx')
fig(f, 'width', 16)
exportfig(gcf, fullfile(FIGURES_PATH, sprintf('p_select%s', data_suffix)), 'format','eps', 'color', 'rgb')

figure
plot(batchsim.t, medP - y )
hold all
plot([0, batchsim.T_max], [0,0], 'k-')
fig(f, 'width', 16)
exportfig(gcf, fullfile(FIGURES_PATH, sprintf('p_select_median_residuals%s', data_suffix)), 'format','eps', 'color', 'rgb')


f = plotBatchProbT(batchsim.t, batchsim.xnLogOdds, 'logOdds' );
fig(f, 'width', 16)
exportfig(gcf, fullfile(FIGURES_PATH, sprintf('log_odds%s', data_suffix)), 'format','eps', 'color', 'rgb')
