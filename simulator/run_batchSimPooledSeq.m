%==== A simulator script for the mapping process
close all; clear all; clc; %<= magic spells

%% add useful functions
depend_path = 'D:\MATLABuserfunctions';
addpath( fullfile(depend_path, 'mtimesx') );
addpath( fullfile(depend_path, 'MinMaxSelection'));
addpath(depend_path);

addpath( fullfile(depend_path, 'cm_and_cb_utilities') );
%= add the parent path
addpath(cd(cd('..')))
addpath('../emission')
FIGURES_PATH = 'sim_figures';
%% == Initialize parameters:
LINKAGE_LOSENING = 1;
xCausativeSNP = 0.2;

n_iter = 1000;
n_crossovers_pop = 50;
chr_length = 1.5;
n_plants = 25;
r_expected = 12;
batchsim = batchSimPoolSeq(xCausativeSNP, n_plants, r_expected, n_crossovers_pop, chr_length);
batchsim.Selection = false;
%% == generate a path 
batchsim.iterate(n_iter, 'grid');
%% == plot the distribution of the waiting times
data_suffix = sprintf('-iter_%u-crossovers_%u-cM_%u-N_%u-r%u', n_iter, n_crossovers_pop,...
    round(100*chr_length), n_plants, r_expected);


f_sel = figure('name', 'logOdds' );
plot( batchsim.n_t, batchsim.xnLogOdds, 'color', [1,1,1]*0.8 )

%% 
xnPselNormalized = bsxfun(@minus, batchsim.xnPsel, calcMarginal(batchsim.xnPsel,1) );

% figure('name', 'P selection')
% surf( xnPselNormalized, 'linestyle', 'none')
% view(0, 90)
% 
% figure('name', 'logOdds')
% surf( batchsim.xnLogOdds , 'linestyle', 'none')
% view(0, 90)


[f_sel, medP, upper_quantile_P] = plotBatchProbT(batchsim.t*100, xnPselNormalized , 'P selection' );
hold all
% a = 9.5; % -log10(batchsim.pop.Pstat(end))
% yL = exp( - a*batchsim.t);
% yR = exp( - a*(batchsim.T_max-batchsim.t));
% y = a * ( yL + yR - 1 );
% plot(batchsim.t, y , 'kx')
fig(f_sel, 'width', 16)
% exportfig(gcf, fullfile(FIGURES_PATH, sprintf('p_select%s', data_suffix)), 'format','eps', 'color', 'rgb')


figure
plot( batchsim.t*100, upper_quantile_P - flipud(upper_quantile_P) )

% f_res = figure;
% plot(batchsim.t, medP - y )
% hold all
% plot([0, batchsim.T_max], [0,0], 'k-')
% fig(f_res, 'width', 16)
% exportfig(gcf, fullfile(FIGURES_PATH, sprintf('p_select_median_residuals%s', data_suffix)), 'format','eps', 'color', 'rgb')

[f_odd, ~, upper_quantile_L] = plotBatchProbT(batchsim.t, batchsim.xnLogOdds, 'logOdds' );

fig(f_odd, 'width', 16)
% exportfig(gcf, fullfile(FIGURES_PATH, sprintf('log_odds%s', data_suffix)), 'format','eps', 'color', 'rgb')

qu = 0.05;
for ii = 8:-1:2
    nn = fix(n_iter / ii);
    P = [batchsim.xnLogOdds(:, 1:nn), flipud(batchsim.xnLogOdds(:, nn+1:2*nn)) ];
    medianP_odds(:,ii) = median( P, 2 );
    upper_quantile_odds(:, ii) = quantile( P, 1-qu, 2 );
end

figure('name', 'upper quantile')
plot( batchsim.t*100, upper_quantile_odds)
cmlines()

figure('name', 'median')
plot( batchsim.t*100, medianP_odds)
cmlines()

figure('name', 'median - flipped median')
plot( batchsim.t*100, medianP_odds - flipud(medianP_odds))
cmlines()

figure
plot( batchsim.t*100, upper_quantile_L - flipud(upper_quantile_L) )

% obj.T(:,:,1).*diag(obj.pop.Pstat)

% batchsim.HMM.T(:,:,1)*diag(batchsim.pop.Pstat)