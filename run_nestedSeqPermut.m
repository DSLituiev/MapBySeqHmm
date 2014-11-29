close all; clear all; clc;
dbclear if warning
tic

GRAPH_EXPORT_TYPE = 'pdf'; % 'eps'
% DATA_PATH = '/media/Processing/seq/data';
DATA_PATH = './raw_data';
USERFNCT_PATH = './dependencies';
FIGURES_PATH = './figures';
REF_PATH = './reference';
addpath(genpath(USERFNCT_PATH));
savepath;
addpath('./utilites');
addpath('./plotting');
addpath('./emission');

%=  provide the reference ID, which is the name of the csv file without '.csv'
%= together with the backround genotype ID
%= known positions of the causative SNP can be also provided here for
%= further visualization
%= + number of individuals in the mapping population (N)
dataID = 'HL7/p889_20110125_HL7_Paired-rmdup-clipOverlap-q20-ems-annotation'; x0 = 5672441; chr0 = 1; N = 50;

linkageLoosening = 1; % 1./(0:0.01:1);

%%
disp(['=======  Processing data from the run ''', dataID, ''' ======='])


%= construct the path to the (primary experimental) data file
dataPath = fullfile(DATA_PATH, [dataID, '.csv'] );
%= extract the refenece reads if the reference ID is given:
dataID = constructOutName( dataID, linkageLoosening );
mkdir(fullfile(FIGURES_PATH, dataID))
copyfile('./ResultsInfo.txt', fullfile(FIGURES_PATH, dataID, 'README.txt'))

% [AR, ~] = readSequencingDataCsv(dataPath, 'noannotation');

[AR, annotation] = subtractBackGroundGenotype(dataPath);
AR.pop = N;
AR.filterFields('r', @(x)(x<quantile(AR.r, 0.98))); % mutant reads
AR.filterFields('q', @(x)(x>7)); % mutant reads
AR.filterFields('f', @(x)(x<.8)); % SNP ratio
AR.calcDxMin;
AR.unmix();

AR.set_cMaxX( fullfile(REF_PATH ,'TAIR10-chr-counts.dat') );
mapPath = fullfile(REF_PATH ,'ChrMap.mat');
AR.setLinkageMap(mapPath)
theta = 1e-2;
lambda1 = 1;
AR.setBetaBinomialEmission( theta, lambda1)

%%
nsp = nestedSeqPermutations(AR);

n_iter = 1000;
nests = 1;
[xnPsel, xnLogOdds] = nsp.iterate( n_iter, nests);
%%
AR.runHMM('silent', true);
%%
FIGURES_PATH = './figures/HL7_permutations';
% mkdir(FIGURES_PATH )
for cc = 1:5
    chr_inds = (nsp.chromosome == cc) ;
    data_suffix = sprintf('-permutations_chr%u_iter%u_nests%u', cc, n_iter, nests);
    %%
    plotBatchProbT(nsp.x(chr_inds), xnPsel(chr_inds, 1:end-1), 'P[selection]', 'nolog', 'Mb');
    hold on
    plot(nsp.x(chr_inds)*1e-6, xnPsel(chr_inds, end), 'g-', 'linewidth', 2)
    if cc == chr0
        plot(x0*1e-6*[1,1], get(gca, 'ylim'), 'g--')
    end
    fig(gcf)
    exportfig(gcf, fullfile(FIGURES_PATH, sprintf('p_sel%s', data_suffix)), 'format','eps', 'color', 'rgb')
    %%
    plotBatchProbT(nsp.x(chr_inds), xnLogOdds(chr_inds, 1:end-1), 'log Odds', 'nolog', 'Mb');
    hold on
    plot(nsp.x(chr_inds)*1e-6, xnLogOdds(chr_inds, end), 'g-', 'linewidth', 2)
    if cc == chr0
        plot(x0*1e-6*[1,1], get(gca, 'ylim'), 'g--')
    end
    fig(gcf)
    exportfig(gcf, fullfile(FIGURES_PATH, sprintf('log_odds%s', data_suffix)), 'format','eps', 'color', 'rgb')
end
% % figure
% hold on
% plot(nsp.x(chr_inds), AR.xLogOdds(chr_inds), 'g-')
% set(gca, 'ylim', [-50, 0]);

%AR.plotChromosomes('xPsel', 'yscale', 'lin', 'norm', true, 'figure', 'new');
%AR.plotChromosomes('xLogOdds', 'yscale', 'lin', 'figure', 'new');


