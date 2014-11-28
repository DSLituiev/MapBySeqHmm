close all; clear all; clc;
dbclear if warning
tic

GRAPH_EXPORT_TYPE = 'pdf'; % 'eps'
DATA_PATH = '/media/Processing/seq/data';
% DATA_PATH = './raw_data';
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
AR.set_cMaxX( fullfile(REF_PATH ,'TAIR10-chr-counts.dat') );

% [AR.xPrior, AR.maxHitGene, AR.maxHitEffect,...
%     AR.positionCDS, AR.effectAA, AR.effectCodone] = constructPriorStr(annotation);

figure
myhist(50, AR.f, 'r')
hold all
myhist(50, AR.f(AR.f>0.05), 'g')


AR.filterFields('q', @(x)(x>7)); % mutant reads

AR.calcDxMin;
f = AR.visualizeStat;
fig(f, 'width', 24)
exportfig(gcf, fullfile(FIGURES_PATH, dataID,'qualityCtrl'), 'format', GRAPH_EXPORT_TYPE, 'color', 'rgb')

% visualizeAnnotationStat(annotation)

AR.filterFields('r', @(x)(x<quantile(AR.r, 0.98))); % mutant reads
AR.filterFields('q', @(x)(x>7)); % mutant reads
AR.filterFields('f', @(x)(x<.8)); % SNP ratio

% AR = AR.filterFields('dx', @(x)(x>120)); % SNP spacing
mixtObj = AR.unmix('plot');
f = plotReadSpacingPDF(AR, mixtObj, 'all');
clear mixtObj
AR.filterDx();

AR.plotChromosomes('dx', 'yscale', 'log', 'figure', 'new', 'plotfun', @(x,y)plot(x,y, 'rx'));
AR.plotChromosomes('dxFiltered', 'yscale', 'log', 'figure', 'old', 'plotfun', @(x,y)plot(x,y, 'm-'));

%% plot SNP ratio
AR.clearPlots;

KERNEL = 0; % 81;
AR.plotSnpRatio(KERNEL)
fig(gcf, 'width', 24)

exportfig(gcf, fullfile('figures', dataID, sprintf('SNP_Ratio_k%u', KERNEL)), 'format',GRAPH_EXPORT_TYPE, 'color', 'rgb')

%% load the linkage map
mapPath = fullfile(REF_PATH ,'ChrMap.mat');
AR.setLinkageMap(mapPath)

%% general experimental constants:
AR.pop = N;
% %%
% if ~AR.flagWT
%     AR.unmix();
%     [theta, lambda1] = runSimpleEM_BetaBinomAndUniform(AR.pop.Pstat, double(AR.q),...
%         double(AR.r), N, 'v',...
%         'contribution', AR.contrib ,'errTol', 1e-5, 'contribution', AR.contrib);
% else
%     [theta, lambda1] = runSimpleEM_BetaBinomAndUniform(AR.pop.Pflat, double([AR.q; AR.qw]),...
%         double([AR.r; AR.rw]), N, 'v',...
%         'contribution', AR.contrib ,'errTol', 1e-4);
% end
% %
theta = 1e-2;
lambda1 = 1;

AR.emissionHandle = @(q, r, study)emissionMixBetaBinomial(q, r, AR.pop, theta, lambda1);
% AR.emissionHandle = @(q, r, study)emissionMixBetaBinomial(q, r, AR.pop.N, theta, gi1(1:AR.Mtot));
% AR.emissionHandle = @(q, r, study)emissionBetaBinomial(q, r, study, theta);
% AR.emissionHandle = @(qq, rr, ff)emissionk0(qq, rr, AR.pop);

AR.Alpha =  linkageLoosening;

AR.plotChromosomes('contrib', 'yscale', 'lin', 'norm', true, 'figure', 'new', 'plotfun', @(x,y)plot(x,y, 'r^', 'markersize', 6, 'linewidth', 2));

AR.runHMM();

AR.normalizeChromosomes;

AR.plotChromosomes('xPsel', 'yscale', 'lin', 'norm', true, 'figure', 'new');

AR.plotChromosomes('xPstat', 'yscale', 'lin', 'norm', true);

z = AR.calcMleZ('plot');
fig(gcf, 'width', 24)
exportfig(gcf, fullfile(FIGURES_PATH, dataID, 'zMLE'),'format',GRAPH_EXPORT_TYPE, 'color', 'rgb')

if exist('chr0', 'var')
    AR.plotChromosome2D(chr0);
end

AR.plotChromosomes('xLogOdds', 'yscale', 'lin', 'figure', 'new', 'yThr', 0);

fig(gcf, 'width', 24)
exportfig(gcf, fullfile(FIGURES_PATH, dataID,'LogLiOdds'), 'format', GRAPH_EXPORT_TYPE, 'color', 'rgb')

%% plot likelihood and posterior
AR.plotStemsLP( 'ylim', [-20,0])
if exist('chr0', 'var')
    ind = (AR.x==x0 & AR.chromosome == chr0);    
    AR.plotChromosomes('xPselNorm', 'yscale', 'lin', 'figure', 'old', 'plotfun',...
        @(x,y)plot(x,y, 'r^', 'markersize', 6, 'linewidth', 2), 'select', ind);
end
fig(gcf, 'width', 24)
exportfig(gcf, fullfile(FIGURES_PATH, dataID, 'LH-Posterior'), 'format',GRAPH_EXPORT_TYPE, 'color', 'rgb')

%%
AR.plotChromosomes('xPselNorm', 'yscale', 'lin', 'figure', 'new', 'yThr', 0);

AR.plotChromosomes('xPosteriorNorm', 'yscale',  'lin','figure', 'old', 'yThr', 0); %, 'ylim', [-10,0]);

% AR.plotChromosomes('contrib', 'yscale', 'lin', 'figure', 'new', 'yThr', 0);
%% print
numberOfHits = 2e6;
[~,dataName,ext] = fileparts(dataID);
outListName = [dataName,ext, '-out.csv'];
outListPath = fullfile(FIGURES_PATH, dataID, outListName);
SO_DICT_PATH = fullfile(REF_PATH , 'SO_terms.csv');
AR.printTopHits(outListPath, numberOfHits, 'xLogOdds', 'cutoffValue', -Inf, 'so', SO_DICT_PATH)

%% plot misc stuff for the selected chromosome: crowdedness
if exist('chr0', 'var')
    figure('name', 'dx and membership')
    scatter(AR.x(AR.chromosome == chr0)*1e-6, AR.dx(AR.chromosome == chr0), 2.5+2.5*AR.contrib(AR.chromosome == chr0), AR.contrib(AR.chromosome == chr0), 'o')
    set(gca, 'yscale', 'log')
    
    figure('name', 'membership')
    plot(AR.x(AR.chromosome == chr0)*1e-6, AR.contrib(AR.chromosome == chr0), 'bx-')
end
