close all; clear all; clc;
dbclear if warning
tic


DATA_PATH = '/media/Processing/seq/data';
% DATA_PATH = '/media/Processing/seq/olddata';
% DATA_PATH = './data';
USERFNCT_PATH = './dependencies';
% addpath(fullfile(USERFNCT_PATH, '/binomial') );
addpath(USERFNCT_PATH);
addpath(fullfile(USERFNCT_PATH, 'mtimesx'));
addpath(fullfile(USERFNCT_PATH, 'MinMaxSelection'));
addpath(fullfile(USERFNCT_PATH, 'newtonraphson'));
 savepath;
addpath('./utilites');
addpath('./plotting');
addpath('./emission');

%=  provide the reference ID, which is the name of the csv file without '.csv'
%= together with the backround genotype ID
%= known positions of the causative SNP can be also provided here for
%= further visualization

%= number of plants:
% dataID = '/MO/20140514.A-MO7-bwa-nm2-ems-annotation'; N = 181; chr0 = 1;

% dataID = '/MO/20140514.A-MO7-bq20-ems-annotation-ecotypeInfo';  N = 181; chr0 = 1;
% dataID = '/MO/20140514.A-MO7-bq20-ems-annotation'; N = 181; chr0 = 1;

dataID = '/MO/20140514.A-MO8-bq20-ems-annotation-ecotypeInfo'; N = 278; chr0 = 5;

% dataID = '/MO/20140514.A-MO8-bq20-ems-annotation'; N = 278; chr0 = 5;
% dataID = '20140514.A-MO7-bwa-rmdup-clipOverlap-q20-freebayes-ems-annotation-rna'; N = 181; chr0 = 1;
% dataID = 'MO8-rmdup-clipOverlap-q20-freebayes-ems-annotation-array-rna'; N = 278;  chr0 = 5;
% dataID = 'MO7-rmdup-clipOverlap-q20-freebayes-ems-annotation-array-rna'; N = 181; chr0 = 1;

%%
%= linkage loosening factor:
ALPHA = 1;
% ALPHA = 1./(0:0.01:1);
%= construct the path to the (primary experimental) data file
dataPath = fullfile(DATA_PATH, [dataID, '.csv'] );
disp(['=======  Processing data from the run ''', dataID, ''' ======='])
%= and genetic distance in cM between them)

dataID = sprintf(strcat(dataID, '_%u'), ALPHA);
mkdir(fullfile('figures',dataID))

clear AR 

[AR] = readSequencingDataCsv2(dataPath);
%% apply filters and visualise statistics
AR = AR.filter('q', @(x)(x>7)); % mutant reads
AR = AR.filter('f', @(x)(x<.8)); % SNP ratio

AR.calcDxMin;
AR = AR.filter('dx', @(x)(x>10), 'chromosome', @(x)(x==chr0)); % 
AR = AR.filter('dx', @(x)(x<1e4)); % 


AR.calcDxMin;
AR.visStat;
fig(gcf, 'width', 24)
exportfig(gcf, fullfile('figures', dataID,'qualityCtrl'), 'format','eps', 'color', 'rgb')

% visualizeAnnotationStat(annotation)

% AR.unmix('plot');
%% load the recombination map 
%= (contains positions of the markers
%= and genetic distance in cM between them)
load('./reference/ChrMap.mat')
AR.chrMap = ChrMap;
clear ChrMap;
%% plot SNP ratio

addpath(fullfile(USERFNCT_PATH, 'fastmedfilt1d'));
KERNEL = 81;

AR.plotSnpRatio(KERNEL, chr0)
% AR.xTemp = AR.snpEcotypesInfo & AR.r > 40 & AR.r > 60;
% AR.plotSnpRatio(KERNEL, chr0, 'xTemp')

fig(gcf, 'width', 24)
exportfig(gcf, fullfile('figures', dataID, sprintf('SNP_Ratio_k%u', KERNEL)), 'format','eps', 'color', 'rgb')

% XRANGE = uint32([23,24]*1e6);
XRANGE = [7 , 8]*1e6;
set(gca, 'xlim', XRANGE)
set(gca,'XTickLabelMode','auto')
set(gca, 'ylim', [-0.05, .55])
set(gca,'XTick',XRANGE(1):1e5:XRANGE(2))
exportfig(gcf, fullfile('figures', dataID, sprintf('SNP_Ratio_k%u-range_%1.1e-%1.1e.eps', KERNEL,XRANGE(1), XRANGE(2))), 'format','eps', 'color', 'rgb')


AR.plotScatterMW()
fig(gcf, 'width', 24)
exportfig(gcf, fullfile('figures', dataID, 'SNP_Ratio_WT_vs_MT.eps'), 'format','eps', 'color', 'rgb')

return
%% set general experimental constants:
AR.pop = N;
%%

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

theta = 1e-2;
lambda1 = 1-1e-3;

% AR.includeRnaPresence()

AR.emissionHandle = @(q, r, study)emissionMixBetaBinomial(q, r, AR.pop, theta, lambda1);
% AR.emissionHandle = @(q, r, study)emissionMixBetaBinomial(q, r, AR.pop.N, theta, gi1(1:AR.Mtot));

%  emissionHandle = @(q, r, study)emissionBetaBinomial(q, r, study, theta);
% AR.emissionHandle = @(qq, rr, ff)emissionk0(qq, rr, AR.pop);
% AR.contrib = ones(size(AR.x));


AR.Alpha = ALPHA;
% 
% load( 'emission_fun.mat')
%  AR.E = E;
%  AR.contrib = contrib;
%  AR.contrib = 1;really
AR.calcEmission;
AR.run();

figure; plot(AR.xPsel)


AR.clearPlots;

AR.normalizeChromosomes;

AR.plotChromosomes('xPsel', 'yscale', 'lin', 'norm', true, 'figure', 'new');

AR.plotChromosomes('xPstat', 'yscale', 'lin', 'norm', true);

% AR.calcLogOdds;
AR.plotChromosomes('xLogOdds', 'yscale', 'lin', 'figure', 'new', 'yThr', 0);

fig(gcf, 'width', 24)
exportfig(gcf, fullfile('figures', dataID,'LogLiOdds'), 'format','eps', 'color', 'rgb')

%% LH and posterior
AR.plotStemsLP(chr0, 'ylim', [-20,0])
fig(gcf, 'width', 24)
exportfig(gcf, fullfile('figures', dataID, 'LH-Posterior'), 'format','eps', 'color', 'rgb')

% XRANGE = [23,24]*1e6;
set(gca, 'xlim', XRANGE)
set(gca,'XTickLabelMode','auto')
set(gca, 'ylim', [-4, 0])
set(gca,'XTick',XRANGE(1):1e5:XRANGE(2))
exportfig(gcf, fullfile('figures', dataID, sprintf('LH-Posterior-range_%1.1e-%1.1e.eps', XRANGE(1), XRANGE(2))), 'format','eps', 'color', 'rgb')
%%

AR.plotChromosomes('xPselNorm', 'yscale', 'lin', 'figure', 'new', 'yThr', 0);

AR.plotChromosomes('xPosteriorNorm', 'yscale',  'lin','figure', 'old', 'yThr', 0); %, 'ylim', [-10,0]);
%% dx
AR.plotOneChromosome(chr0, 'dx', 'yscale', 'log', 'figure', 'new', 'plotfun', @(x,y,z)plot(x,y,'r.'));

fig(gcf, 'width', 24)
exportfig(gcf, fullfile('figures', dataID, 'log_dx'), 'format','eps', 'color', 'rgb')

XRANGE = [23,24]*1e6;
set(gca, 'xlim', XRANGE)
set(gca,'XTickLabelMode','auto')
set(gca,'XTick',XRANGE(1):1e5:XRANGE(2))

set(gca, 'ylim', [0, 5e3])
exportfig(gcf, fullfile('figures', dataID, sprintf('log_dx-range_%1.1e-%1.1e.eps', XRANGE(1), XRANGE(2))), 'format','eps', 'color', 'rgb')
%%
%%

% AR.plotChromosomes('contrib', 'yscale', 'lin', 'figure', 'new', 'yThr', 0);


figure
plot(AR.x(AR.chromosome == chr0)*1e-6, AR.dx(AR.chromosome == chr0), 'rx')
set(gca, 'yscale', 'log')


% 
% x0 = AR.x(AR.chromosome == chr0)*1e-6 ;
% figure; 
% plot(x0,  [fm81, fw81])
% hold all
% [~, peakInd] = maxk(fm81 - fw81, 10 );
% clutteredSubs = find(abs(diff(peakInd))<KERNEL);
% plot(x0(peakInd)*[1, 1], [0, 0.5])
fig(gcf, 'width', 24)
exportfig(gcf, fullfile('figures', dataID, 'median filtered SNP Ratio'), 'format','eps', 'color', 'rgb')

%%
%%

% AR.plotStems('xPosteriorNorm', 'figure', 'new');
% AR.plotStems('xPselNorm', 'yThr', 0);
numberOfHits = 2e6;


[~,dataName,ext] = fileparts(dataID);
AR.printTopHits(fullfile('figures', dataID, [dataName,ext, '-out.csv']), numberOfHits, 'xLogOdds', 'cutoffValue', 0)

% 
% E = AR.E;
% contrib = AR.contrib;
% save('emission_obj', 'E', 'contrib')

finishing_beep(1, .2, 3)
