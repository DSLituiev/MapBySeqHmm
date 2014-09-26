close all; clear all; clc;
dbclear if warning
tic


DATA_PATH = '/media/Processing/seq/data';
% DATA_PATH = '/media/Processing/seq/olddata';
% DATA_PATH = './data';
USERFNCT_PATH = '/media/Processing/MATLABuserfunctions';
% addpath(fullfile(USERFNCT_PATH, 'MATLABuserfunctions/binomial') );
addpath(USERFNCT_PATH);
addpath(fullfile(USERFNCT_PATH, 'mtimesx'));
addpath(fullfile(USERFNCT_PATH, 'MinMaxSelection'));
addpath(fullfile(USERFNCT_PATH, 'newtonraphson'));
addpath(fullfile(USERFNCT_PATH, 'fastmedfilt1d'));
savepath;
addpath('./utilites');
addpath('./plotting');
addpath('./emission');

%=  provide the reference ID, which is the name of the csv file without '.csv'
%= together with the backround genotype ID
%= known positions of the causative SNP can be also provided here for
%= further visualization

%= number of plants:
dataID = 'CH/73-bwa-rmdup-clipOverlap-q20-nm2-ems-annotation'; N = 100;
% dataID = 'CH/73-bwa-rmdup-clipOverlap-q20-nm2-ems-annotation'; N = 100;
% dataID = 'CH/91-bwa-rmdup-clipOverlap-q20-nm2-ems-annotation'; N = 100;
% dataID = 'CH/91-ngm-rmdup-clipOverlap-q20-ems-annotation'; N = 100;
%  dataID = 'CH/39-ngm-rmdup-clipOverlap-q20-nm2-ems-annotation'; N = 100;
% dataID = 'CH/39-ngm-rmdup-clipOverlap-q20-ems-annotation'; N = 100;

% dataID = '91-ngm-rmdup-clipOverlap-q20-ems-annotation'; N = 100;
% dataID = 'ABD/20140605.X-ABD192-bwa-rmdup-clipOverlap-q20-nm2-ems-annotation';  N = 100;

% dataList = {'ABD/20120427-ABD159-ngm-rmdup-clipOverlap-q20-nm2-ems-annotation',   100,  2,  17521246 ;...
%             'ABD/20120427-ABD159-bwa-rmdup-clipOverlap-q20-nm2-ems-annotation',   100,  2,  17521246 ;...
%             'ABD/20130426.B-ABD173-ngm-rmdup-clipOverlap-q20-nm2-ems-annotation',  50,  3, 1619248; ...
%             'ABD/20130426.B-ABD173-bwa-rmdup-clipOverlap-q20-nm2-ems-annotation',  50 , 3, 1619248; ...
%             'HL10/HL10-rmdup-clipOverlap-nm2-q20-temp-ems-annotation', 50, 3 , 16473265; ...
%             'HL10/HL10-rmdup-clipOverlap-nm2-q20-ems-annotation',  50, 3 , 16473265; ...
%             'HL10/HL10-rmdup-clipOverlap-q20-freebayes-ems-annotation',  50, 3 , 16473265; ...
%             'HL7/p889_20110125_HL7_Paired-rmdup-clipOverlap-nm2-q20-temp-ems-annotation', 50, 1,  5672441; 
%             'HL7/p889_20110125_HL7_Paired-rmdup-clipOverlap-nm2-q20-ems-annotation', 50, 1,  5672441; 
%             'HL7/p889_20110125_HL7_Paired-rmdup-clipOverlap-q20-freebayes-ems-annotation', 50, 1, 5672441};
            
           
% dataID = 'ABD/20120427-ABD159-ngm-rmdup-clipOverlap-q20-nm2-ems-annotation';   N = 100; chr0 = 2; x0 = 17521246;
% dataID = 'ABD/20120427-ABD159-bwa-rmdup-clipOverlap-q20-nm2-ems-annotation'; chr0 = 2; x0 = 17521246;  N = 100;

%  dataID = 'ABD/20130426.B-ABD173-ngm-rmdup-clipOverlap-q20-nm2-ems-annotation'; chr0 = 3 ; x0 = 1619248;  N = 50;  % bkgrID = 'ABD241-rmdup-clipOverlap-q20-freebayes';
%  dataID = 'ABD/20130426.B-ABD173-bwa-rmdup-clipOverlap-q20-nm2-ems-annotation'; chr0 = 3 ; x0 = 1619248;  N = 50;  % bkgrID = 'ABD241-rmdup-clipOverlap-q20-freebayes';
 
% dataID = 'MO8_mutant_pool-rmdup-clipOverlap-q20-freebayes-ems-annotation-repfilt';

% dataID = 'HL4/p889_20110501_HL4_pairing-rmdup-clipOverlap-nm2-q20-ems-annotation';  chr0 =  3 ;  x0 =  16473265;N = 50;
% dataID = 'HL10/HL10-rmdup-clipOverlap-nm2-q20-temp-ems-annotation';  chr0 =  3 ;  x0 =  16473265;N = 50;
% dataID = 'HL10/HL10-rmdup-clipOverlap-nm2-q20-ems-annotation';  chr0 =  3 ;  x0 =  16473265;N = 50;
% dataID = 'HL10/HL10-rmdup-clipOverlap-q20-freebayes-ems-annotation';  chr0 =  3 ;  x0 =  16473265;N = 50;

% dataID = 'HL7/p889_20110125_HL7_Paired-rmdup-clipOverlap-nm2-q20-temp-ems-annotation'; x0 = 5672441; chr0 = 1; N = 50;
% dataID = 'HL7/p889_20110125_HL7_Paired-rmdup-clipOverlap-nm2-q20-ems-annotation'; x0 = 5672441; chr0 = 1; N = 50;
% dataID = 'HL7/p889_20110125_HL7_Paired-rmdup-clipOverlap-q20-freebayes-ems-annotation'; x0 = 5672441; chr0 = 1; N = 50;

    
ALPHA = 1;
% AR.Alpha = 1./(0:0.01:1);

%%
disp(['=======  Processing data from the run ''', dataID, ''' ======='])

%= and genetic distance in cM between them)

%=       load the recombination map (contains positions of the markers
%= and genetic distance in cM between them)

%= construct the path to the (primary experimental) data file
dataPath = fullfile(DATA_PATH, [dataID, '.csv'] );
%= extract the refenece reads if the reference ID is given:

dataID = sprintf(strcat(dataID, '_%u'), ALPHA);
mkdir(fullfile('figures',dataID))

% [AR, ~] = readSequencingDataCsv(dataPath, 'noannotation');

[AR, annotation] = subtractBackGroundGenotype(dataPath);
%
% [AR.xPrior, AR.maxHitGene, AR.maxHitEffect,...
%     AR.positionCDS, AR.effectAA, AR.effectCodone] = constructPriorStr(annotation);


figure
myhist(50, AR.f, 'r')
hold all
myhist(50, AR.f(AR.f>0.05), 'g')

AR = AR.filter('q', @(x)(x>7)); % mutant reads

AR.calcDxMin;
AR.visStat;
fig(gcf, 'width', 24)
exportfig(gcf, fullfile('figures', dataID,'qualityCtrl'), 'format','eps', 'color', 'rgb')

% visualizeAnnotationStat(annotation)


AR = AR.filter('r', @(x)(x<quantile(AR.r, 0.98))); % mutant reads
AR = AR.filter('q', @(x)(x>7)); % mutant reads

% AR = AR.filter('f', @(x)(x<1)); % SNP ratio
% AR = AR.filter('f', @(x)(x<.8)); % SNP ratio

 AR = AR.filter('dx', @(x)(x>120)); % SNP ratio
AR.unmix('plot');

AR.filterDx();

AR.plotChromosomes('dx', 'yscale', 'log', 'figure', 'new', 'plotfun', @(x,y)plot(x,y, 'rx'));
AR.plotChromosomes('dxFiltered', 'yscale', 'log', 'figure', 'old', 'plotfun', @(x,y)plot(x,y, 'm-'));

% AR = AR.filter('dxFiltered', @(x)(x>70491)); % SNP ratio
%% SNP ratio

AR.clearPlots;

KERNEL = 27;
AR.plotSnpRatio(KERNEL)
fig(gcf, 'width', 24)

% KERNEL = 81;
% AR.plotSnpRatio(KERNEL)
% fig(gcf, 'width', 24)
exportfig(gcf, fullfile('figures', dataID, sprintf('SNP_Ratio_k%u', KERNEL)), 'format','eps', 'color', 'rgb')

%% linkage map
load('./reference/ChrMap.mat')
AR.chrMap = ChrMap;
clear ChrMap;

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

% AR.includeRnaPresence()

AR.emissionHandle = @(q, r, study)emissionMixBetaBinomial(q, r, AR.pop, theta, lambda1);
% AR.emissionHandle = @(q, r, study)emissionMixBetaBinomial(q, r, AR.pop.N, theta, gi1(1:AR.Mtot));

%  emissionHandle = @(q, r, study)emissionBetaBinomial(q, r, study, theta);
% AR.emissionHandle = @(qq, rr, ff)emissionk0(qq, rr, AR.pop);
% AR.contrib = ones(size(AR.x));

% AR.Alpha = 1./(0:0.01:1);
AR.Alpha =  ALPHA;

%
% load( 'emission_fun.mat')
%  AR.E = E;
%  AR.contrib = contrib;
%  AR.contrib = 1;really
AR.calcEmission;
AR.run();

% AR.xPsel = calcMarginal(AR.xkPflat(:,floor(end/2):end),2);

AR.normalizeChromosomes;

AR.plotChromosomes('xPsel', 'yscale', 'lin', 'norm', true, 'figure', 'new');

AR.plotChromosomes('xPstat', 'yscale', 'lin', 'norm', true);

z = AR.calcMleZ('plot');
fig(gcf, 'width', 24)
exportfig(gcf, fullfile('figures', dataID, 'zMLE'), 'format','eps', 'color', 'rgb')

if exist('chr0', 'var')
    AR.plotChromosome2D(chr0)
end


% AR.calcLogOdds;
AR.plotChromosomes('xLogOdds', 'yscale', 'lin', 'figure', 'new', 'yThr', 0);

fig(gcf, 'width', 24)
exportfig(gcf, fullfile('figures', dataID,'LogLiOdds'), 'format','eps', 'color', 'rgb')
%%
AR.plotStemsLP( 'ylim', [-20,0])
if exist('chr0', 'var')
    ind = (AR.x==x0 & AR.chromosome == chr0);    
    AR.plotChromosomes('xPselNorm', 'yscale', 'lin', 'figure', 'old', 'plotfun',...
        @(x,y)plot(x,y, 'r^', 'markersize', 6, 'linewidth', 2), 'select', ind, 'ylim', [-20,0]);
end
fig(gcf, 'width', 24)
exportfig(gcf, fullfile('figures', dataID, 'LH-Posterior'), 'format','eps', 'color', 'rgb')

%%
AR.plotChromosomes('xPselNorm', 'yscale', 'lin', 'figure', 'new', 'yThr', 0);

AR.plotChromosomes('xPosteriorNorm', 'yscale',  'lin','figure', 'old', 'yThr', 0); %, 'ylim', [-10,0]);

% AR.plotChromosomes('contrib', 'yscale', 'lin', 'figure', 'new', 'yThr', 0);
%% print
numberOfHits = 2e6;
[~,dataName,ext] = fileparts(dataID);
AR.printTopHits(fullfile('figures', dataID, [dataName,ext, '-out.csv']), numberOfHits, 'xLogOdds', 'cutoffValue', -Inf)
%%
if exist('chr0', 'var')
    figure
    plot(AR.x(AR.chromosome == chr0)*1e-6, AR.dx(AR.chromosome == chr0), 'rx')
    set(gca, 'yscale', 'log')
end

%%

% AR.plotStems('xPosteriorNorm', 'figure', 'new');
% AR.plotStems('xPselNorm', 'yThr', 0);

%
% E = AR.E;
% contrib = AR.contrib;
% save('emission_obj', 'E', 'contrib')

finishing_beep(1, .2, 3)