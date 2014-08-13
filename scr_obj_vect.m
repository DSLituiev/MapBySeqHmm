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
 savepath;
addpath('./utilites');
% addpath('./betabinomial');
addpath('./plotting');
addpath('./emission');
% addpath('.\@chromProb');
% addpath('.');

%=       provide the reference ID, which is the name of the csv file without '.csv'
%= together with the backround genotype ID
%= known positions of the causative SNP can be also provided here for
%= further visualization

%= number of plants:

dataID = 'MO8-rmdup-clipOverlap-q20-freebayes-ems-annotation-rna'; N = 278;  chr0 = 5;
% dataID = 'MO7-rmdup-clipOverlap-q20-freebayes-ems-annotation-rna'; N = 181; chr0 = 1;
% dataID = '73-ngm-rmdup-clipOverlap-q20';
% dataID = '39-ngm-rmdup-clipOverlap-q20-ems-annotation';
% dataID = 'ABD159-rmdup-clipOverlap-q20-freebayes-ems-annotation'; chr0 = 2; x0 = 17521246;  N = 100; % bkgrID = 'ABD241-rmdup-clipOverlap-freebayes';
% dataID = 'ABD173-rmdup-clipOverlap-q20-freebayes-ems-annotation'; chr0 = 3 ; x0 = 1619248;  N = 50;  % bkgrID = 'ABD241-rmdup-clipOverlap-q20-freebayes';

% dataID = 'MO8_mutant_pool-rmdup-clipOverlap-q20-freebayes-ems-annotation-repfilt';

% dataID = 'HL10-rmdup-clipOverlap-q20-freebayes-ems-annotation';  chr0 =  3 ;  x0 =  16473265;N = 50;

%  dataID = 'HL7_Paired-rmdup-clipOverlap-q20-freebayes-ems-annotation'; x0 = 5672441; chr0 = 1; N = 50;
% dataID = 'HL7_Paired-ngm-rmdup-clipOverlap-freebayes-ems-annotation-repfilt'; x0 = 5672441; chr0 = 1;
% dataID = 'HL7_Paired-rmdup-clipOverlap-mpileup-ems-annotation-repfilt'; x0 = 5672441; chr0 = 1;

disp(['=======  Processing data from the run ''', dataID, ''' ======='])
%= and genetic distance in cM between them)

%%

mkdir(fullfile('figures',dataID))
%=       load the recombination map (contains positions of the markers
%= and genetic distance in cM between them)

%= construct the path to the (primary experimental) data file
dataPath = fullfile(DATA_PATH, [dataID, '.csv'] );
%= extract the refenece reads if the reference ID is given:
clear AR 

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
AR = AR.filter('f', @(x)(x<1)); % SNP ratio

AR.calcDxMin;
AR.visStat;
fig(gcf, 'width', 24)
exportfig(gcf, fullfile('figures', dataID,'qualityCtrl'), 'format','eps', 'color', 'rgb')

% visualizeAnnotationStat(annotation)



% AR.unmix('plot');

load('./reference/ChrMap.mat')
AR.chrMap = ChrMap;
clear ChrMap;
%% general experimental constants:

AR.pop = N;
%%
if ~AR.flagWT
    AR.unmix();
    [theta, lambda1] = runSimpleEM_BetaBinomAndUniform(AR.pop.Pstat, double(AR.q),...
        double(AR.r), N, 'v',...
        'contribution', AR.contrib ,'errTol', 1e-5, 'contribution', AR.contrib);
else
    [theta, lambda1] = runSimpleEM_BetaBinomAndUniform(AR.pop.Pflat, double([AR.q; AR.qw]),...
        double([AR.r; AR.rw]), N, 'v',...
        'contribution', AR.contrib ,'errTol', 1e-4);
end
% 
% theta = 1e-2;
% lambda1 = 1;

% AR.includeRnaPresence()

AR.emissionHandle = @(q, r, study)emissionMixBetaBinomial(q, r, AR.pop, theta, lambda1);

%  emissionHandle = @(q, r, study)emissionBetaBinomial(q, r, study, theta);
% AR.emissionHandle = @(qq, rr, ff)emissionk0(qq, rr, AR.pop);
% AR.contrib = ones(size(AR.x));

% AR.Alpha = 1./(0:0.05:1);
AR.Alpha = 17;
% 
% load( 'emission_fun.mat')
%  AR.E = E;
%  AR.contrib = contrib;
%  AR.contrib = 1;really
AR.calcEmission;
AR.run();

AR.clearPlots;

AR.normalizeChromosomes;

AR.plotChromosomes('xPsel', 'yscale', 'lin', 'norm', true, 'figure', 'new');

AR.plotChromosomes('xPstat', 'yscale', 'lin', 'norm', true);



% AR.calcLogOdds;
AR.plotChromosomes('xLogOdds', 'yscale', 'lin', 'figure', 'new', 'yThr', 0);

fig(gcf, 'width', 24)
exportfig(gcf, fullfile('figures', dataID,'LogLiOdds'), 'format','eps', 'color', 'rgb')

AR.plotStemsLP( 'ylim', [-20,0])
fig(gcf, 'width', 24)
exportfig(gcf, fullfile('figures', dataID, 'LH-Posterior'), 'format','eps', 'color', 'rgb')



AR.plotChromosomes('xPselNorm', 'yscale', 'lin', 'figure', 'new', 'yThr', 0);

AR.plotChromosomes('xPosteriorNorm', 'yscale',  'lin','figure', 'old', 'yThr', 0); %, 'ylim', [-10,0]);


AR.plotChromosomes('f', 'yscale', 'lin', 'figure', 'new', 'yThr', 0, 'plotfun', @(x,y)plot(x,y, 'rx'));
if AR.flagWT
    AR.plotChromosomes('fw', 'yscale', 'lin', 'figure', 'old', 'yThr', 0, 'plotfun', @(x,y)plot(x,y, 'gx'));
end
set(AR.prevAxes(:,end), 'box','on')
fig(gcf, 'width', 24)
exportfig(gcf, fullfile('figures', dataID, 'SNP Ratio'), 'format','eps', 'color', 'rgb')

% AR.plotChromosomes('contrib', 'yscale', 'lin', 'figure', 'new', 'yThr', 0);


figure
plot(AR.x(AR.chromosome == chr0)*1e-6, AR.dx(AR.chromosome == chr0), 'rx')
set(gca, 'yscale', 'log')

%%
%%

% AR.plotStems('xPosteriorNorm', 'figure', 'new');
% AR.plotStems('xPselNorm', 'yThr', 0);
numberOfHits = 2e6;

AR.printTopHits(fullfile('figures',dataID, [dataID, '-out.csv']), numberOfHits, 'xLogOdds', 'cutoffValue', 0)

% 
% E = AR.E;
% contrib = AR.contrib;
% save('emission_obj', 'E', 'contrib')

finishing_beep(1, .2, 3)
