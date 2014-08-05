close all; clear all; clc;
dbclear if warning
tic


DATA_PATH = '/media/Processing/seq/data';
% DATA_PATH = './data';
USERFNCT_PATH = '/media/Processing/MATLABuserfunctions';
% addpath(fullfile(USERFNCT_PATH, 'MATLABuserfunctions/binomial') );
addpath(USERFNCT_PATH);
addpath(fullfile(USERFNCT_PATH, 'mtimesx'));
addpath(fullfile(USERFNCT_PATH, 'MinMaxSelection'));
addpath(fullfile(USERFNCT_PATH, 'newtonraphson'));
 savepath;
addpath('./utilites');
addpath('./betabinomial');
addpath('./plotting');
addpath('./emission');
% addpath('.\@chromProb');
% addpath('.');

%=       provide the reference ID, which is the name of the csv file without '.csv'
%= together with the backround genotype ID
%= known positions of the causative SNP can be also provided here for
%= further visualization

%= number of plants:

% dataID = 'MO8-rmdup-clipOverlap-q20-freebayes-ems-annotation'; N = 278;
dataID = 'MO7-rmdup-clipOverlap-q20-freebayes-ems-annotation'; N = 181;
% dataID = '73-ngm-rmdup-clipOverlap-q20';
% dataID = '39-ngm-rmdup-clipOverlap-q20-ems-annotation';
% dataID = 'ABD159-rmdup-clipOverlap-q20-freebayes-ems-annotation-repfilt'; chr0 = 2; x0 = 17521246; % bkgrID = 'ABD241-rmdup-clipOverlap-freebayes';
% dataID = 'ABD173-rmdup-clipOverlap-q20-freebayes-ems-annotation-repfilt'; chr0 = 3 ; x0 = 1619248; % bkgrID = 'ABD241-rmdup-clipOverlap-q20-freebayes';

% dataID = 'MO8_mutant_pool-rmdup-clipOverlap-q20-freebayes-ems-annotation-repfilt';

% dataID = 'HL10-rmdup-clipOverlap-q20-freebayes-ems-annotation-repfilt';  chr0 =  3 ;  x0 =  16473265;

%  dataID = 'HL7_Paired-rmdup-clipOverlap-freebayes-ems-annotation-repfilt'; x0 = 5672441; chr0 = 1;
%  dataID = 'HL7_Paired-rmdup-clipOverlap-q20-freebayes-ems-annotation-repfilt'; x0 = 5672441; chr0 = 1;
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

% AR.unmix();

% AR.unmix('plot');

load('./reference/ChrMap.mat')
AR.chrMap = ChrMap;
clear ChrMap;
%% general experimental constants:

AR.pop = N;
%%
% theta = 1e-2;
% lambda1 = 1;
% 
[theta, lambda1, ~] = runSimpleEM_BetaBinomAndUniform(AR.pop, double(AR.q),...
    double(AR.r), N, 'v',...
    'contribution', AR.contrib ,'errTol', 1e-4);% , 'contribution', AR.contrib);

AR.emissionHandle = @(q, r, study)emissionMixBetaBinomial(q, r, AR.pop, theta, lambda1);

%  emissionHandle = @(q, r, study)emissionBetaBinomial(q, r, study, theta);
% emissionHandle = @(qq, rr, ff)emissionk0(qq, rr, study);
% AR.contrib = ones(size(AR.x));

% AR.Alpha = 1./(0:0.05:1);
AR.Alpha = 17;
% 
% load( 'emission_fun.mat')
%  AR.E = E;
%  AR.contrib = contrib;
%  AR.contrib = 1;
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

AR.plotStemsLP( 'ylim', [-10,0])
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

AR.plotChromosomes('dx', 'yscale', 'log', 'figure', 'new',  'plotfun', @(x,y)plot(x,y, 'rx'));

figure
plot(AR.x(AR.chromosome == 5)*1e-6, AR.dx(AR.chromosome == 5), 'rx')
set(gca, 'yscale', 'log')

figure
plot(AR.x(AR.chromosome == 1)*1e-6, AR.dx(AR.chromosome == 1), 'rx')
set(gca, 'yscale', 'log')
%%
%%

% AR.plotStems('xPosteriorNorm', 'figure', 'new');
% AR.plotStems('xPselNorm', 'yThr', 0);
numberOfHits = 2000;

AR.printTopHits(fullfile('figures',dataID, [dataID, '-out.csv']), numberOfHits)

% 
% E = AR.E;
% contrib = AR.contrib;
% save('emission_obj', 'E', 'contrib')

return
% ff =  0:(1/2/N):(1);
q1 = 25;
q2 = 60;
figure
plot(f, (binopdf(25,100,f)) , 'g-' , 'linewidth', 2)
hold all
plot(f, 10.^logBetaBinomialThetaMu0(q1 ,100,f, theta), 'm-', 'linewidth', 2)
plot(f, (1-lambda1)*1/N + lambda1*10.^logBetaBinomialThetaMu0(q1 ,100,f, theta), 'r-', 'linewidth', 2)
plot(f, (binopdf(q2 ,100,f)) , 'g:', 'linewidth', 2)
plot(f, 10.^logBetaBinomialThetaMu0(q2 ,100,f, theta) , 'm:', 'linewidth', 2)
plot(f, (1-lambda1)*1/N + lambda1*10.^logBetaBinomialThetaMu0(q2,100,f, theta), 'r:', 'linewidth', 2)
legend({'binomial', '\beta-binomial', '\beta-binomial + uniform'})
fig(gcf, 'height', 12, 'width', 16)

% exportfig(gcf, fullfile('../figures', ['demo-beta-mix']), 'format','eps', 'color', 'rgb')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%% reading specific line from a text file  
% 
% linenum = 3; % Or whichever line you wish to read
% textscan(fid, '%s', 1, 'delimiter', '\n', 'headerlines', linenum-1);
