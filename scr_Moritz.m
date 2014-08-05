close all; clear all; clc;
dbclear if warning
tic

DATA_PATH = '/media/Processing/seq/data';
% DATA_PATH = './data';
USERFNCT_PATH = '/media/Processing/';
addpath(fullfile(USERFNCT_PATH, 'MATLABuserfunctions/binomial') );
addpath(fullfile(USERFNCT_PATH, 'MATLABuserfunctions/mtimesx'));
addpath(fullfile(USERFNCT_PATH, 'MATLABuserfunctions/MinMaxSelection'));
addpath(fullfile(USERFNCT_PATH, 'MATLABuserfunctions/newtonraphson'));
addpath(fullfile(USERFNCT_PATH, 'MATLABuserfunctions')); savepath;
addpath('./utilites');
addpath('./betabinomial');
% addpath('..\@chromProb');

%=       provide the reference ID, which is the name of the csv file without '.csv'
%= together with the backround genotype ID
%= known positions of the causative SNP can be also provided here for
%= further visualization
% dataID = 'ABD159-rmdup-clipOverlap-q20-freebayes'; chr0 = 2; x0 = 17521246; % bkgrID = 'ABD241-rmdup-clipOverlap-freebayes';
% dataID = 'ABD173-rmdup-clipOverlap-freebayes'; chr0 = 3 ; x0 = 1619248; % bkgrID = 'ABD241-rmdup-clipOverlap-q20-freebayes';


% dataID = 'HL10-rmdup-clipOverlap-q20-freebayes';  chr0 =  3 ;  x0 =  16473265;
 dataID = 'MO7_mutant_pool-rmdup-clipOverlap-q20-freebayes-ems-annotation-repfilt';
 
% dataPath = '/media/Processing/seq/data/MO8_mutant_pool-rmdup-clipOverlap-q20-freebayes-ems-annotation-repfilt.csv';

% dataPath = '/media/Processing/seq/data/MO8_WT_pool-rmdup-clipOverlap-q20-freebayes-ems-annotation-repfilt.csv';
%  dataID = 'HL7_Paired-rmdup-clipOverlap-q20-freebayes'; x0 = 5672441; chr0 = 1;
% dataID = 'HL7_Paired-ngm-rmdup-clipOverlap-freebayes'; x0 = 5672441; chr0 = 1;
% dataID = 'HL7_Paired-rmdup-clipOverlap-mpileup'; x0 = 5672441; chr0 = 1;

disp(['=======  Processing data from the run ''', dataID, ''' ======='])

%%
% mkdir(fullfile('figures',dataID))
%=       load the recombination map (contains positions of the markers
%= and genetic distance in cM between them)

%= construct the path to the (primary experimental) data file
dataPath = fullfile(DATA_PATH, [dataID, '.csv'] );

%= extract the refenece reads if the reference ID is given:
clear AR1 AR z y
% [AR, annotation] = subtractBackGroundGenotype(dataPath);
 [AR, ~] = readSequencingDataCsv(dataPath, 'noannotation');
%% general experimental constants:
%= number of plants:
N = 50;
AR.pop = N;
%= and genetic distance in cM between them)
load ChrMap
AR.chrMap = ChrMap;
clear ChrMap;


AR.filter('q', @(x)(x>7))

AR.filter('f', @(x)(x>.01))

AR.visStat;

AR = calcDxMin(AR);

% AR = unmix(AR);



AR.Mtot
 theta = 2e-2;
 lambda1 = 1;

AR.emissionHandle = @(q, r, study)emissionMixBetaBinomial(q, r, AR.pop, theta, lambda1);

AR.Alpha = 17;

AR.calcEmission;

AR.run();

AR.clearPlots;

AR.normalizeChromosomes;

AR.plotChromosomes('xPsel', 'yscale', 'lin', 'norm', true, 'figure', 'new');

AR.plotChromosomes('xPstat', 'yscale', 'lin', 'norm', true);

% AR.calcLogOdds;
AR.plotChromosomes('xLogOdds', 'yscale', 'lin', 'figure', 'new', 'yThr', 0);


AR.plotChromosomes('f', 'yscale', 'lin', 'figure', 'new', 'yThr', 0.25);




return

figure;
myhist100(50, AR.f,  'r')
hold all
myhist100(50, AR.f(AR.q>5),  'b')
myhist100(50, AR.f(AR.q>7),  'c')
myhist100(50, AR.f(AR.q>10),  'g')


figure;
myhist(50, AR.f,  'r')
hold all
myhist(50, AR.f(AR.q>5),  'b')
myhist(50, AR.f(AR.q>7),  'c')
myhist(50, AR.f(AR.q>10),  'g')
set(gca, 'yscale', 'log')

figure;
myhist100(50, AR.f,  'r')
hold all
myhist100(50, AR.f(AR.f>.05),  'g')

figure;
hist(log10(AR.r), 50)


figure;
hist(AR.f(AR.q>10& AR.r <=100), 50)

for chr =5:-1:1
figure('name', sprintf('chr %u', chr));
plot(AR.x(AR.qual>.0 & AR.q>3 & AR.f<0.75 & AR.chromosome == chr), AR.f(AR.qual>.0 & AR.q>3 & AR.f<0.75 & AR.chromosome == chr))
hold on
plot([min(AR.x), max(AR.x)], .5*[1,1], 'r-')
ylim([0,1])
end
