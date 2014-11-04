function run_obj_hmm(data, ALPHA )

DATA_PATH = '/media/Processing/seq/data';
[dataID, N, chr0, x0] = data{:};

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
AR = AR.filter('f', @(x)(x<.8)); % SNP ratio

 AR = AR.filter('dx', @(x)(x>120)); % SNP ratio
AR.unmix('plot');

AR.filterDx();

AR.plotChromosomes('dx', 'yscale', 'log', 'figure', 'new', 'plotfun', @(x,y)plot(x,y, 'rx'));
AR.plotChromosomes('dxFiltered', 'yscale', 'log', 'figure', 'old', 'plotfun', @(x,y)plot(x,y, 'm-'));

% AR = AR.filter('dxFiltered', @(x)(x>70491)); % SNP ratio
%% SNP ratio

AR.clearPlots;

KERNEL = 81;
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