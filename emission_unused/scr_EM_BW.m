close all; clear all; clc;
dbclear if warning
tic

USERFNCT_PATH = '/media/Processing/';
% addpath(fullfile(USERFNCT_PATH, 'MATLABuserfunctions/binomial') );
addpath(fullfile(USERFNCT_PATH, 'MATLABuserfunctions/mtimesx'));
addpath(fullfile(USERFNCT_PATH, 'MATLABuserfunctions/MinMaxSelection'));
addpath(fullfile(USERFNCT_PATH, 'MATLABuserfunctions/newtonraphson'));
addpath(fullfile(USERFNCT_PATH, 'MATLABuserfunctions')); savepath;
addpath('../utilites');
addpath('../betabinomial');
% addpath('..\@chromProb');
addpath('..');

%=       provide the reference ID, which is the name of the csv file without '.csv'
%= together with the backround genotype ID
%= known positions of the causative SNP can be also provided here for
%= further visualization
% dataID = 'ABD159-rmdup-clipOverlap-q20-freebayes'; chr0 = 2; x0 = 17521246; % bkgrID = 'ABD241-rmdup-clipOverlap-freebayes';
% dataID = 'ABD173-rmdup-clipOverlap-freebayes'; chr0 = 3 ; x0 = 1619248; % bkgrID = 'ABD241-rmdup-clipOverlap-q20-freebayes';


% dataID = 'HL10-rmdup-clipOverlap-q20-freebayes';  chr0 =  3 ;  x0 =  16473265;

 dataID = 'HL7_Paired-rmdup-clipOverlap-freebayes'; x0 = 5672441; chr0 = 1;
%  dataID = 'HL7_Paired-rmdup-clipOverlap-q20-freebayes'; x0 = 5672441; chr0 = 1;
% dataID = 'HL7_Paired-ngm-rmdup-clipOverlap-freebayes'; x0 = 5672441; chr0 = 1;
% dataID = 'HL7_Paired-rmdup-clipOverlap-mpileup'; x0 = 5672441; chr0 = 1;

disp(['=======  Processing data from the run ''', dataID, ''' ======='])
%= and genetic distance in cM between them)

%%
% mkdir(fullfile('figures',dataID))
%=       load the recombination map (contains positions of the markers
%= and genetic distance in cM between them)

%= construct the path to the (primary experimental) data file
dataPath = fullfile('../data', [dataID, '-ems-annotation-repfilt.csv'] );
%= extract the refenece reads if the reference ID is given:
clear AR1 AR 
[AR, annotation] = subtractBackGroundGenotype(dataPath);

AR = calcDxMin(AR);

AR = unmix(AR);

load ChrMap
AR.chrMap = ChrMap;
clear ChrMap;
%% general experimental constants:
%= number of plants:
N = 50;
AR.pop = N;
%%

[theta, lambda1, ~] = runSimpleEM_BetaBinomAndUniform(double(AR.q),...
    double(AR.r), N, 'v',...
    'contribution', AR.contrib ,'errTol', 1e-6);% , 'contribution', AR.contrib);
 
emissionHandle = @(q, r, study)emissionMixBetaBinomial(q, r, AR.pop, theta, lambda1);
%  emissionHandle = @(q, r, study)emissionBetaBinomial(q, r, study, theta);
% emissionHandle = @(qq, rr, ff)emissionk0(qq, rr, study);

% AR.contrib = ones(size(AR.x));
[AR.E, AR.c] = wrapEmissionMatrix(AR.q, AR.r, AR.pop, emissionHandle, AR.contrib);
Alpha = 17;

chr = 1;
AR.calcT(chr, Alpha);
AR.run(1, Alpha, 'plot')


reads = chromContainer(dataID, study);

for chr = AR.chrNumber:-1:1
  
    reads.addChromosome(getChromosome(AR, chr), chr)
    
    %= calculate Transition matrices
    reads.chromosome{chr}.calcT(study, ChrMap(chr), Alpha);
    
%   [theta(chr), lambda(chr), reads.chromosome{chr}.c] = runSimpleEM_BetaBinomAndUniform_obj(reads.chromosome{chr}, study);
    
    reads.chromosome{chr}.calcEmission( study, theta, lambda1);
    
%     [~, fh] = reads.chromosome{chr}.run('plot');
      [~, fh] = reads.chromosome{chr}.run();

end
%%
reads.plotChromosomes('Psel', 'exp10', false, 'yscale', 'lin', 'ylim' , 'loose')

reads.plotChromosomes('contrib',  'yscale', 'lin', 'ylim' , [0,1])

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
