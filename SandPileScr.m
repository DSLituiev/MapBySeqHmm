% close all; clear all;
% clc
% addpath('D:\MATLABuserfunctions\mtimesx');
% addpath('D:\MATLABuserfunctions\MinMaxSelection');
% addpath('D:\MATLABuserfunctions'); savepath;
% addpath('.\utilites');
% 
load('./data/HL7-AR.mat')

filt.cThr = .5;
filt.f_max = .65;
filt.r_max  = 1e3;
%== discard reads in repeat regions
  RR = AR;
% [ RR ] = selectFieldIndices( AR , ~AR.notaRepeat);
% [ RR ] = selectFieldIndices( AR , ...
%     AR.notaRepeat & AR.contrib> filt.cThr  &...
%     AR.f < filt.f_max & AR.r < filt.r_max  );
 

%% filter
filt.r_max = quantile(RR.r, filt.r_up_quantile);
inds = ( RR.qual>=10 &...
    RR.q >= filt.q_min &...
    RR.r < filt.r_max &...
    RR.f < filt.f_max );
%= r
RR.contrib(~inds) = 0.1;
RR.contrib( RR.contrib < filt.cThr) = 0;
%%
% clear AR;
markerSz = 4;
%= plot the SNP ratio
[f, spl] = plotAllChrNt(RR, 'f', 'exp10', false, 'ylim', [0 1],...
    'plotfun',@(x,y)plot(x,y,'MarkerEdgeColor','r', ...
    'MarkerSize', markerSz, 'Color', 'r'),...
    'FigName', 'SNP ratio', 'yscale', 'lin' );

%% general experimental constants:
%= number of plants:
study.N = 50;
%= number of potentially positive plants by a SNP:
study.kvect = uint32(0: 1: study.N)';
%= number of chromosomes:
study.chrnum = length(ChrMap);
%== cutoff for emission calculation
CUT_READ_NUM_BINOM = 120;
%% approximate binomial by Gaussian for high numbers
[EmMatrix, limReads] = initializeLimEmissionMatrix(study, RR, CUT_READ_NUM_BINOM);
%% flags
flag.prior = isfield(RR, 'logPrior'); %== use prior if available

%% initialization of the Calculation loop
M = zeros(study.chrnum, 1);  % number of reads on the chromosome
%= log-likelihood of the chromosome to be stationary
statChrLogProb = zeros(study.chrnum, 1);
%= log-likelihood of the chromosome to be under selection
slctChrLogProb = zeros(study.chrnum, 1);
%= log-likelihood of the locus to be under selection (within chromosome)
RR.logPoSlct = zeros(size(inds));
%=
Alpha = 5;
%%   Calculation cycle
for chr = chr0 % 1: study.chrnum;
    ticInit = tic;
    [ inds, M(chr), x, q, r,  contrib ] = recastChromosomeReads( RR, chr );
    
    [  RR.logPoSlct(inds), slctChrLogProb(chr), statChrLogProb(chr),...
        Es{chr}, Ts{chr}, taus{chr} ]  = ...
        estimateLikelihoodsOnChromosome( ChrMap(chr), EmMatrix, limReads, study.N,...
        x, q, r, contrib, Alpha);
    %== check for NaNs
    
    fprintf('\t The iterations for the chr.#\t%u\t took \t%4.2f\t s \t SNPs:\t%u \n', chr, toc(ticInit), sum(inds) )
    %== Ready!
end

if ~exist( 'chr0', 'var'); chr0 = 1; end
inds = RR.chromosome== chr0 ;

xMb = double(RR.x)*1e-6;

figure; surf(xMb(inds) , 0:50, Es{chr0}, 'linestyle', 'none')

figure; surf(  Es{chr0}, 'linestyle', 'none')
% figure; plot(xMb, sum( Es{chr0}) )

p0 = sum( inds ) *log10(  study.N+1 );

stat = StationaryDistr( study.N )'* Es{chr0};
logMmeanStat0 = log10(mean(stat));
RR.logPoSlctCorr( inds)  = RR.logPoSlct( inds) + log10(RR.contrib(inds));

figure
scatter( xMb(inds), RR.logPoSlct( inds ) + p0, 'g')
hold on
scatter( xMb(inds),  RR.logPoSlctCorr( inds)  + p0, 'b')
% plot( xMb(inds), RR.logPoSlctChr( inds) + log10(Es{chr}(end,:))' + p0, 'c-')
scatter( xMb(inds& RR.x == x0), RR.logPoSlct( inds & RR.x == x0 )+ p0, 'r')
plot(get(gca, 'xlim'), statChrLogProb(chr)*[1,1] + p0, 'g-')
plot(get(gca, 'xlim'), [0, 0], 'k-')
plot( xMb(inds), log10(Es{chr}(end,:)) , 'r:')
plot( get(gca, 'xlim'), [1, 1]*log10(mean(Es{chr}(end,:))), 'r--')
plot( get(gca, 'xlim'), [1, 1]*logMmeanStat0(1), 'g--')
plot( xMb(inds), stat, 'g:')
% x = RR.x( inds );
% logPoSlctChr = RR.logPoSlctChr( inds );
% weighted = RR.logPoSlctChr( inds) + log10(RR.contrib(inds));
% save('./data/HL7-result-soft-filtering.mat', 'x', 'logPoSlctChr', 'weighted')
 
% load('./data/HL7-result-soft-filtering.mat', 'x', 'logPoSlctChr', 'weighted')
% 
% figure
% scatter(x*1e-6, logPoSlctChr, 'g')
% hold on
% scatter(x*1e-6, weighted, 'b')
% plot(get(gca, 'xlim'), statChrLogProb(chr)*[1,1], 'g-')



% 
% %% calculate P(other chrs are stationary)
% I = logical(eye(study.chrnum));
% othersStatLogProb = -Inf(study.chrnum, 1);
% % thisSelLogProb = -Inf(study.chrnum, 1);
% thisSelLogProb0 = -Inf(study.chrnum, 1);
% totIntChrSum = -Inf(study.chrnum, 1);
% RR.logPoSlctNC = -Inf(size(RR.logPoSlctChr));
% RR.logPoSlct = -Inf(size(RR.logPoSlctChr));
% 
% for chr = study.chrnum:-1:1
%     othersStatLogProb(chr) = sum( statChrLogProb(~I(:,chr)) );
%     inds = (RR.chromosome == chr);
%     thisSelLogProb0(chr) = calcMarginal( RR.logPoSlctChr(inds) );
% end
% 
% % LOD = statChrLogProb - thisSelLogProb0;
% 
% totIntChrSum = calcMarginalDim( [statChrLogProb, thisSelLogProb0], 2 );
% % 
% % aSt = statChrLogProb - totIntChrSum;
% % aSl = thisSelLogProb0 - totIntChrSum;
% % 
% % for chr = study.chrnum:-1:1
% %     c(chr) = sum( aSt(~I(:,chr)) )+ sum( aSl(I(:,chr)) );
% %     d(chr) = log10( 1 - 10.^sum( aSt(I(:,chr)) ) ) + sum( log10( 1 - 10.^ aSl(~ I(:,chr)) ) );
% % end
% % 
% % calcMarginalDim( [a,b], 2)
% % [ 10.^a , 10.^b]
% 
% 
% % normFactChr = othersStatLogProb - calcMarginal(thisSelLogProb);
% % othersStatLogProb - calcMarginal(othersStatLogProb)
% 
% for chr = study.chrnum:-1:1
%     inds = (RR.chromosome == chr);
% %     RR.logPoSlctNC(inds) = RR.logPoSlctChr(inds) + normFactChr(chr);
%     RR.logPoSlct(inds) = RR.logPoSlctChr(inds) - totIntChrSum(chr);
% end
% 
% %== dot product:
% % thisChr = slctChrLogProb + othersStatLogProb;
% % normFactor = calcMarginal(thisChr);
% 
% normFactor = calcMarginal(RR.logPoSlct);
% RR.logPobsSelectionNorm  = RR.logPoSlct - normFactor;
% % thisSelLogProb(chr) = selChrLogProb(chr)
% 
% %% plot the likelihood
% markerSz = 4;
% plotcutoff = -10;
% %== not normalized !!!! (only comparable within each chromosome)
% 
% 
% [f, spl] = plotAllChrNt(RR, 'logPobsSelectionNorm','exp10',false,'ylim',[plotcutoff 0],...
%     'plotfun',@(x,y)stem(x,y,'MarkerEdgeColor','g',  'MarkerSize', markerSz, ...
%     'Color', 'g', 'BaseValue',(plotcutoff-2)),...
%     'FigName', 'Likelihood', 'yscale', 'lin');
% 
% if isfield(RR, 'logPrior')
%         z = calcMarginal(RR.logPoSlct + RR.logPrior);
%         RR.logPost  = RR.logPoSlct + RR.logPrior - z;
% %     RR.logPrior = multiplyAndNormalizeLogProbEachChr(  RR.logPrior, 0, RR.chromosome );
% %     RR.logPost  = multiplyAndNormalizeLogProbEachChr( RR.logPobsSelectionNorm, RR.logPrior, RR.chromosome, statChrLogProb );
%     RR.logPost  = RR.logPost  - calcMarginal(RR.logPost);
%     %     figure
%     %     plot(RR.logPost  -  RR.logPoSlct)
%     [fn, spln] = plotAllChrNt(RR, 'logPost','exp10',false,'ylim',[plotcutoff 0],...
%         'plotfun',@(x,y)stem(x,y,'MarkerEdgeColor','b',  'MarkerSize', markerSz,...
%         'Color', 'b', 'BaseValue',(plotcutoff-2)),...
%         'FigName', 'Likelihood', 'yscale', 'lin',...
%         'OldFig', f, 'OldSp', spl  );
%     
%     for ii = numel(spl):-1:1
%         v = get(spl(ii),'Children');
%         % v = allchild(spl(ii));
%         set(spl(ii),'Children',flipud(v))
%     end
%     
% end
% 
% axis(spl);
% 
% if exist( 'x0', 'var')
%     [f, spl] = plotAllChrNt(RR, 'logPost','exp10',false,'ylim',[plotcutoff 0],...
%         'plotfun',@(x,y)stem(x,y,'MarkerEdgeColor', 'r',  'MarkerSize', markerSz+2, ...
%         'Color', 'g', 'BaseValue',(plotcutoff-2)),...
%         'FigName', 'Likelihood', 'yscale', 'lin',...
%         'OldFig', f, 'OldSp', spl,...
%         'select', (RR.x == x0& RR.chromosome == chr0) );
% end
% % dbclear all; 