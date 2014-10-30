%=== script for general descriptive analysis of read data 


close all; clear all; clc;
tic

% runID = 'ABD159';
runID = 'HL7';
disp(['=======  Processing data from the run ''', runID, ''' ======='])
load ChrMap


dataPath = fullfile('data', [runID, '-ems-annotation-repfilt.csv'] );
[OR , ChrNumber] = readSequencingDataCsv(dataPath);
%%
OR = calcDx(OR);

%%
OR = unmixRepeats( OR, 'plot' );

%%

dxAll = OR.dx( ~isnan(OR.dx) );

dxNonRep = diff(double(OR.x(OR.notaRepeat)));
dxNonRep = dxNonRep(dxNonRep>0);

nBins = 25;

[fDxAll, xDx] = hist(log10(dxAll), nBins);
[fDxNonRep, xDx] = hist(log10(dxNonRep), xDx);


figure
b = barstairs(10.^xDx, fDxAll, 'b', 'EdgeColor','b');
alpha(b, 0.4)
hold all
b = barstairs(10.^xDx, fDxNonRep, 'r', 'EdgeColor','r');
alpha(b, 0.4)
set(gca, 'xscale', 'log')
%%
nBins = 25;
figure
[fAll, xf] = hist( log10(double(OR.r ) ) , nBins);
b = bar(xf, fAll, 'b');
hold all
%---
[fIncLow, xf] = hist( log10(double(OR.r( OR.notaRepeat | strcmpi( OR.repeat, 'Low_complexity')) ) ) , xf , 'facecolor', 'r');
b = bar(xf, fIncLow , 'g','EdgeColor','w');
%---
[fNonRep, xf] = hist( log10(double(OR.r( OR.notaRepeat) ) ) , xf , 'facecolor', 'r');
b = bar(xf, fNonRep ,'r','EdgeColor','w');

%%
figure
%--
% barstairs(xf, fIncLow./sum(fIncLow),  'g', 'EdgeColor','w')
%---
b = barstairs(xf, fNonRep./sum(fNonRep), 'r', 'EdgeColor','r');
 alpha(b, 0.4)
%---
b = barstairs(xf, fAll./sum(fAll), 'b');
 alpha(b, 0.4)

%%
figure
scatter(OR.q, OR.qAltQual )
hold all
scatter(OR.p, OR.pRefQual )

scatter(OR.r , OR.qual )
set(gca, 'xscale','log', 'yscale','log')

figure
scatter(OR.x(OR.chromosome ==1) , OR.r(OR.chromosome  ==1) , 'b')
hold all
scatter(OR.x(OR.chromosome ==1 & OR.notaRepeat) , OR.r(OR.chromosome  ==1& OR.notaRepeat) , 'r')



quals = 10.^(- OR.qual );

figure
hist((quals), 50)

pErrThr = 0.1;
QThr = -10*log10(pErrThr);

fi = fields(OR);
for ii = numel(fi):-1:1
NewObservationReads.(fi{ii})= OR.(fi{ii})(OR.qual> QThr);
end
%% = sums of alleles
figure
hist( log10(double(NewObservationReads.r) ) , nBins)

%% = individual alleles
allVariantCounts = [NewObservationReads.r - NewObservationReads.q; NewObservationReads.q] ;
allVariantCounts(allVariantCounts == 0)  = [];

figure
hist( log10(double(allVariantCounts) ) , nBins)

%%
NewObservationReads.p  = NewObservationReads.r-NewObservationReads.q;
pNonZeros = double(NewObservationReads.p(NewObservationReads.p >0));

parmhat = nbinfit(pNonZeros);
m = mean( ( pNonZeros ));
v = std( ( pNonZeros )  );

% mu = log(m^2/sqrt(v + m^2));
mu = mean( ( pNonZeros ) );
% sigma = sqrt( log(1 + v/m^2) );
sigma = std( ( pNonZeros )  ) ;
% 2*abs(mu - m); % std( log10( pNonZeros )  ) ;
x = round(10.^(0.1:0.1:5));
%  mu^2/(sigma^2 - mu)
% beta = mu/sigma;
% alpha = mu* beta;
 sigma = 15.2;
 r = parmhat(1); % mu^2/(sigma^2 - mu)
 p = parmhat(2); % mu/sigma 
 mode = (r-1)* p/(1-p)
% r = 4;
% g = pdf('gamma', (x), mu*beta ,  1/beta);
% g = pdf('logn', (x), mu ,  sigma);
%  g = pdf('nbin', (x), r ,  p);
t = .03; 
% g = t * pdf('poisson', (x), 10^0.85) + (1 - t) * pdf('poisson', (x), 10^1.88)  ;

mu1 =  10^0.85; 
mu2 = 10^2;
beta = .04;
g = t * pdf('poisson', (x), mu1) + (1 - t) * pdf('gamma', (x), mu2*beta, beta^-1)  ;
 
[yh, xh] = hist( log10( pNonZeros ), nBins );

figure
bar( xh , yh)
hold all
plot( log10(x), 2* g* sum(pNonZeros)/40, 'r')

% set(gca, 'xscale', 'log')

