function [E, c] = emissionBetaBinomial(q, r, study, theta)
% emissionk returns emission probabilities given the emission look-up matrix
%
% mins = [~, min_r, min_q]
%% P[q|r,k]
M = numel(q);
f = double(study.kvect)./(2*study.N);
% E = binopdf(qq, rr, ff);
E = 10.^ logBetaBinomialThetaMu0(double(q), double(r), f', theta);

cT = sum(E,2);
%% cF emission from 'false' states
fF =  double(study.N + study.kvect(2:end))./(2*study.N);
% REMOVE SUMMATION
EF = 10.^ logBetaBinomialThetaMu0( double(q), double(r), fF', theta);

cF = sum(EF,2);
cTotal = cT + cF;
c = cT./(cTotal);
% E(:, c>.5) = flipud(E(:, c>.5));
% c = ones(size(cT));
% figure; plot([cTotal; cT; cF]')
% figure; surf(E, 'linestyle', 'none')
% figure; hist( log10(1-c(1-c>0)), 40 )

% E(:, sum(E,1) == 0 ) = 1/size(E,1);
 E = bsxfun(@rdivide, E, cTotal);