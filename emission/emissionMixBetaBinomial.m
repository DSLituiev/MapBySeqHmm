function [E, c] = emissionMixBetaBinomial(q, r, study, theta, lambda1)
% emissionk returns emission BetaBinomial-Uniform mixture probabilities
%

if ~isreal(theta)|| ~isreal(lambda1)
    error('emissionMixBetaBinomial:complexParams', 'complex theta or lambda!')
end
%% P[q|r,k]
% M = numel(q);
f = double(study.kvect)./(2*study.N);
% E = binopdf(qq, rr, ff);
% E = 10.^ logBetaBinomialThetaMu0(double(q), double(r), f', theta);

N = study.N;

[pBB, pU] = conditBetaBinomStat(q, r, theta, N, f, 1./N);

gi1 = lambda1*pBB./(lambda1*pBB + (1-lambda1)*pU );
 
% E = (1-lambda1)/study.N+ lambda1*10.^logBetaBinomialThetaMu0(q, r, f, theta);
E = bsxfun(@plus, (1-gi1)./1 , bsxfun(@times, gi1, 10.^logBetaBinomialThetaMu0(q,r,f,theta)));
c = 1;
% E(:, c>.5) = flipud(E(:, c>.5));
% c = ones(size(cT));
% figure; plot([cTotal; cT; cF]')
% figure; surf(E, 'linestyle', 'none')
% figure; hist( log10(1-c(1-c>0)), 40 )

% E(:, sum(E,1) == 0 ) = 1/size(E,1);

%   E = bsxfun(@rdivide, E, sum(E,2));