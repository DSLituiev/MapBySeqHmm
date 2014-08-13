function [E, c] = emissionk0(q, r, study)
% emissionk returns emission probabilities given the emission look-up matrix
%
% mins = [~, min_r, min_q]
%% P[q|r,k]
M = numel(q);
qq = double(repmat(q,[ 1, study.N+1]));
rr = double(repmat(r,[ 1, study.N+1]));
ff = repmat(double(study.kvect)./(2*study.N), [ M, 1]);
E = binopdf(qq, rr, ff);
cT = sum(E,2);
%% cF emission from 'false' states
ffF = repmat( double(study.N + study.kvect(2:end))./(2*study.N), [M, 1]);
EF = binopdf(qq(:,2:end), rr(:,2:end), ffF);
cF = sum(EF, 2);
cTotal = cT + cF;
c = (cT./(cTotal));
% E(:, c>.5) = flipud(E(:, c>.5));
% c = ones(size(cT));
% figure; plot([cTotal; cT; cF]')
% figure; surf(E, 'linestyle', 'none')
% figure; hist( log10(1-c(1-c>0)), 40 )

% E(:, sum(E,1) == 0 ) = 1/size(E,1);
%  E = bsxfun(@rdivide, E, cTotal);