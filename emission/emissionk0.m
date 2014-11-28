function [E, c] = emissionk0(q, r, study)
% emissionk returns emission probabilities
% mins = [~, min_r, min_q]
%% P[q|r,k]
assert(numel(q) == numel(r), 'the number of elements in `q` and `r` differs!')
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