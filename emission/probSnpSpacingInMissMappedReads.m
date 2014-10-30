function [ P ] = distrMissMappedReads( t, T, N)
%DISTRMISSMAPPEDREADS -- returns theoretical distribution of log10(inter-SNP distance)
%
% INPUTS:
% t  --  distance (independent variable)
% T  --  read length (parameter)
% N  --  number of SNPs per read (parameter)

P = LikelihoodMissMappedReads(t, T, N);

normFactor = zeros(size(T));
for ii = numel(T):-1:1
    normFactor(ii) = sum(LikelihoodMissMappedReads(1:T(ii), T(ii), N));
end

P = bsxfun(@rdivide, P, normFactor);
%% set all entries with (t > T) to zero:
% dim_t = (size(t)>1);
% S.type = '()';
% S.subs = repmat( {':'}, numel(size(P)), 1 );
% S.subs{dim_t} = bsxfun(@gt, t, T);
% P = subsasgn(P, S, 0);
zeroInds = bsxfun(@and, bsxfun(@gt, t, T), N>=0);
P(zeroInds) = 0;
% P( t./T > 1, :) = 0;

    function LL = LikelihoodMissMappedReads(t, T, N)
        t_rat = bsxfun(@rdivide, t, T);
        LL = bsxfun(@times,  bsxfun(@power, (1 - t_rat), N));
    end

end

