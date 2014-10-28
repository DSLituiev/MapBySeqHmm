function [ P ] = probSnpLogSpacingInMissMappedReads( t, L, M, varargin)
%PROBSNPLOGSPACINGINMISSMAPPEDREADS-- returns a theoretical distribution of log10(inter-SNP distance)
%
% Assumes read length `t` is continuous (a real number)
%
% ## INPUTS:
% t  --  distance                 (independent variable)
% L  --  read length              (parameter)
% N  --  number of SNPs per read  (parameter)
% [prior] -- 'tent' or'uniform'   (prior type)
% 
% ## OUTPUT:
% P[\log_{10} t |M, L, $ ], where $ is event of absence of any SNP within `t`
%%
p = inputParser;
addRequired(p, 't', @isnumeric);
addRequired(p, 'L', @isnumeric);
addRequired(p, 'M', @isnumeric);
%
addOptional(p,     'prior',  'tent',  @ischar);

parse(p, t, L, M, varargin{:});
%%
t_rat = bsxfun(@rdivide, t, L);
Likelihood = bsxfun(@power, (1 - t_rat), M);

% P = log(10) * bsxfun(@times, M + 1, bsxfun(@times, t_rat, bsxfun(@power, (1 - t_rat), M)) );

switch p.Results.prior
    case {'tent', 't'}
        Prior = tent( t, L );
        normFactor = 2*(1-2.^(-M-1))./ ((1 + M)*(2 + M));
    case {'uniform', 'u', 'flat', 'f'}
        Prior = 1./L;
        normFactor = 1 ./ (1 + M);
    otherwise
         error('probSnpLogSpacingInMissMappedReads:unknownPrior', 'unknown prior type')      
            
end

P = log(10) * bsxfun(@times, bsxfun(@times, t, Likelihood), bsxfun(@rdivide, Prior, normFactor));
%% set all entries with (t > T) to zero:

zeroInds = bsxfun(@and, bsxfun(@gt, t, L), M>=0);
P(zeroInds) = 0;

end

