function [ P ] = logProbSnpSpacingInMissMappedReads( t, T, N)
%LOGPROBSNPSPACINGINMISSMAPPEDREADS -- returns theoretical distribution of log10(inter-SNP distance)
%
% Assumes read length `t` is continuous (a real number)
%
% ## INPUTS:
% t  --  distance (independent variable)
% T  --  read length (parameter)
% N  --  number of SNPs per read (parameter)
%
% ## OUTPUT:
% P

t_rat = bsxfun(@rdivide, t, T);
P = log(10) * bsxfun(@times, M + 1, bsxfun(@times, t_rat, bsxfun(@power, (1 - t_rat), M)) );

%% set all entries with (t > T) to zero:

zeroInds = bsxfun(@and, bsxfun(@gt, t, T), N>=0);
P(zeroInds) = 0;

end

