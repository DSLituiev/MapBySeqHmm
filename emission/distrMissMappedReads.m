function [ P ] = distrMissMappedReads( t, T, N)
%DISTRMISSMAPPEDREADS -- returns theoretical distribution of log10(inter-SNP distance)
% 
% INPUTS:
% t  --  distance (independent variable)
% T  --  read length (parameter)
% N  --  number of SNPs per read (parameter)


P = log(10) * bsxfun(@times,  N + 1, bsxfun(@times, t/T, bsxfun(@power, (1 - t./T), N)) );

%% set all entries with (t > T) to zero:
dimT = size(t)>1;
S.type = '()';
S.subs = repmat( {':'}, numel(size(P)), 1 );
S.subs{dimT} = (t > T) ;
P = subsasgn(P, S, 0);

% P( t./T > 1, :) = 0;
 
end

