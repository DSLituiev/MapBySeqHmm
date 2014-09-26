function [ sumLogP ] = calcMarginal( logP, varargin )
%calcMarginal calculates log-marginal given log-probabilities as input
% taking care of numeric overflow
%
%   two modes of operation:
%  (1) all values in similar range ('acceptableRange' set to 7)
%      takes the maximal value out of parantheses
%
%  (2) the values vary tremendously
%        the following trick is used:
%
% 10^(log(a)-log(a)) + 10^(log(b)-log(a)) = (a+b)/a
% a*(a+b)/a = a*(1+b/a) = a*( 1 + 10^(log(b)-log(a)) )
% log(a+b) = log(a) + log(1+ 10^(log(b)-log(a)) )
%
% To see how it works, try:
%    a = 5e-4; b = 12e-4; a+b
%    10^(log10(a) + log10(1+ 10^(log10(b) - log10(a) )))

acceptableRange = 7;

if nargin>1 && ~ logical(varargin{1} - floor(varargin{1}) )
    dim = varargin{1};
else
    dim = 1;
end

if all( abs( max(logP, [], dim) - min(logP, [], dim) ) < acceptableRange )
    maxLogP = max(logP, [], dim);
    scaledLogP = bsxfun(@minus, logP, maxLogP);
    sumLogP = log10(sum(10.^scaledLogP, dim)) + maxLogP;
else
    maxLogP = max(logP, [], dim);
    scaledLogP = bsxfun(@minus, sort(logP, dim, 'descend'), maxLogP );
    
    S.type = '()';     S.subs = cell(1, ndims( logP ));
    S.subs(:) = {':'};
    S.subs{dim} = 1;    
    sumLogP = subsref(scaledLogP, S);

    for ii = 2:size(scaledLogP, dim)
        if ~isinf(scaledLogP(ii))
            S.subs{dim} = ii;
            logDiff = bsxfun(@minus, sumLogP, subsref(scaledLogP, S) );
            sumLogP = sumLogP +...
                log10(1 + 10.^( -logDiff ));
%         else
%             warning('logP = -Inf!');
        end
        if any(isinf(sumLogP))
           warning('logP = -Inf!');
        end
    end
   sumLogP = sumLogP + maxLogP;
end

if any(isnan(sumLogP))
     warning('calcMarginal:NaN', '%u entries out of %u are NaN', sum(isnan(sumLogP(:))), numel(sumLogP));
%    warning('calcMarginal:NaN', 'the sum is NaN in the entries %u', find(isnan(sumLogP))) 
end

end

