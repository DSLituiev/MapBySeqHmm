function [sumdLdTheta, varargout] = maxLlThetaBetaBinomial(k, n, mu, theta, varargin)
% maxLlThetaBetaBinomial calculates log-likelihood (10-based) of
% Beta-Binomial distribution with given parameters mu and theta
% mu = alpha/(alpha + beta)
% theta = 1/(alpha + beta)
%
% References:
% Smith, D.M. "Algorithm AS 189: Maximum Likelihood Estimation  of the
% Parameters  of the Beta Binomial Distribution". Journal of the Royal
% Statistical Society. Series C (Applied Statistics),  Vol. 32, No. 2
% (1983), pp. 196-204;  URL: http://www.jstor.org/stable/2347299

if nargin>4 && ~isempty(varargin{1})
    flagLL = true;
else
    flagLL = false;
end

dimK = find(size(k)>1);
if isempty(dimK)
    dimK = find(size(theta)==1, 1, 'first');
end

dims = size(bsxfun(@plus, theta, k));
S.type = '()';
S.subs = cellstr(repmat(':', [numel(dims),1]));
if size(k)~=size(n)
    error('maxLlThetaBetaBinomial:k_and_n_mismatch', 'k and n must have the same size')
end


%% shortcut
% Smith, 1983
rr = permute(0:(max(n)-1), [1,3,2]);

S_k = cumS(k, rr, dimK);
S_n_k = cumS(n-k, rr,dimK);
S_n = cumS(n, rr, dimK);

rrTheta = bsxfun(@times, rr, theta);

rrTmu  = bsxfun(@plus ,     mu, rrTheta) ;
rrT1mu = bsxfun(@plus , 1 - mu, rrTheta);
rrT1   = 1 + rrTheta;
dLdTheta = sum( bsxfun(@times, rr, bsxfun(@minus,...
     bsxfun(@rdivide,  S_k ,  rrTmu    ) + ...
    bsxfun(@rdivide, S_n_k,  rrT1mu   ), ...
    bsxfun(@rdivide,  S_n ,  rrT1     ))  ), 3 );

if flagLL    
    LL =  bsxfun(@plus ,...
        sum( ...
        bsxfun(@times, S_k  , log10(rrTmu ) ) + ...
        bsxfun(@times, S_n_k, log10(rrT1mu) ) - ...
        bsxfun(@times, S_n  , log10(rrT1) ) ,...
        3 ),...
        log10(exp(1))* sum(gammaln(n + 1)-gammaln(k + 1)-gammaln(n - k + 1) , dimK ) ) ;
    %   sum( log10(binomial(n,k)), dimK ) ) ;
end



sumdLdTheta = sum(dLdTheta, dimK);

varargout{1} = dLdTheta;
if flagLL
    varargout{2} = LL;
end

%% subfunction
    function S_x = cumS(x, rr, dimK)
        F_x = uint32(bsxfun(@le, rr, x-1));
        S_x = sum(F_x, dimK);
        %         S_x = (cumsum(g_x));
    end

end
