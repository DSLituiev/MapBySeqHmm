function LL = logBetaBinomialThetaMu(k, n, mu, theta)
% logBetaBinomialThetaMu calculates log-likelihood (10-based) of
% Beta-Binomial distribution with given parameters mu and theta
% mu = alpha/(alpha + beta)
% theta = 1/(alpha + beta)
%
% References: 
% Smith, D.M. "Algorithm AS 189: Maximum Likelihood Estimation  of the 
% Parameters  of the Beta Binomial Distribution". Journal of the Royal 
% Statistical Society. Series C (Applied Statistics),  Vol. 32, No. 2 
% (1983), pp. 196-204;  URL: http://www.jstor.org/stable/2347299


dimK = find(size(k)>1);
if isempty(dimK)||numel(dimK)>1
    dimK = find(size(theta)==1, 1, 'first');
end

dims = size(bsxfun(@plus, theta, k));
S.type = '()';
S.subs = cellstr(repmat(':', [numel(dims),1]));
if size(k)~=size(n)
    error('maxLlThetaBetaBinomial:k_and_n_mismatch', 'k and n must have the same size')
end

rr = permute(0:(max(n)-1), [1,3,2]); % in the third dimension

S_k   = cumS(k, rr, dimK);
S_n_k = cumS(n-k, rr,dimK);
S_n   = cumS(n, rr, dimK);

rrTheta = bsxfun(@times, rr, theta);

LL =  bsxfun(@plus ,...
    sum( ...
    bsxfun(@minus, ...
    bsxfun(@times, S_k  , log10( bsxfun(@plus ,     mu, rrTheta) ) ) + ...
    bsxfun(@times, S_n_k, log10( bsxfun(@plus , 1 - mu, rrTheta) ) ) , ...
    bsxfun(@times, S_n  , log10( 1                      + rrTheta) ) ),...
    3 ),...
    log10(exp(1))* sum(gammaln(n + 1)-gammaln(k + 1)-gammaln(n - k + 1) , dimK ) ) ;
%   sum( log10(binomial(n,k)), dimK ) ) ;


%% subfunction
    function S_x = cumS(x, rr, dimK)
        F_x = uint32(bsxfun(@le, rr, x-1));
        S_x = sum(F_x, dimK);
    end

end
