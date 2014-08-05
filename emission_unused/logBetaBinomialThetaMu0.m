function [ LL] = logBetaBinomialThetaMu0(k, n, mu, theta)
% maxLlThetaBetaBinomial(k, n, mu, theta):
% calculates log-likelihood (10-based) of Beta-Binomial distribution 
% with given parameters mu and theta
% mu = alpha/(alpha + beta)
% theta = 1/(alpha + beta)
%

rr =  permute((0:1:(max(n(:))-1))', [1,3,2]);

n_k = n-k;

cumLL_n  =  mySeries(rr,             1,      theta(:)');
cumLL_k  =  mySeries(rr( 1:max(k(:))   ),  mu(:)',     theta(:)');
cumLL_n_k = mySeries(rr( 1:max(n_k(:)) ) , 1 - mu(:)',  theta(:)');

LL_n = cumLL_n(n);
LL_k   = zeros(numel(k), max(numel(theta), numel(mu)) );
LL_n_k = zeros(numel(k), max(numel(theta), numel(mu)));

LL_k(k>0,:)   = cumLL_k(k(k>0),:);
LL_n_k(n_k>0, :) = cumLL_n_k(n_k(n_k>0),:);

%%
LL  = bsxfun(@plus, ...
    bsxfun(@minus, LL_k   + LL_n_k , LL_n ), ...
    log10(exp(1))* (gammaln(n + 1)-gammaln(k + 1)-gammaln(n - k + 1) ) );
% dLdTheta = Res_k + Res_n_k - Res_n;
% L   = bsxfun(@times, binomial(n(:), k(:)), L_k.*L_n_k ./ L_n );

    function [s] = mySeries(r, m, th)
        %         R  = sum( bsxfun(@rdivide, x, bsxfun(@plus, m, bsxfun(@times, x, t))), 3);
        s = cumsum( log10( bsxfun(@plus, m, bsxfun(@times, r, th))  ) , 1);
        %         L  = prod( (  m + bsxfun(@times, x, t) ) , 3);
    end
end