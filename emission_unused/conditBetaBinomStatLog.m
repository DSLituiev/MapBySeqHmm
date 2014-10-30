function [pBB, pU, f] = conditBetaBinomStatLog(k, n, theta, N, varargin)

if nargin>5
    f  = varargin{1};
    Pf = varargin{2};
else
    f = 0:1/(2*N):1/2;
    %= P(f| stationary model)
    Pf = log10(StationaryDistr(N)');
end
%= P(k,n | f, theta)
LL = logBetaBinomialThetaMu0(k, n, f, theta);

%= P(k,n|theta, beta-binomial, stationary model)
pBB = calcMarginal( bsxfun(@plus, Pf, LL), 2);

%= P(k,n|theta, uniform, stationary model)
pU = -log10(N)+ calcMarginal(Pf, 2);

end