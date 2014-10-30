function [pBB, pU, f] = conditBetaBinomStat(k, n, theta, N, varargin)
% 
%   P[ k, n | theta,  model* ] = P[ k, n | f, theta] * P[ f | model* ]
%
%   * model: no selection
% 

if nargin>5
    f  = varargin{1};
    Pf = varargin{2};
else
    f = 0:1/(2*N):1;
    %= P(f| stationary model)
    Pf = StationaryDistr(N)';
end
%= P(k,n | f, theta)
LL = logBetaBinomialThetaMu0(k, n, f, theta);

%= P(k,n|theta, beta-binomial, stationary model)
pBB = sum(bsxfun(@times, Pf, 10.^LL(:, f<1/2)), 2);

%= P(k,n|theta, uniform, stationary model)

U = 1/(2*N);
pU = sum(10.^(LL + log(U)) , 2);
% sum(Pf) .* sum(10.^LL, 2)./(2*N);

end