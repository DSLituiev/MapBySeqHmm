function [ P ] = logPoisson1( x, mu , varargin)
%LOGPOISSON1(x, mu, ['inverse']) -- computes log-Poisson distribution for k = 1
% which is the logarithm of the interarrival time of the process.
%
%   P[y = log(x)] = log(10) * mu * x * exp( - mu * x )
%
%   Instead of mu, an inverse quantity (1/mu) can be supplied, with a flag
%   'inverse' or just 'i'
%

if nargin>2 && (varargin{1}(1) == 'i')
    mu_x = bsxfun(@rdivide, x, mu);
else
    mu_x = bsxfun(@times, x, mu);    
end    

P = log(10) * mu_x .* exp( - mu_x);

end

