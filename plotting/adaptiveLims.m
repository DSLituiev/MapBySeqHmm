function [yLims, plotcutoff]  = adaptiveLims(logY, varargin)

if nargin>1 && isscalar(varargin{1})
    qu = varargin{1};
else
    qu = .95;
end
yLims = [floor(quantile(logY(~isinf(logY)), qu)), ceil(max(logY)) ];
plotcutoff = yLims(1);
