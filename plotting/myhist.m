function varargout = myhist(x,y, varargin)

yLims = get(gca, 'ylim');

[pdx, x] = hist(y,  x );

if nargout>0
    varargout{1} = barstairs(x, pdx, varargin{:});
else
    barstairs(x, pdx, varargin{:});
end

set(gca, 'ylim', [0, max(yLims(2), ceil(.1*max(pdx))*10)] )