function varargout = myhist100(x,y, varargin)


if isempty(get(gca, 'Children'))
    yLims = [NaN, NaN];    
else
    yLims = get(gca, 'ylim');
end

[pdx, x] = hist(y,  x );
pdx = pdx'./sum(pdx);

if nargout>0
    varargout{1} = barstairs(x, pdx, varargin{:});
else
    barstairs(x, pdx, varargin{:});
end

set(gca, 'ylim',  [0, nanmax(yLims(2), ceil(10*max(pdx))/10)] )