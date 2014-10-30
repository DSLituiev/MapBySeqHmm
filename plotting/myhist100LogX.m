function varargout = myhist100LogX(x,y, varargin)


if isempty(get(gca, 'Children'))
    yLims = [NaN, NaN];
else
    yLims = get(gca, 'ylim');
end

[pdx, x] = hist(y,  x );
pdx = pdx'./sum(pdx);

if nargout>0
    varargout{1} = barstairs(10.^x, pdx, varargin{:});
else
    barstairs(10.^x, pdx, varargin{:});
end

if nargout>1
    [~, modeI] = max(pdx);
    if modeI< numel(x)
        varargout{2} = 10.^( 0.5*(x(modeI)+x(modeI+1))  );
    else
        varargout{2} = 10.^( (x(modeI))  );
    end
end

set(gca, 'ylim',  [0, nanmax(yLims(2), ceil(10*max(pdx))/10)] , 'xscale', 'log')