function varargout = myhistlog(x,y, varargin)


% p = inputParser;
% parse(p, x,y, varargin{:});

pdx = hist(y,  x );
[xb,yb] = stairs(x, pdx');
pp = patch([xb(:,1); 1], [yb(:,1); 0], varargin{:});
if nargout>0
    varargout{1} = pp;
end