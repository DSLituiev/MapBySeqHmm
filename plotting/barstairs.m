function varargout = barstairs(x, pdx, color, varargin)
% BARSTAIRS plots a bar graph based on MATLAB 'stairs' function
% allows usage of transparency
%
% ========= INPUT ===========
% 
% - x
% - pdf(x)
% - color
% 
% - optional 

[xb, yb] = stairs(x, pdx);
if  all( abs(diff(x)-diff(x(1:2))) < 1e-13)
    dx = diff(x(1:2));
    xb = xb- dx/2;
end
    
minY = min(-10, min(yb(:)));

if nargin>2
    pp = patch([xb(:); xb(end)+diff(xb(end-2:end-1));xb(end)+diff(xb(end-2:end-1)); xb(1)], [yb(:); yb(end);minY; minY],...
        color, 'edgecolor', color , varargin{:});
else
    
    pp = patch([xb(:); xb(end)+diff(xb(end-2:end-1));xb(end)+diff(xb(end-2:end-1)); xb(1)], [yb(:); yb(end);minY; minY],...
        color, 'edgecolor', color);
end
% ylim([min([yb(:);currYLims(1)]), max([yb(:);currYLims(2)])])

if nargout>0
    varargout = {pp, xb};
end