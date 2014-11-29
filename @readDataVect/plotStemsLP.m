function plotStemsLP(obj, varargin)

if nargin>2 && isscalar(varargin{1}) && isnumeric(varargin{1})
    chr0 = varargin{1};
    fieldPlotFun = @(x)obj.plotOneChromosome(chr0, x{:});
    args = varargin(2:end);
else
    chr0 = (1:obj.chrNumber)';
    fieldPlotFun = @(x)obj.plotChromosomes(x{:});
    args = varargin;
end
    
markerSz = 5;

%             plotChromosomes(obj, 'xPosteriorNorm', ...
%                 'plotfun', @(x,y)stem(x,y, 'MarkerSize', markerSz, 'MarkerEdgeColor', 'r', 'linewidth', 0.8),... % 'MarkerEdgeColor', 'r',
%                 'exp10',false, 'yscale', 'lin', 'figure', 'new', varargin{:});
%
%             set(obj.prevLine(:, end), 'BaseValue',( obj.prevYLims(1)-2) );

fieldPlotFun({'xPselNorm', ...
    'plotfun', @(x,y)stem(x,y, 'MarkerSize', markerSz-1, 'MarkerEdgeColor', 'b', 'linewidth', 0.8),... % 'MarkerEdgeColor', 'r',
    'exp10',false, 'yscale', 'lin', 'yThr', 0, 'figure', 'new', args{:}});

set(obj.prevLine(chr0, end), 'BaseValue',( obj.prevYLims(1)-2) );

hits = (obj.xPosteriorNorm > obj.xPselNorm );
%  colors = bsxfun(@times, [0 .8 0], hits) + bsxfun(@times, [.2 1 .2], 0.5 * ~hits) ;
% 
%  fieldPlotFun({ 'xPselNorm', ...
%      'plotfun', @(x,y,z)scatter(x,y, markerSz, colors(z)),... % 'MarkerEdgeColor', 'r',
%      'exp10',false, 'yscale', 'lin',  args{:}});

colors = [.8 0.3 0.3]; %colors = bsxfun(@times, [1 0.3 0.3], hits) + bsxfun(@times, [1 0 0], ~hits) ;
fieldPlotFun({'xPosteriorNorm', ...
    'plotfun', @(x,y,z)scatter(x,y, markerSz*5, colors, 'o', 'fill'),... % 'MarkerEdgeColor', 'r',
    'exp10',false, 'yscale', 'lin',  args{:}});


set(obj.prevAxes(chr0, end), 'TickDir', 'out')
end
