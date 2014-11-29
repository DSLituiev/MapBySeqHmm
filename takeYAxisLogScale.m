function [ ax ] = takeYAxisLogScale( varargin )
%TAKEYAXISLOGSCALE Summary of this function goes here
%   Detailed explanation goes here

if nargin>0 && ishandle(varargin{1})
   ax0 = varargin{1};
else
    ax0 = gca;
end
    new_y_ticks = 10.^ get(ax0,'ytick');

if ~any( diff(new_y_ticks) <= 0)
    ax = axes('position', get(ax0, 'position'), 'ylim',  10.^get(ax0, 'ylim'),...
        'Color', 'none', 'xtick', [], 'yscale', 'log');
    set(ax, 'ytick', new_y_ticks )
    set(ax0, 'ytick', [])
else
    ax = ax0;
end

end

