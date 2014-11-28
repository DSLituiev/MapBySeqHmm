function [ ax ] = takeYAxisLogScale( varargin )
%TAKEYAXISLOGSCALE Summary of this function goes here
%   Detailed explanation goes here

if nargin>0 && ishandle(varargin{1})
   ax0 = varargin{1};
else
    ax0 = gca;
end
    
ax = axes('position', get(ax0, 'position'), 'ylim',  10.^get(ax0, 'ylim'),...
    'Color', 'none', 'xtick', [], 'yscale', 'log');
set(ax, 'ytick', 10.^ get(ax0,'ytick') )
set(ax0, 'ytick', [])

end

