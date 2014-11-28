function [ f, medianP, upper_quantile ] = plotBatchProbT( t, P, fig_name, varargin )
%PLOTBATCHPROB Summary of this function goes here
%   Detailed explanation goes here

if nargin>3 && isscalar(varargin{1})
    qu = varargin{1};
else
    qu = 0.05;
end
marginal = calcMarginal( P , 2 ) - log10(size( P , 2 ));
medianP = median( P, 2 );
upper_quantile = quantile( P, 1-qu, 2 );
f = figure('name', fig_name );
plot(t, P, 'color', [1,1,1]*0.8 )
hold all
% p_arithm = plot(t, marginal, ':', 'color', [0,1,0]*0.8, 'linewidth', 2.5 );
p_geom   = plot(t, mean( P, 2 ), 'b--', 'linewidth', 2 );
p_median = plot(t, medianP, 'r-', 'linewidth', 2 );
p_quart  = plot(t, quantile( P, qu, 2 ), '-', 'color', [1,0.3,0.3]*0.8, 'linewidth', 2 );
plot(t, upper_quantile, '-', 'color', [1,0.3,0.3]*0.8, 'linewidth', 2 )
%% y-axis limits
ylims = [floor(min(quantile( P, 0.05, 2 ))), ceil(max( max(P(:)), 0 ))];
set(gca, 'ylim', ylims )
set(gca, 'ytick', 2*floor(ylims(1)/2):2:ylims(2) )
%% make log-scale axes
ax0 = gca;
ax1 = takeYAxisLogScale( ax0 );
set(ax1, 'box', 'on', 'tickdir', 'out');

ylabel(ax1, fig_name)
xlabel(ax0, 'genetic distance, cM')

leg = legend([p_geom, p_median, p_quart], {'geometric mean', 'median', sprintf('%u%% bounds', qu*100)});
set(leg, 'location', 'southwest', 'EdgeColor','white')

% axes(ax0);
end

