function [ f, medianP ] = plotBatchProbT( t, P, fig_name )
%PLOTBATCHPROB Summary of this function goes here
%   Detailed explanation goes here

marginal = calcMarginal( P , 2 ) - log10(size( P , 2 ));
medianP = median( P, 2 );

f = figure('name', fig_name );
plot(t, P, 'color', [1,1,1]*0.8 )
hold all
p_arithm = plot(t, marginal, ':', 'color', [0,1,0]*0.8, 'linewidth', 2.5 );
p_geom   = plot(t, mean( P, 2 ), 'b--', 'linewidth', 2 );
p_median = plot(t, medianP, 'r-', 'linewidth', 2 );
p_quart  = plot(t, quantile( P, 0.25, 2 ), '-', 'color', [1,0.3,0.3]*0.8, 'linewidth', 1.2 );
plot(t, quantile( P, 0.75, 2 ), '-', 'color', [1,0.3,0.3]*0.8, 'linewidth', 1.2 )
%% y-axis limits
set(gca, 'ylim', [min(quantile( P, 0.05, 2 )), max( max(P(:)), 0 )])
%% make log-scale axes
ax0 = gca;
ax1 = takeYAxisLogScale( ax0 );
set(ax1, 'box', 'on');

leg = legend([p_arithm, p_geom, p_median, p_quart], {'arithmetic mean', 'geometric mean', 'median', 'quartiles'});

set(leg, 'location', 'southwest')

axes(ax0);
end

