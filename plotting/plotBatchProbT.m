function [ f, medianP, upper_quantile ] = plotBatchProbT( t, P, fig_name, varargin )
%PLOTBATCHPROB Summary of this function goes here
%   Detailed explanation goes here

if nargin>3 && isscalar(varargin{1})
    qu = varargin{1};
else
    qu = 0.05;
end

if  nargin>3 && any(strcmpi('nolog', varargin))
    flagNoLog = true;
else
    flagNoLog = false;
end

if  nargin>3 && any(strcmpi('cM', varargin))
    t = t*100;
    flag_cM = true;
    flag_Mb = false;
else
    flag_cM = false;
end

if  nargin>3 && any(strcmpi('Mb', varargin))
    t = t*1e-6;
    flag_cM = false;
    flag_Mb = true;
else
    flag_Mb = false;
end

% marginal = calcMarginal( P , 2 ) - log10(size( P , 2 ));
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
% ylims = get(gca, 'ylim'); %
ylims = [floor(min(quantile( P, 0.05, 2 ))), ceil(max( max(P(:)), 0 ))];
if diff(ylims) < 20
    set(gca, 'ylim', ylims )
    set(gca, 'ytick', 2*floor(ylims(1)/2):2:ylims(2) )
end
%% make log-scale axes
ax0 = gca;
if ~flagNoLog
    ax1 = takeYAxisLogScale( ax0 );
    set(ax1, 'box', 'on', 'tickdir', 'out');
    ylabel(ax1, fig_name)
else
    ylabel(ax0, fig_name)
end

if flag_cM
    xlabel(ax0, 'genetic distance, cM')
end
if flag_Mb
    xlabel(ax0, 'chromosomal position, Mb')
end


leg = legend([p_geom, p_median, p_quart], {'geometric mean', 'median', sprintf('%u%% bounds', qu*100)});
set(leg, 'location', 'southwest', 'EdgeColor','white')

axes(ax0);
end

