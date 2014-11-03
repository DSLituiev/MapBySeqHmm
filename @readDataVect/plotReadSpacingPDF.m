function [ f ] = plotReadSpacingPDF(obj, mixtObj,  varargin)
%UNTITLED  plots empirical and theoretical PDFs for read spacing

%% check the input parameters
p = inputParser;
addRequired(p, 'obj', @(x)(isstruct(x) || isobject(x)) );
addRequired(p, 'mixtObj', @(x)(isstruct(x) || isobject(x)) );
%
addOptional(p, 'dxMode', '', @(x)(ischar(x)|isempty(x)) );
addOptional(p, 'plot_theor_1', '', @(x)(ischar(x)|isempty(x)) );
%
addParamValue(p,     'READ_LENGTH',          130, @isscalar);
addParamValue(p,     'NUMBER_OF_SNPS_PER_READ',          2, @isscalar);
addParamValue(p,     'N_BINS',          25, @isscalar);
parse(p, obj, mixtObj, varargin{:});
%% process the input flags
if isempty(p.Results.plot_theor_1)
    flagPlotTheor1 = false;
elseif ~isempty(regexpi(p.Results.plot_theor_1), 'plot')
    flagPlotTheor1 = true;
end

if isempty(p.Results.dxMode) || ~isempty(regexpi(p.Results.dxMode, 'min'))
    flagDxMin = true;
    fprintf('plotting min(dx)\n')
elseif ~isempty(regexpi(p.Results.dxMode, 'all'))
    flagDxMin = false;
    fprintf('plotting all dx\n')
end
%% parameter for theoretical distribution
muInvPoisson = sum(obj.cMaxX)./ obj.Mtot ;

%% indices of true and spurious modes:
[~, iTrue] = max(mixtObj.mu(:,1)); % dx
[~, iFalse] = min(mixtObj.mu(:,1)); % dx

%% decide if the dx for repeat sequences is to be plotted
flagRepeats = sum(~obj.notaRepeat)./numel(obj.notaRepeat) > 0.01;

%% assign the variable to plot [dx or min(dx)]
if flagDxMin
    log_x = log10( obj.dx);
else
    dxRaw = diff(double(obj.x));
    dxRaw(dxRaw < 0) = NaN;
    log_x = log10( dxRaw);
end

log_x = log_x(~isnan(log_x));

XLIM = [0, 0.1*ceil(max(10*log_x))];
%===========================================================
%% plotting
f= figure('name', 'mixture: empirical and theoretical PDFs');
[yh, xh] = hist(  log_x, p.Results.N_BINS );
yMax = max(yh);
b(1) = barstairs( xh, yh./sum(yh), [0.2, 0.4, .8]);
hold all
legendText = {'empirical PDF'};
%---------------------
if  flagRepeats
    [yhR, xhR] = hist( log_x(obj.notaRepeat), xh);
    b(2) = barstairs( xhR, yhR./sum(yh), [0.1, 0.9, .3], 'faceColor', 'none', 'linewidth', 2);
    %---------------------
    [yhN, xhN] = hist( log_x(~obj.notaRepeat),  xh);
    b(3) = barstairs( xhN, yhN./sum(yh), [1, 0.15, 0.1], 'faceColor', 'none', 'linewidth', 2);
    legendText = {'empirical PDF (all)',  'non-repeatitive', 'repeatitive'};
end
%-------------------
dx = diff(xh(2:3));
xh = 0:dx/5:max(xh);
if flagPlotTheor1
    logP(:,2) = probSnpLogSpacingInMissMappedReads( 10.^xh', p.Results.READ_LENGTH, p.Results.NUMBER_OF_SNPS_PER_READ, 'uniform') .* ...
        mixtObj.PComponents(iFalse);
end
logP(:,1) = logPoisson1(10.^xh', muInvPoisson , 'i') .* mixtObj.PComponents(iTrue) ;
logP = logP .* dx;
logP = sum(logP, 2);
b(4) = plot(xh', logP, ':', 'color' , [.8, 0.3, .8], 'linewidth', 2);
b(5) = plot(xh', pdf(mixtObj, xh').* dx, '--', 'color' , 0.15*[1 1 1], 'linewidth', 2);


legendText = [legendText, {'theoretical PDF', 'estimated log-normal PDF'}];
fprintf('parameters of the theoretical distributions:\n')
muClustered =  p.Results.READ_LENGTH/ p.Results.NUMBER_OF_SNPS_PER_READ;
fprintf('L_R = %u\t M_R = %u\t mu_R = %u\n', p.Results.READ_LENGTH, p.Results.NUMBER_OF_SNPS_PER_READ, muClustered)
fprintf('L_G = %u\t M_G = %u\t mu_G = %u\n',  sum(obj.cMaxX), obj.Mtot, round(muInvPoisson))
if  flagRepeats
    legend(b, legendText ,'Location', 'NorthEast' )
else
    legend(b([1,4:end]), legendText ,'Location', 'NorthEast' )
end
%     b(5) = plot(xh, logP, ':', 'color' , [.8, 0.3, .8], 'linewidth', 2);

%-------------------
%     xl = xlabel( firstParam );
xl = xlabel( '${\Delta}x$', 'interpreter', 'latex');
set(xl, 'Units', 'Normalized');
pos = get(xl, 'Position');
set(xl, 'Position', pos + [0, -0.03, 0]);

%     ylabel('count')
ylabel('frequency')
%         set(sp, 'ylim', [ floor(min(logX2)), ceil(max(logX2)) ])
doubleMarker0(mixtObj.mu(1), 'w')
doubleMarker0(mixtObj.mu(2), 'w')
doubleMarker0(log10(muInvPoisson), 'm', '+')
if flagPlotTheor1
    doubleMarker0(log10(muClustered), 'm', '+')
end

doubleMarker0(median(mixtObj.mu(1:2)),'w','d', 5)


yMax = 1/50*ceil(50*max([0.12, yMax./sum(yh), logP(:)'] ));
axo = gca;
set(axo, 'ylim', [0,  yMax])
set(axo, 'xlim', XLIM)

%     xTicks = (get(gca, 'xticklabel'));
%     xTicksNew = cell(size(xTicks));
%     for ii=1:numel(xTicks)
%         xTicksNew{ii} = ['10^{', num2str(xTicks(ii)), '}'];
%     end
%
set(axo, 'xticklabel','')
axn = axes('xlim', 10.^get(gca,'xlim'),'tickdir', 'out', 'xscale', 'log', 'ytick', [],  'Color', 'None');
uistack(axn, 'down')

    function doubleMarker0(x, color, varargin)
        if nargin < 3
            markerType = 'x';
        else
            markerType = varargin{1};
        end
        if nargin < 4
            MarkerSize = 10;
        else
            MarkerSize  = varargin{2};
        end
        LineWidth = 0.4 * MarkerSize;
        
        plot(x, 0, markerType, 'Color', color, 'MarkerSize', MarkerSize, 'LineWidth', LineWidth)
        plot(x, 0, markerType, 'Color', 'k',   'MarkerSize', 0.8*MarkerSize, 'LineWidth', 0.5*LineWidth)
    end
end

