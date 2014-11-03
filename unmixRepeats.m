function [ obj, mu, iTrue, mixtObj, varargout ] = unmixRepeats( obj, firstParam, varargin )
% UNMIXREPEATS represents a bi-variate distribution of the given
% log-scaled feature and read log-number of reads per each locus as a
% mixture of Gaussians (GMM) and performs logistic classification
% based on the submitted assumption about the localization of the 'true' peak.
%
%     The Bayesian information criterium (BIC) is used to chose the number of
% mixing modes with the lowest allowed number of modes of 2.
%
% The logarithmic base of 10 is used throughout this script
%
% ====================== INPUT  ======================
%
% obj       -- a structure with the obligatory field 'r'
%                    (read number in the locus), which is used as
%                     the second feature for classification
%                     and  the field 'notaRepeat' (logical) used for
%                     visualization purposes
% firstParam       -- the name of the field to be used as
%                     the first feature for classification
%
% ____________ optional variables and flags ____________
% 'plot'           -- visualize the results
% 'modeNum': 4     -- maximum number of modes to be tested
% 'repeateEM':  20 -- number of iterations of GMM algorithm
% 'peakVar1': 'max'-- function to use for finding the 'true' peak (feature #1)
% 'peakVar2': 'min'-- function to use for finding the 'true' peak (feature #2)
%              possible values:
%             'max', 'min', 'median', 'none'
%
% ====================== OUTPUT  ======================
%
% obj       -- the structure given in the input with a new field
%                    'contrib' for responsivity/contribution/ believe
%                     that each read comes from the 'true' distribution
% mu               -- a structure of means with two vector fields
%                     corresponding to the features used for classification
% iTrue2           -- the subscipt of the 'true' mean in the vector 'mu'
%

%% check the input parameters
p = inputParser;
addRequired(p, 'OR', @(x)(isstruct(x) || isobject(x)) );
addRequired(p, 'firstParam', @(x)(isfield(obj, x) || isprop(obj, x)) );
%
addOptional(p,'plotFlag', '', @(x)(ischar(x)|isempty(x)) );
%
addParamValue(p,     'modeNum',             4, @isscalar);
addParamValue(p,     'repeateEM',          20, @isscalar);
addParamValue(p,     'peakVar1',            'max',...
    @(x)(isa(eval(['@',x]), 'function_handle')) || strcmpi(x, 'none') );
addParamValue(p,     'peakVar2',            'max', ...
    @(x)(isa(eval(['@',x]), 'function_handle')) || strcmpi(x, 'none') );

addParamValue(p,     'READ_LENGTH',          130, @isscalar);

parse(p, obj, firstParam, varargin{:});

%% assign optional variables
repeateEM = p.Results.repeateEM ;
modeNum = p.Results.modeNum;
if modeNum<2
    modeNum = 2;
end
%% log-transform the mixture variable:
%= select the first feature according to the input variable
x = double(obj.(firstParam));
%= combine the log-scaled feature vectors into one array
log_x = log10( x(~isnan(x)) );

% st.mu = [2 6; 5 1];
% st.Sigma = cat(3,[.2 0;0 1],[.1 0;0 1]);
% st.PComponents = ones(1,2)/2;
% distribObj = gmdistribution.fit(XX ,k, 'Start', st);

%% calculate GMM for a number of mixture components ii = 2...modeNum
dbclear if warning
warning('off','stats:gmdistribution:FailedToConverge')
for ii = modeNum:-1:1
    distribObj0 = gmdistribution.fit(log_x , ii, 'Replicates', repeateEM);
    mixt(ii).distribObj = distribObj0;
    mixtBICs(ii) = distribObj0.BIC;
end
%= find the BIC-optimal number of the mixture components ( >=2 )
[~, iMaxObj] = min( mixtBICs );
%= use the BIC-optimal GMM parameters
mixtObj = mixt( iMaxObj ).distribObj;

clear ii distribObj0 mixt iMaxObj

%=====================
if strcmpi(p.Results.peakVar2 , 'none')
    iTrue = NaN;
    iFalse = NaN;
else
    if strcmpi(p.Results.peakVar2 , 'max')
        [~, iTrue] = max(mixtObj.mu(:,1)); % dx
        [~, iFalse] = min(mixtObj.mu(:,1)); % dx
    elseif strcmpi(p.Results.peakVar2 , 'min')
        [~, iTrue] = min(mixtObj.mu(:,1)); % r
        [~, iFalse] = max(mixtObj.mu(:,1)); % r
    elseif strcmpi(p.Results.peakVar2 , 'median')
        iTrue = floor(numel(mixtObj.mu(:,1))/2)+1; % dx
        iFalse = NaN;
    end
end

%% calculate membership for each data point
%=   P[x| component_true ]
piTrue = mixtObj.PComponents(iTrue);
%=   P[x| component_true ]
PrTrue = mvnpdf(log_x, mixtObj.mu(iTrue,:), mixtObj.Sigma(:,:,iTrue));
%=  Membership = P[x| component_true ] * P[component_true] / P[x| mixture ]
Membership = (piTrue)* PrTrue./pdf(mixtObj, log_x);

%% calculate the odds that the SNP locus comes from the 'true' distribution
% Odds = (1-piTrue)* PrTrue./(pdf(mixtObj, XX) - piTrue * PrTrue );
%  Pr2 = mvnpdf(X0, distribObj.mu(1,:), distribObj.Sigma(:,:,1));
%  Pr1 = mvnpdf(X0, distribObj.mu(2,:), distribObj.Sigma(:,:,2));
%  Odds0 = Pr2./Pr1;

% figure; hist(Membership)
obj.contrib = Membership./max(Membership);

if any(obj.contrib > 1)
    warning('CumMatr:ENegInput', 'contrib > 1 !')
end

if isempty(obj.notaRepeat)
    obj.notaRepeat = true(numel(log_x),1) ;
end

%% assign the mean values to the output object
mu.(firstParam) =  mixtObj.mu(:,1) ;
%% parameter for theoretical distribution
muInvPoisson = sum(obj.cMaxX)./ obj.Mtot ;

%% visualize the results if requested
if strcmpi(p.Results.plotFlag, 'plot')
    
    log_r = log10(double(obj.r(~isnan(x)) ) );
    %=======================================
    f(1) = figure('name', 'mixture: estimated PDFs for non-repetitive regions');
    scatter( log_x(obj.notaRepeat),  log_r(obj.notaRepeat) , 3,...
        pdf(mixtObj, log_x(obj.notaRepeat)) );
    hold all
    %     plot(model.m(1,:), model.m(1,:), 'rx', 'MarkerSize', 7)
    xlabel( firstParam )
    ylabel('\it r')
    % title('non-repeatitive')
    ax(1) = gca;
    %=======================================
    f(2) = figure('name', 'contribution assigned to the components');
    scatter( log_x,  log_r, 3,  obj.contrib);
    xLim = get(gca, 'xlim');
    yLim = get(gca, 'ylim');
    nPoints = 50;
    cLevels = [1, 5, 10, 25, 50, 75, 90, 95, 99];
    hold all
    plot( mixtObj.mu(:,1), median(log_r)*ones(size(mixtObj.mu(:,1))), 'kx', 'MarkerSize', 8, 'LineWidth', 2)
    [X, Y] = meshgrid(linspace(xLim(1), xLim(2), nPoints), linspace(yLim(1), yLim(2), nPoints) );
    
    h = contour(X, Y, reshape(pdf(mixtObj, [X(:)]), [nPoints,nPoints]), cLevels );
    %     h = ezcontour(@(x,y)pdf(distribObj, [x y]), [xLim, yLim]);
    if exist('h', 'var') && all(ishandle(h(:)))
        set(h, 'linewidth', 2);
        set(h, 'LevelList', [0.5, .75, .9, .95]);
    end
    xlabel( firstParam )
    ylabel('\it r')
    title('')
    ax(2) = gca;
    
    %% plot scatter plots for repeat and unique sequences separately    
    flagRepeats = sum(~obj.notaRepeat)./numel(obj.notaRepeat) > 0.01;
    if  flagRepeats
        f(3) = figure('name', 'mixture: estimated PDFs, scatter plots');
        %---------------------
        sp(1) = subplot(3,1,1);
        scatter(  log_x,  log_r, 3);
        %     scatter(  logX1,  logX2, 3,...
        %         pdf(mixtObj, [logX1]) );
        hold all
        plot(mixtObj.mu(:,1), median(log_r), 'wx', 'MarkerSize',10, 'LineWidth', 4)
        plot(mixtObj.mu(:,1), median(log_r), 'kx', 'MarkerSize',8, 'LineWidth', 2)
        %  xlabel( firstParam )
        ylabel('\it r')
        title('all')
        %---------------------
        sp(2) = subplot(3,1,2);
        scatter( log_x(obj.notaRepeat), log_r(obj.notaRepeat) , 3)
        %     scatter( logX1(obj.notaRepeat), logX2(obj.notaRepeat) , 3,...
        %         pdf(mixtObj, [logX1(obj.notaRepeat) ]) );
        hold all
        plot(mixtObj.mu(:,1), median(log_r), 'wx', 'MarkerSize',10, 'LineWidth', 4)
        plot(mixtObj.mu(:,1),  median(log_r), 'kx', 'MarkerSize', 8, 'LineWidth', 2)
        %   xlabel( firstParam )
        ylabel('\it r')
        title('non-repeatitive')
        %---------------------
        sp(3) = subplot(3,1,3);
        scatter(   log_x(~obj.notaRepeat),  log_r(~obj.notaRepeat) , 3);
        %     scatter(   logX1(~obj.notaRepeat),  logX2(~obj.notaRepeat) , 3,...
        %         pdf(mixtObj, [ logX1(~obj.notaRepeat) ]) );
        hold all
        plot(mixtObj.mu(:,1), median(log_r), 'wx', 'MarkerSize',10, 'LineWidth', 4)
        plot(mixtObj.mu(:,1),  median(log_r), 'kx', 'MarkerSize', 8, 'LineWidth', 2)
        xlabel( firstParam )
        ylabel('\it r')
        title('repeatitive')
        %-------------------
        axis(sp, 'equal')
        set(sp, 'xlim', [ .5*floor(2*min(log_x)), .5*ceil(2*max(log_x)) ])
        set([sp(:); ax(:)], 'ylim', [ floor(min(log_r)), ceil(quantile(log_r, .99)) ])
        for ii = 1:3;  box(sp(ii), 'on') ; end
        %         set(sp, 'ylim', [ floor(min(logX2)), ceil(max(logX2)) ])
    end
    %% empirical and theoretical PDFs
    
    N_BINS = 25;
    f(4) = figure('name', 'mixture: empirical and theoretical PDFs');
    [yh, xh] = hist(  log_x, N_BINS );
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
    NUMBER_OF_SNPS_PER_READ = 2;
    
    dx = diff(xh(2:3));
    xh = 0:dx/5:max(xh);
    logP(:,1) = probSnpLogSpacingInMissMappedReads( 10.^xh', p.Results.READ_LENGTH, NUMBER_OF_SNPS_PER_READ, 'uniform') .* ...
        mixtObj.PComponents(iFalse);
    logP(:,2) = logPoisson1(10.^xh', muInvPoisson , 'i') .* mixtObj.PComponents(iTrue) ;
    logP = logP .* dx;
    logP = sum(logP, 2);
    b(4) = plot(xh', logP, ':', 'color' , [.8, 0.3, .8], 'linewidth', 2);
    b(5) = plot(xh', pdf(mixtObj, xh').* dx, '--', 'color' , 0.15*[1 1 1], 'linewidth', 2);
    
    
    legendText = {legendText{:}, 'theoretical PDF', 'estimated log-normal PDF'};
    fprintf('parameters of the theoretical distributions:\n')
    fprintf('L_R = %u\t M_R = %u\t mu_R = %u\n', p.Results.READ_LENGTH, NUMBER_OF_SNPS_PER_READ, p.Results.READ_LENGTH/ NUMBER_OF_SNPS_PER_READ)
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
    doubleMarker0(log10(p.Results.READ_LENGTH/NUMBER_OF_SNPS_PER_READ), 'm', '+')
    
    doubleMarker0(median(mixtObj.mu(1:2)),'w','d', 5)
    
    
    yMax = 1/50*ceil(50*max([0.12, yMax./sum(yh), logP(:)'] ));
    axo = gca;
    set(axo, 'ylim', [0,  yMax])
    set(axo, 'xlim', [0, 6])
    
    %     xTicks = (get(gca, 'xticklabel'));
    %     xTicksNew = cell(size(xTicks));
    %     for ii=1:numel(xTicks)
    %         xTicksNew{ii} = ['10^{', num2str(xTicks(ii)), '}'];
    %     end
    %
    set(axo, 'xticklabel','')
    axn = axes('xlim', 10.^get(gca,'xlim'),'tickdir', 'out', 'xscale', 'log', 'ytick', [],  'Color', 'None');
    uistack(axn, 'down')
    
    %% empirical CDFs
    %     f(5) = figure('name', 'mixture: empirical CDFs');
    %     b(1) = barstairs( xh, cumsum(yh)./sum(yh), 'b');
    %     hold all
    %     doubleMarker0(mixtObj.mu(1), 'w')
    %     doubleMarker0(mixtObj.mu(2), 'w')
    %     doubleMarker0( median(mixtObj.mu(1:2)), 'w', '+')
    %     set(gca, 'ylim', [0, 1])
    %     set(gca, 'xlim', [0, 6])
    %     xlabel( firstParam )
    %     %     ylabel('count')
    %     ylabel('empirical CDF')
    
    if nargout>4
        varargout = {f};
    end
end

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