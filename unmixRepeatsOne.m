function [ dataStruct, mu, iTrue2, mixtObj, varargout ] = unmixRepeatsOne( dataStruct, firstParam, varargin )
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
% dataStruct       -- a structure with the obligatory field 'r'
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
% 'C_CONTR' : 1/16 -- 'rigidness' of the margin between true and spurious reads
% 'B_CONTR' : 0    -- the lower 'contribution' value for spurious reads
% 'peakVar1': 'max'-- function to use for finding the 'true' peak (feature #1)
% 'peakVar2': 'min'-- function to use for finding the 'true' peak (feature #2)
%              possible values:
%             'max', 'min', 'median', 'none'
%
% ====================== OUTPUT  ======================
%
% dataStruct       -- the structure given in the input with a new field
%                    'contrib' for responsivity/contribution/ believe
%                     that each read comes from the 'true' distribution
% mu               -- a structure of means with two vector fields
%                     corresponding to the features used for classification
% iTrue2           -- the subscipt of the 'true' mean in the vector 'mu'
%

%% check the input parameters
p = inputParser;
addRequired(p, 'OR', @(x)(isstruct(x) || isobject(x)) );
addRequired(p, 'firstParam', @(x)(isfield(dataStruct, x) || isprop(dataStruct, x)) );
%
addOptional(p,'plotFlag', '', @(x)(ischar(x)|isempty(x)) );
%
addParamValue(p,     'modeNum',             4, @isscalar);
addParamValue(p,     'repeateEM',          20, @isscalar);
addParamValue(p,     'C_CONTR',          1/16, @isscalar);
addParamValue(p,     'B_CONTR',             0, @isscalar);
addParamValue(p,     'peakVar1',            'max',...
    @(x)(isa(eval(['@',x]), 'function_handle')) || strcmpi(x, 'none') );
addParamValue(p,     'peakVar2',            'max', ...
    @(x)(isa(eval(['@',x]), 'function_handle')) || strcmpi(x, 'none') );

parse(p, dataStruct, firstParam, varargin{:});

%% assign optional variables
repeateEM = p.Results.repeateEM ;
modeNum = p.Results.modeNum;
if modeNum<2
    modeNum = 2;
end
%%
%= select the first feature according to the input variable
x1 = double(dataStruct.(firstParam));
%= combine the log-scaled feature vectors into one array
    logX1 = log10( x1(~isnan(x1)) );
    logX2 = log10(double(dataStruct.r(~isnan(x1)) ) );

XX = [logX1];
%= select non-repeatitive reads
ingXX = XX( dataStruct.notaRepeat, : );

% st.mu = [2 6; 5 1];
% st.Sigma = cat(3,[.2 0;0 1],[.1 0;0 1]);
% st.PComponents = ones(1,2)/2;
% obj = gmdistribution.fit(XX ,k, 'Start', st);

%% calculate GMM for a number of mixture components ii = 2...modeNum
dbclear if warning
warning('off','stats:gmdistribution:FailedToConverge')
for ii = modeNum:-1:1
    obj0 = gmdistribution.fit(XX , ii, 'Replicates', repeateEM);
    mixt(ii).obj = obj0;
    mixtBICs(ii) = obj0.BIC;
end
%= find the BIC-optimal number of the mixture components ( >=2 )
[~, iMaxObj] = min( mixtBICs );
%= use the BIC-optimal GMM parameters
mixtObj = mixt( iMaxObj ).obj;

clear ii obj0 mixt iMaxObj

%=====================
if strcmpi(p.Results.peakVar2 , 'none')
    iTrue2 = NaN;
    iFalse2 = NaN;
else
    if strcmpi(p.Results.peakVar2 , 'max')
        [~, iTrue2] = max(mixtObj.mu(:,1)); % dx
        [~, iFalse2] = min(mixtObj.mu(:,1)); % dx
    elseif strcmpi(p.Results.peakVar2 , 'min')
        [~, iTrue2] = min(mixtObj.mu(:,1)); % r
        [~, iFalse2] = max(mixtObj.mu(:,1)); % r
    elseif strcmpi(p.Results.peakVar2 , 'median')
        iTrue2 = floor(numel(mixtObj.mu(:,1))/2)+1; % dx
        iFalse2 = NaN;
    end
end
%=====================
%= calculate the odds that the SNP locus comes from the 'true' distribution

    piTrue = mixtObj.PComponents(iTrue2);
    PrTrue = mvnpdf(XX, mixtObj.mu(iTrue2,:), mixtObj.Sigma(:,:,iTrue2));
    Odds = (1-piTrue)* PrTrue./(pdf(mixtObj, XX) - piTrue * PrTrue );
    Contr = (piTrue)* PrTrue./pdf(mixtObj, XX);
%  Pr2 = mvnpdf(X0, obj.mu(1,:), obj.Sigma(:,:,1));
%  Pr1 = mvnpdf(X0, obj.mu(2,:), obj.Sigma(:,:,2));
%  Odds0 = Pr2./Pr1;

% calculate 'contribution'/'responsivity' based on odds
% figure; hist(Contr)
dataStruct.contrib = Contr./max(Contr);% B_CONTR + (1 - B_CONTR)./(1 + C_CONTR./Odds);

if any(dataStruct.contrib > 1)
    warning('CumMatr:ENegInput', 'contrib > 1 !')
end

if isempty(dataStruct.notaRepeat)
   dataStruct.notaRepeat = true(numel(logX1),1) ;
end

%% assign the mean values
mu.(firstParam) =  mixtObj.mu(:,1) ;

%% visualize the results if requested
if strcmpi(p.Results.plotFlag, 'plot')
    
    %==
    f(1) = figure('name', 'mixture: estimated PDFs for non-repetitive regions');
    scatter( logX1(dataStruct.notaRepeat),  logX2(dataStruct.notaRepeat) , 3,...
        pdf(mixtObj, [logX1(dataStruct.notaRepeat)]) );
    hold all
    %     plot(model.m(1,:), model.m(1,:), 'rx', 'MarkerSize', 7)
    xlabel( firstParam )
    ylabel('\it r')
    title('non-repeatitive')
    ax(1) = gca;
    %=======================================
    f(2) = figure('name', 'contribution assigned to the components');
    scatter( logX1,  logX2, 3,  dataStruct.contrib);
    xLim = get(gca, 'xlim');
    yLim = get(gca, 'ylim');
    nPoints = 50;
    cLevels = [1, 5, 10, 25, 50, 75, 90, 95, 99];
    hold all
    plot( mixtObj.mu(:,1), median(logX2)*ones(size(mixtObj.mu(:,1))), 'kx', 'MarkerSize', 8, 'LineWidth', 2)
    [X, Y] = meshgrid(linspace(xLim(1), xLim(2), nPoints), linspace(yLim(1), yLim(2), nPoints) );
    
    h = contour(X, Y, reshape(pdf(mixtObj, [X(:)]), [nPoints,nPoints]), cLevels );
    %     h = ezcontour(@(x,y)pdf(obj, [x y]), [xLim, yLim]);
    if exist('h', 'var') && all(ishandle(h(:)))
        set(h, 'linewidth', 2);
        set(h, 'LevelList', [0.5, .75, .9, .95]);
    end
    xlabel( firstParam )
    ylabel('\it r')
    title('')
    ax(2) = gca;

    %=======================================
    f(3) = figure('name', 'mixture: estimated PDFs, scatter plots');
    %---------------------
    sp(1) = subplot(3,1,1);    
    scatter(  logX1,  logX2, 3);
%     scatter(  logX1,  logX2, 3,...
%         pdf(mixtObj, [logX1]) );
    hold all
    plot(mixtObj.mu(:,1), median(logX2), 'wx', 'MarkerSize',10, 'LineWidth', 4)
    plot(mixtObj.mu(:,1), median(logX2), 'kx', 'MarkerSize',8, 'LineWidth', 2)
 %  xlabel( firstParam )
    ylabel('\it r')
    title('all')
    %---------------------
    sp(2) = subplot(3,1,2);
    scatter( logX1(dataStruct.notaRepeat), logX2(dataStruct.notaRepeat) , 3)
%     scatter( logX1(dataStruct.notaRepeat), logX2(dataStruct.notaRepeat) , 3,...
%         pdf(mixtObj, [logX1(dataStruct.notaRepeat) ]) );
    hold all
    plot(mixtObj.mu(:,1), median(logX2), 'wx', 'MarkerSize',10, 'LineWidth', 4)
    plot(mixtObj.mu(:,1),  median(logX2), 'kx', 'MarkerSize', 8, 'LineWidth', 2)
 %   xlabel( firstParam )
    ylabel('\it r')
    title('non-repeatitive')
    %---------------------
    sp(3) = subplot(3,1,3);
    scatter(   logX1(~dataStruct.notaRepeat),  logX2(~dataStruct.notaRepeat) , 3);
%     scatter(   logX1(~dataStruct.notaRepeat),  logX2(~dataStruct.notaRepeat) , 3,...
%         pdf(mixtObj, [ logX1(~dataStruct.notaRepeat) ]) );
    hold all
    plot(mixtObj.mu(:,1), median(logX2), 'wx', 'MarkerSize',10, 'LineWidth', 4)
    plot(mixtObj.mu(:,1),  median(logX2), 'kx', 'MarkerSize', 8, 'LineWidth', 2)
    xlabel( firstParam )
    ylabel('\it r')
    title('repeatitive')
    %-------------------
    axis(sp, 'equal')
     set(sp, 'xlim', [ .5*floor(2*min(logX1)), .5*ceil(2*max(logX1)) ])
     set([sp(:); ax(:)], 'ylim', [ floor(min(logX2)), ceil(quantile(logX2, .99)) ])
     for ii = 1:3;  box(sp(ii), 'on') ; end
%         set(sp, 'ylim', [ floor(min(logX2)), ceil(max(logX2)) ])

%=======================================
nBins = 25;
f(4) = figure('name', 'mixture: empirical PDFs');
    [yh, xh] = hist(  logX1, nBins );
    yMax = max(yh);
    b(1) = barstairs( xh, yh./sum(yh), 'b');
    hold all
    %---------------------
    [yhR, xhR] = hist( logX1(dataStruct.notaRepeat), xh);
    b(2) = barstairs( xhR, yhR./sum(yh), 'r', 'faceColor', 'none', 'linewidth', 3);
    %---------------------
    [yhN, xhN] = hist( logX1(~dataStruct.notaRepeat),  xh);
    b(3) = barstairs( xhN, yhN./sum(yh), 'g', 'faceColor', 'none', 'linewidth', 3);
    %-------------------
    xlabel( firstParam )
%     ylabel('count')
    ylabel('frequency')
    legend(b, {'all', 'non-repeatitive', 'repeatitive' } ,'Location', 'NorthWest' )
%         set(sp, 'ylim', [ floor(min(logX2)), ceil(max(logX2)) ])
    plot(mixtObj.mu(1), 0, 'wx',  'MarkerSize', 10, 'LineWidth', 4)
    plot(mixtObj.mu(2), 0, 'wx',  'MarkerSize', 10,  'LineWidth', 4)
    hold all
    plot(mixtObj.mu(1), 0, 'kx', 'MarkerSize', 8, 'LineWidth', 2)
    plot(mixtObj.mu(2), 0, 'kx', 'MarkerSize', 8,  'LineWidth', 2)    
    plot(median(mixtObj.mu(1:2)), 0, 'k+', 'MarkerSize', 8,  'LineWidth', 2)
%     set(gca, 'ylim', [0, yMax])
   set(gca, 'ylim', [0, max(0.12, yMax./sum(yh) )])
   set(gca, 'xlim', [0, 6])

   
f(5) = figure('name', 'mixture: empirical CDFs');
    b(1) = barstairs( xh, cumsum(yh)./sum(yh), 'b');
    hold all
    plot(mixtObj.mu(1), 0, 'kx', 'MarkerSize', 8, 'LineWidth', 2)
    plot(mixtObj.mu(2), 0, 'kx', 'MarkerSize', 8,  'LineWidth', 2)
    plot(median(mixtObj.mu(1:2)), 0, 'k+', 'MarkerSize', 8,  'LineWidth', 2)
   set(gca, 'ylim', [0, 1])
   set(gca, 'xlim', [0, 6])
   xlabel( firstParam )
%     ylabel('count')
    ylabel('empirical CDF')    
    
if nargout>4
   varargout = {f};
end
end

end

