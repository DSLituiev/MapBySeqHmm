function f_out = plotScatterMW(obj, varargin)
% PLOTSCATTERMW plots scatter plots statistics of mutant vs wild type reads
% 
% ## INPUT
% chromosome (optional)  -- to plot data only for one chromosome
%

if ~obj.flagWT
    return
end
if nargin>1 && isscalar(varargin{1}) && isnumeric(varargin{1})
    chr0 = varargin{1};
    inds = obj.chromosome == chr0 ;
else
    inds = true(size(obj.x));
end

f_out = figure('name', 'mt vs wt SNP ratio');
da(1) = scatter(obj.f(inds), obj.fw(inds), 'b.');

hold on

da(2) =  scatter(obj.fmMedianF(inds ), obj.fwMedianF(inds ), 4, [0.7,.3, .3], 'o');
da(3) =  scatter(obj.fmMeanF(inds ), obj.fwMeanF(inds ), 4, [0.2,.8, .2], 'o');
daInfo = {'raw [all]', 'median', 'mean'};

dl(1) = plot([0,1], [0,1], 'm--');
dl(2) = plot([0,.5], [0.5, 0], 'r-');
dlInfo = {'m = w (errors)', 'm + w = 1/2 (expected)'};

deInfo = {};
if ~isempty(obj.snpEcotypesInfo)
    de(1) =  scatter(obj.f(inds & obj.snpEcotypesInfo), obj.fw(inds & obj.snpEcotypesInfo), 'g.');
    deInfo = 'BG ecotype';
    uistack(da(2:3), 'top')
    uistack(dl, 'top')
end

axis equal
xlabel('f_{mt}'); ylabel('f_{wt}')
xlim([0,1]); ylim([0,1])
if ~isempty(obj.snpEcotypesInfo)
    legend([da, dl, de], {daInfo{:}, dlInfo{:}, deInfo}, 'location', 'northwest')
else
    legend([da, dl], {daInfo{:}, dlInfo{:}}, 'location', 'northwest')
end
end