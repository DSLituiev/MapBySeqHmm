function globalStatSeq(ChrReads, varargin)
fontSize = 12;
% [f, spl]  = plotAllChrNt(ChrReads, 'f','xname','r','ylim',[0,1],...
%     'plotfun',@(x,y)scatter(x,y,5,'rx'), 'xscale', 'log', 'yscale', 'lin');
% plotAllChrNt(ChrReads, 'f','xname','r','ylim',[0,1],...
%     'plotfun',@(x,y)scatter(x,y,3,'bo'), 'xscale', 'log', 'yscale', 'lin',...
%     'xlim', [0, filt.r_max],...
%     'select', 'EMSqr','OldSP', spl, 'OldFig', f);
% legend({'raw','filtered'})
% 
% 
% [f, spl]  = plotAllChrNt(ChrReads, 'f','ylim',[0,1],...
%     'plotfun',@(x,y)scatter(x,y,5,'rx'), 'yscale', 'lin');
% plotAllChrNt(ChrReads, 'f','ylim',[0,1],...
%     'plotfun',@(x,y)scatter(x,y,3,'bo'), 'yscale', 'lin',...
%     'select', 'EMSqr','OldSP', spl, 'OldFig', f);
% 
% [f, spl]  = plotAllChrNt(ChrReads, 'Coverage','ylim',10.^[3,8],'yscale','log');
% [f, ~]    = plotAllChrNt(ChrReads, 'meanCov',...
%     'OldSP', spl, 'OldFig', f,'linestyle', ':','ylim',10.^[3,8],'plotfun',@(x,y)plot(x,y,'k:'));
% fig(f,'units','centimeters','height',20, 'width',20)
% set(f, 'color', 'none');


figure('name','distribution of the read frequencies')
df = 0.02;
fx = df/2:df:(1-df/2);
% hf = hist( fv, fx);
% hfEMS = hist( fvEMS, fx);
% stairs(fvEMS, fx, 'linewidth', 2)
 myhist100(fx,  ChrReads.f, 'r', 'edgecolor', 'none');
 hold on
 alpha(bb, 0.4)
 set(gca, 'xtick', 0:0.125:1, 'Layer', 'top', 'XGrid', 'on');
 xlim([0,1])
title('distribution of the read frequencies')
ylabel('P, %')
xlabel('frequency of mutant reads, %')
legend({'raw','filtered'})
%% 
ll = 2.^(0:.05:6)';
fx = flipud(1./ll);
hf = hist( ChrReads.f, fx)';
hfEMS = hist( fvEMS, fx)';
xticks = (2).^(0:1:6)';
% xticks = sort([xticks; 3/4],'descend')
% ll = (1:1:32)';
% fx = flipud( 1./ll);
% hf = hist( fv, fx)';
% hfEMS = hist( fvEMS, fx)';

figure('name','distribution of the read frequencies inv')
barstairs(flipud(1./fx), (flipud(hf))./sum(hf), 'r')
hold on
 bb = barstairs(flipud(1./fx), (flipud(hfEMS))./sum(hfEMS), 'g');
 alpha(bb, 0.4)
 set(gca, 'xlim', [min(ll), max(ll)], 'xscale', 'log')
set(gca, 'xtick', (xticks) )
set(gca,'XTickLabel',[])
xpos = log2(xticks)./diff(log2(get(gca, 'xlim')));
xtlab = strcat ('$', cellfun(@(x)sprintf('\\frac{1}{%u}',x), num2cell(xticks),'un',0) , '$');
text( xpos , zeros(size(xpos))-0.05,...
    xtlab , ...   %# create text at same locations
    'Interpreter','latex', ...                   %# specify tex interpreter
    'VerticalAlignment','Top', ...          %# v-align to be underneath
    'HorizontalAlignment','center',...         %# h-aligh to be centered
    'units', 'normalized', 'Rotation',0,...
    'FontSize',fontSize);
set(gca, 'Layer', 'top', 'XGrid', 'on');
title('distribution of the read frequencies')
ylabel('P, %')
legend({'raw','filtered'})
%% cumulative
ll = 2.^(0:.05:6)';
fx = flipud(1./ll);
hf = hist( fv, fx)';
hfEMS = hist( fvEMS, fx)';
xticks = (2).^(0:1:6)';

figure('name','cumulative distribution of the read frequencies inv')
barstairs(flipud(1./fx), 100*cumsum(flipud(hf))./sum(hf), 'r')
hold on
 bb = barstairs(flipud(1./fx), 100*cumsum(flipud(hfEMS))./sum(hfEMS), 'g');
 alpha(bb, 0.4)
 set(gca, 'xlim', [min(ll), max(ll)], 'xscale', 'log')
set(gca, 'xtick', (xticks) )
set(gca,'XTickLabel',[])
% xpos = get(gca, 'xtick')./diff(get(gca, 'xlim')); % linear x-axis
xpos = log2(xticks)./diff(log2(get(gca, 'xlim')));
xtlab = strcat ('$', cellfun(@(x)sprintf('\\frac{1}{%u}',x), num2cell(xticks),'un',0) , '$');
text( xpos , zeros(size(xpos))-0.05,...
    xtlab , ...   %# create text at same locations
    'Interpreter','latex', ...                   %# specify tex interpreter
    'VerticalAlignment','Top', ...          %# v-align to be underneath
    'HorizontalAlignment','center',...         %# h-aligh to be centered
    'units', 'normalized', 'Rotation',0,...
    'FontSize',fontSize);
set(gca, 'Layer', 'top', 'XGrid', 'on');
title('distribution of the read frequencies')
ylabel('P, %')
legend({'raw','filtered'})

%%

% figure  
% cm = jet(8);
% for chr = 1:5    
%     fv = cat(1,ChrReads(chr).f);
%     EMSqr = cat(1,ChrReads(chr).EMSqr);
%     hf = hist( fv, fx);
%     hfEMS = hist( fvEMS, fx);
%     stairs(fx', hf'./sum(hf), 'color', cm(chr,:),'linestyle', '-')
%     hold on
%     stairs(fx', hfEMS'./sum(hfEMS), 'color', cm(chr,:), 'linestyle', ':', 'linewidth',2)
%     title(sprintf('distribution of the read frequencies, chr#%u', chr))
%     ylabel('P, %')
%     xlabel('frequency of mutant reads, %')
%     legend({'raw','filtered'})
% end

rv = cat(1,ChrReads(:).r);
rvEMS = rv(EMSqr);

figure('name', 'scatter: f vs r ')
scatter(rv , fv, 5,'rx')
hold on
scatter(rvEMS , fvEMS, 5,'bo')
set(gca, 'xscale', 'log')
quntf = quantile(fv, 0.99); %%%%%%% <---
xlim(10.^[0 4])
xlabel('number of reads, \itr')
ylabel('frequency of mutant reads, \itf')
legend({'raw','filtered'})

dx = cellfun(@diff,{ChrReads(:).x},'un',0);
dx = single( cat(1, dx{:} ) );
EMSqr = cat(1,ChrReads(:).EMSqr);
dxEMS = dx(EMSqr);

ndx = 10.^(0:.1:5)';

pdx = hist(dx,  ndx );
pdxEMS = hist(dxEMS,  ndx );
[xb, yb] = stairs(ndx, [pdx', pdxEMS']);
% ndx = [1:9,10:10:90 ,100:100:900, 1000:1000:9000, 1e4:1e4:9e4,1e5 ];
% ndx = logspace(3, log10(max(dx)), 60);

figure('name' , 'distance')
myhistlog(ndx, dx, 'r', 'edgecolor', 'none')
% patch([xb(:,1); 1], [yb(:,1); 0], 'r', 'edgecolor', 'none')
hold on
% patch([xb(:,2); 1], [yb(:,2); 0], 'g', 'edgecolor', 'none')
myhistlog(ndx, dxEMS, 'g', 'edgecolor', 'none')
legend({'raw','filtered'})
set(gca, 'xscale', 'log', 'tickdir', 'out')

%%
xrv = 1:5001;
figure
hrv = hist(rv, xrv);
hrvEMS = hist(rvEMS, xrv);
plot(xrv, [hrv',hrvEMS'])
xlim([1 5000])
hold all
plot(xrv, length(rv)   * pdf('poisson', xrv, double(mean(rv))) ,'b:','linewidth',2 )
plot(xrv, length(rvEMS)* pdf('poisson', xrv, mean(rvEMS)),      'g:','linewidth',2 )



mean(rv)
mean(rvEMS)
% ah = area(ndx, [pdx', pdxEMS'-pdx']);

% legend({'filtered','raw'})

% export_fig images/AB1Coverage -r600 -eps
