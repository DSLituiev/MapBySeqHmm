function f_out = visualizeStat(obj)
% VISUALIZESTAT  -- plots descriptive statistics on the 'variant called' data
%
% ## Output: a figure handle

%% initialize plotting parameters:
% nBins = max( floor(numel(obj.x)/25), 7 );
%==
df = 0.025;  % 0.01*ceil(100/nBins);
fx = df/2:df:(1-df/2);
%==
dr = 0.05;
rx = 0:dr:(log10(max(obj.r))-dr/2);
%==
dlogx = 0.2;
dlogx_x = 0:dlogx:(log10(max(obj.dx))-dlogx/2);
%%
f_out = figure('name', 'descriptive statistics');
%%
subplot(2,2,1)
myhist100(fx,  obj.f, 'r', 'edgecolor', 'none');
set(gca, 'tickDir', 'out')
title('SNP ratio [f]')
xlabel('f');
%%
subplot(2,2,2)
[~, modeX] = myhist100LogX(rx,  log10(obj.r), 'r', 'edgecolor', 'none');
set(gca, 'tickDir', 'out')
title( sprintf('coverage / read number per locus [r], \n median[r] = %u, mode[r] = %u', floor(median(obj.r)), floor(modeX)) );
xlabel('r [log_{10}-scaled]');
%%
subplot(2,2,3)
scatter(log10(obj.r), obj.f, 'r');
xlabel('log_{10}(r)'); ylabel('f')
set(gca, 'tickDir', 'out')
title('read number [r] vs SNP ratio [f]')
set(gca, 'ylim', [0,1])
%%
subplot(2,2,4)
myhistLogX(dlogx_x , log10(obj.dx), 'r')
set(gca, 'tickDir', 'out', 'xscale', 'log')
title('nt to closest neighbour SNP [${\Delta} x$]', 'interpreter', 'latex')
xlabel('${\Delta} x$ [log$_{10}$-scaled]', 'interpreter', 'latex');
end