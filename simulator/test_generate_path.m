%==== A simulator script for the mapping process
close all; clear all; clc; %<= magic spells

%% add useful functions
addpath('D:\MATLABuserfunctions\mtimesx');
addpath('D:\MATLABuserfunctions\MinMaxSelection');
addpath('D:\MATLABuserfunctions');
%= add the parent path
addpath(cd(cd('..')))
addpath('../emission')

%% == Initialize parameters:
LINKAGE_LOSENING = 1;
xCausativeSNP = 0.2;
T_max = 7;

par = simPooledSeq(xCausativeSNP, 50, 12, 100, T_max);
par.Selection = false;
%% == generate a path 
N_iter = 700;
for ii = N_iter:-1:1
    [tChr{ii}, kChr{ii}] = par.generatePath();
end

figure
for ii = N_iter:-1:1
    plot(tChr{ii}, kChr{ii}, '-', 'color', [1,1,1]*0.75)
    hold on
end
plot([0, T_max], par.pop.N/2* [1, 1], 'k:')
ylim([0,par.pop.N])

%%
t0 = 0;

for ii = N_iter:-1:1
    k_t0(ii) = kChr{ii}( find(tChr{ii} > t0, 1, 'first') );
end



[freq] = hist(k_t0, par.pop.kvect);
epdf = freq./sum(freq);

tpdf = par.pop.Pstat;

figure
plot(tpdf.^freq)


figure; 
bar(par.pop.kvect, epdf, 'linestyle', 'none' )
hold on
plot(par.pop.kvect, tpdf, 'r-', 'linewidth', 2)
xlim([0,par.pop.N])

chi2 = sum( ( freq - tpdf*N_iter).^2 ./ (tpdf*N_iter) )

y = chi2pdf(chi2, par.pop.N)


