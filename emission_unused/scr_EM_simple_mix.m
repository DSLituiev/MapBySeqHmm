close all; clear all; clc;
dbclear if warning
tic

addpath('D:\MATLABuserfunctions\binomial');
addpath('D:\MATLABuserfunctions\mtimesx');
addpath('D:\MATLABuserfunctions\MinMaxSelection');
addpath('D:\MATLABuserfunctions\newtonraphson');
addpath('D:\MATLABuserfunctions'); savepath;
addpath('..\utilites');
addpath('..\betabinomial');
addpath('..');

%=       provide the reference ID, which is the name of the csv file without '.csv'
%= together with the backround genotype ID
%= known positions of the causative SNP can be also provided here for
%= further visualization
% dataID = 'ABD159-rmdup-clipOverlap-q20-freebayes'; chr0 = 2; x0 = 17521246; % bkgrID = 'ABD241-rmdup-clipOverlap-freebayes';
% dataID = 'ABD173-rmdup-clipOverlap-freebayes'; chr0 = 3 ; x0 = 1619248; % bkgrID = 'ABD241-rmdup-clipOverlap-q20-freebayes';


% dataID = 'HL10-rmdup-clipOverlap-q20-freebayes';  chr0 =  3 ;  x0 =  16473265;

 dataID = 'HL7_Paired-rmdup-clipOverlap-freebayes'; x0 = 5672441; chr0 = 1;
%  dataID = 'HL7_Paired-rmdup-clipOverlap-q20-freebayes'; x0 = 5672441; chr0 = 1;
% dataID = 'HL7_Paired-ngm-rmdup-clipOverlap-freebayes'; x0 = 5672441; chr0 = 1;
% dataID = 'HL7_Paired-rmdup-clipOverlap-mpileup'; x0 = 5672441; chr0 = 1;

disp(['=======  Processing data from the run ''', dataID, ''' ======='])

% mkdir(fullfile('figures',dataID))
%=       load the recombination map (contains positions of the markers
%= and genetic distance in cM between them)

%= construct the path to the (primary experimental) data file
dataPath = fullfile('../data', [dataID, '-ems-annotation-repfilt.csv'] );
%= extract the refenece reads if the reference ID is given:
[AR, annotation] = subtractBackGroundGenotype(dataPath);

AR = calcDx(AR,  @(x)nanmin(x, [], 2));
[AR, mu, iT, ~] = unmixRepeatsOne( AR, 'dx', '','modeNum', 2);
%%

% [AR, mu, ~, mixtObj, fh] = unmixRepeatsOne( AR, 'dx', 'plot','modeNum', 2);
% fig(fh(3), 'height', 15, 'width', 12)
% exportfig(fh(3), fullfile('../figures', [dataID, '-scatter']), 'format','eps', 'color', 'rgb')
% 
% fig(fh(4), 'height', 12, 'width', 12)
% exportfig(fh(4), fullfile('../figures', [dataID, '-dx_hist']), 'format','eps', 'color', 'rgb')
% 
% figure('name','true')
% hist( log10( double(AR.r(AR.contrib>0.5)) ), 20)
% 
% figure
% hist(log10(double(AR.r(AR.contrib<0.5)) ) , 20)
% 
% 
% 
% figure('name','true')
% [yh,xh] = hist( (AR.f(AR.r>50)), 20);
% bar(xh,yh./sum(yh))
% hold on
% plot([1,1]*0.5, [0,1]*0.3, 'r-')
% xlim([0,1]); ylim([0,1]*0.3)
% fig(gcf, 'height', 8, 'width', 12)
% exportfig(gcf, fullfile('../figures', [dataID, '-f_hist']), 'format','eps', 'color', 'rgb')
% 
% return
%%
inds =  ( AR.chromosome~=chr0); %
q = double(AR.q(inds));
r = double(AR.r(inds));

N = 50;

% 
% [theta, lambda1, gi1] = runSimpleEM_BetaBinomAndUniform(q,r, N, 'v',...
%      'errTol', 1e-6);

 [theta, lambda1, gi1] = runSimpleEM_BetaBinomAndUniform(q,r, N, 'v',...
     'contribution', AR.contrib, 'errTol', 1e-6);

[thetSimple, ~, ~] = runSimpleEM_BetaBinomAndUniform(q,r, N, 's',...
    'contribution', AR.contrib, 'errTol', 1e-6);


% [theta, lambda1, gi1] = runSimpleEM_BetaBinomAndUniform(q,r, N, 'v',...
%     'errTol', 1e-6);

figure
plot(gi1);


inds = AR.r>20;
[yh, xh] = hist(AR.f(inds), 30);
[yhT] = hist(AR.f(inds& AR.contrib > 0.5), xh);
[yhF] = hist(AR.f(inds& AR.contrib < 0.5), xh);

figure
barstairs(xh, yh./sum(yh), 'r')
hold on
b = barstairs(xh, yhT./sum(yhT), 'g'); alpha(b, 0.4)
ylim([0,1]*0.25)
xlim([0,1])




f = 0:(1/2/N):(1/2);
Pf = StationaryDistr(N);

[pBB, pU] = conditBetaBinomStat(q, r, theta, N, f, Pf');

gi1 = lambda1*pBB./(lambda1*pBB + (1-lambda1)*pU );

figure
surf( 1:numel(q),f, bsxfun(@plus, (1-gi1)/N , bsxfun(@times, gi1, 10.^logBetaBinomialThetaMu0(q,r,f,theta)))', 'linestyle', 'none');
view(0, 90)
zlim([0, 0.5]); % clim([0, 0.5])


% ff =  0:(1/2/N):(1);
q1 = 25;
q2 = 60;
figure
plot(f, (binopdf(25,100,f)) , 'g-' , 'linewidth', 2)
hold all
plot(f, 10.^logBetaBinomialThetaMu0(q1 ,100,f, theta), 'm-', 'linewidth', 2)
plot(f, (1-lambda1)*1/N + lambda1*10.^logBetaBinomialThetaMu0(q1 ,100,f, theta), 'r-', 'linewidth', 2)
plot(f, (binopdf(q2 ,100,f)) , 'g:', 'linewidth', 2)
plot(f, 10.^logBetaBinomialThetaMu0(q2 ,100,f, theta) , 'm:', 'linewidth', 2)
plot(f, (1-lambda1)*1/N + lambda1*10.^logBetaBinomialThetaMu0(q2,100,f, theta), 'r:', 'linewidth', 2)
legend({'binomial', '\beta-binomial', '\beta-binomial + uniform'})
fig(gcf, 'height', 12, 'width', 16)

% exportfig(gcf, fullfile('../figures', ['demo-beta-mix']), 'format','eps', 'color', 'rgb')
