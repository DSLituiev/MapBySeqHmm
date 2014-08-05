clear all ; clc
addpath('D:\MATLABuserfunctions\NewtonRaphson');
addpath('D:\MATLABuserfunctions\binomial');
addpath('D:\MATLABuserfunctions');
addpath('.\betabinomial');

n0 = 100;
k = (0:1:n0)';
n = n0*ones(size(k));
% 
% k = 10;
% n = 40;

mu = .25;

theta = 0.125;

b = (1-mu)./theta;
a = mu./theta;
 
% Var =  bsxfun(@rdivide, ...
%     bsxfun(@times, n0, alpha.*beta).*(bsxfun(@plus, (alpha+beta), n0) ),...
%     ( (alpha+beta).^2 .* (alpha+beta+1) ) );



Lb = likelihoodBetaBinomial(k, n, b, mu);

[ Res0,  Residue, LL ] = maxLlThetaBetaBinomial(k, n, mu, theta, true);


[ LL] = logBetaBinomialThetaMu0(k, n, mu, theta)

% sum(Lb.*(k - sum(k.*Lb)).^2)
% sum(LL.*(k - sum(k.*LL)).^2)

figure
plot(k, Lb)
hold all
plot(k, 10.^LL, '--')
% plot(k, L0, ':')
plot(  sum(k.*Lb), 0, 'bx')
plot(  sum(k.*10.^LL), 0.001, 'gx')
% plot(  sum(k.*L0), 0.002, 'rx')


figure
plot(k, Res0)

10.^(LL0)


prod(Lb)
10.^(L)


%%
mu = .25;

k = n0.*mu + 5;newtonraphson
n = n0;
% 
% k = 10;
% n = 40;

% 

k = 100*round(rand(1,200));
n = 3.4*k+ round(40*rand(size(k)));% 

% mEmp = nanmean(k./n);
% sEmp = nanvar(k./n);

theta = (0:0.005:2)';

b = (1-mu)./theta;
a = mu./theta;

[ L, logLtot, Var, phi, Lb0   ] = likelihoodBetaBinomial(k, n, b, mu);

[sumdLdTheta,  Residue, LL ] = maxLlThetaBetaBinomial(k, n, mu, theta, true);

figure
% plot(theta, L, '-')
hold all
plot(theta, 10.^LL, ':')
plot(theta, 10.^mean(LL,2),  'k-', 'linewidth', 2)
% plot(  sum(bsxfun(@times,theta, 10.^LL)), 0.001, 'x')

figure('name', 'dlnL/dTheta')
plot(theta, abs(Residue))
hold all
plot(theta, abs(sumdLdTheta), 'k-', 'linewidth', 2)
% figure('name', 'dlnL/dTheta')
% plot(theta(2:end)-diff(theta(1:2))/2, abs( bsxfun(@rdivide, diff(LL,1,1), diff(theta) ) ), ':')

plot(theta(2:end)-diff(theta(1:2))/2,...
    abs( bsxfun(@rdivide, diff(sum(LL,2), 1), diff(theta) ) ), 'r-', 'linewidth', 2)
set(gca, 'yscale', 'log')

options = optimset('TolX',1e-10); % set TolX
[theta1Next, ~, ~] = newtonraphson(@(x) maxLlThetaBetaBinomial(k, n, mu, x), [0.5], options);

fitBetaBinomialStat(k', n', 0.25)


fitBetaBinomialStat(k', n', 0.25, 0.5719)


thetaEst = fitBetaBinomialStat( single(AR.q), single(AR.r), 0.25 )