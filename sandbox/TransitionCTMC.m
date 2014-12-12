addpath('..')
clear all

clear T
pop = population(10);

tau = 1;
theta0 = 1;

T0 = expm(tau * theta0 * pop.Q);

theta = theta0 + permute(2.^-(0:1:10), [3, 1, 2]);

T = zeros(pop.Np, pop.Np, numel(theta));

for ii = 1:numel(theta)
    T(:,:, ii) = expm(tau * theta(ii) * pop.Q);
end

dT =  bsxfun(@minus, T, T0);
dT_dtheta = bsxfun(@rdivide, dT, theta);

n = squeeze(sum(sum( bsxfun(@minus, dT_dtheta, bsxfun(@times, pop.Q, tau)).^2, 1),2)) ;

figure
plot(squeeze(theta-theta0), n)
set(gca, 'xscale', 'log', 'yscale', 'log')


figure
pcolor(pop.Q)

tt = numel(theta-3);
figure
pcolor(dT_dtheta(:,:,tt))


figure
pcolor(T(:,:,tt))