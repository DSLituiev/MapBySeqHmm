
mu = 10.^(0:0.5:6);
t = 10.^(0:0.5:6)';

z = bsxfun(@times, mu, exp(-bsxfun(@times,mu,t)));

figure
surf(t,mu,z)
set(gca, 'yscale', 'log', 'xscale', 'log', 'zscale', 'log')

t = 0:0.1:4;
figure
plot(t, 1-exp(-t/2))


M = 2000;
G = 5000;
G/M
mu = M/G;
tu = G*sort(rand(M,1));

d2x = (bsxfun(@minus, tu, tu')).^2 ;
hits = any(d2x < 1 - eye(M));

sum(hits)
dxRaw = diff(tu(hits));
dxPair = [dxRaw, circshift(dxRaw,-1)];
dx = nanmin(dxPair,[],2);




figure
hist( log10(diff(tu)), 100 )
hold all
plot( log10(G/M*[1,1]), [0,1], 'rx')


y = (-3:0.05:2)';
t = 10.^y;
plogt = log(10).*t.*mu.*exp(-mu.*t);

figure
plot(y, plogt, 'g:')
hold on
plot( log10(G/M)*[1,1], [0,1], 'r-x')

%% for miss-mapping

y = (-8:0.2:1)';
t = 2.^y;

N = ceil(2.^(0:1:15));
T = 1;
P = bsxfun(@times,  N + 1, bsxfun(@times, t/T, bsxfun(@power, (1 - t./T), N)) );
P( t./T > 1, :) = 0;
 
figure
plot(t, P)
set(gca, 'xscale', 'log')
sum(P, 2)

figure
surf( t, N, P')
set(gca, 'xscale', 'log', 'yscale', 'log')
view(0, 90)
ylabel('N')
xlabel('\Delta{}x')

hold all
plot3(T./(N+1), N, ones(size(N)), 'r-', 'linewidth', 2)
