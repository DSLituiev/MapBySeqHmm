
addpath('D:\MATLABuserfunctions\mtimesx');

N = 100;
kvect = 0:1:N;
tvect = 0:0.01:1;
s = 4;

100*2^-s

T = transition3D_expm(N, tvect);
gxs = permute(  ((1+exp(-2*tvect))./2).^(s-1), [1,3,2]);

squeeze(1-gxs)

Pst = StationaryDistr(N);
Psel = zeros(N+1,1);
Psel(end) = 1;
Pflat = ones(N+1,1)./(N+1);

Tcis = zeros(N+1);
Tcis(:,1) = 1;

Tf = bsxfun(@times, T, gxs) + bsxfun(@times, Tcis, 1-gxs);

sum(Tf,2)
% 
% Px = mtimesx(Pst', T);
% Px = mtimesx(Pst', Tf);

Px = mtimesx(Psel', Tf);


figure
pcolor( tvect, kvect, squeeze(Px))
hold all
plot( tvect,  squeeze(sum(bsxfun(@times, Px, kvect),2)), 'g-', 'linewidth', 2)
plot( tvect,  squeeze(sum(bsxfun(@times, Px, kvect),2)), 'r-', 'linewidth', 1)

