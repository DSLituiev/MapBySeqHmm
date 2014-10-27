
addpath('D:\MATLABuserfunctions\NewtonRaphson');

bb = 0:0.5:50;


[Ltot, V, L] = likelihoodBetaBinomial(double(AR.q), double(AR.r), bb, 0.25);

figure;surf(bb, 1:numel(AR.r) ,L, 'linestyle', 'none')
view( 0, 90 )
%
figure;
plot(bb, Ltot)

t0 = 0:0.02:1;

[Ltot, V, L, phi] = likelihoodBetaBinomial(double(AR.q), double(AR.r), bb, permute(t0, [1,3,2]) );

figure;
surf( t0, bb,  squeeze(Ltot), 'linestyle', 'none')
view( 0, 90 )


figure;
surf( t0, bb,  squeeze(phi)./(1+squeeze(phi)), 'linestyle', 'none')
view( 0, 90 )

figure;
surf( t0, bb,  squeeze(nanmean(V,1)), 'linestyle', 'none')
view( 0, 90 )

set(gca, 'clim', [-4 0])


%%
pi01 = 0.99;
f1 = @(x,y,b,z)likelihoodBetaBinomial(x, y, b, 0.25, z);
f2 = @(x,y,b)(1/51);
b0 = 10;

for ii = 1:100
    [bNext, piNext1] = myEM( double(AR.q), double(AR.r), pi01, f1, f2, b0);
    piNext1Mean = mean(piNext1);
    fprintf('delta theta: %g \n', abs(b0 - bNext ))
    fprintf('theta: %g \t pi %g\n', bNext,  piNext1Mean)
    if abs(bNext-b0) < 1e-9
        break
    end
    b0 = bNext;
    pi01 = piNext1Mean;
end


