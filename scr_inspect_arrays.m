AR.xLogOdds = AR.xPsel - AR.xPflat

AR.plotChromosomes('xLogOdds', 'yscale', 'lin', 'figure', 'new', 'plotfun', @(x,y)plot(x,y, 'r-'));


AR.xLogOdds = AR.xPstat - AR.xPflat

AR.plotChromosomes('xPflat', 'yscale', 'lin', 'figure', 'new', 'plotfun', @(x,y)plot(x,y, 'r-'));

figure
p = pcolor( AR.x(AR.chromosome == 5), AR.pop.kvect,  10.^AR.xkPflat(AR.chromosome == 5,:)' );
set(p, 'linestyle', 'none')
set(gca, 'ydir', 'normal')
hold all
plot( AR.x(AR.chromosome == 5), AR.fmFiltered(AR.chromosome == 5,:).*AR.pop.N*2 )


figure
plot( AR.x(AR.chromosome == 5), calcMarginal( AR.xkPflat(AR.chromosome == 5, 1:end),2) )


figure
p = pcolor( AR.x(AR.chromosome == 5), AR.pop.kvect,  log10(AR.E(AR.chromosome == 5,:)') );
set(p, 'linestyle', 'none')
set(gca, 'ydir', 'normal')
view(0,90)
hold all
plot( AR.x(AR.chromosome == 5), AR.fmFiltered(AR.chromosome == 5,:).*AR.pop.N*2)

figure
p = pcolor( AR.x(AR.chromosome == 5), AR.pop.kvect,  (AR.E(AR.chromosome == 5,:)') );
set(p, 'linestyle', 'none')
set(gca, 'ydir', 'normal')
view(0,90)
hold all
plot3( AR.x(AR.chromosome == 5), AR.fmFiltered(AR.chromosome == 5,:).*AR.pop.N*2, 20*ones(sum(AR.chromosome == 5)), 'r-' )


figure
plot( AR.x(AR.chromosome == 5),  log10( sum(AR.E(AR.chromosome == 5,:)')) )

% 
% 
% figure
% plot( obj.x(obj.chromosome == 5),  log10( sum(obj.E(obj.chromosome == 5,:)')) )
