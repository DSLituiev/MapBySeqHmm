AR.pop = N;

r = 25;
q = (0:r)';
f = 0.25;

E = emissionk0(q, repmat(r, r+1, 1), AR.pop);
% p = binopdf(q,r, AR.pop.f);
p = AR.pop.Pstat*E';

figure

AR.filterFields('r', @(x)(x>=r)); % mutant reads

df = 0.05;  % 0.01*ceil(100/nBins);
fx = df/2:df:(1-df/2);
fy = AR.f(AR.chromosome~=1);

figure
myhist100(fx,  fy, 'r', 'edgecolor', 'none');
set(gca, 'tickDir', 'out')
xlabel('SNP ratio [\it f \rm]')
ylabel('PDF')

hold all

plot(q./r, p, 'linewidth', 2)

legend({'empirical PDF, $\;\;r \ge 25$', 'theoretical PDF, $r = 25$'}, 'interpreter', 'latex')
xlim([0,1])
ylim([0,0.22])
fig(gcf, 'width', 16)


exportfig(gcf, 'overdispersion', 'format','eps', 'color', 'rgb')


