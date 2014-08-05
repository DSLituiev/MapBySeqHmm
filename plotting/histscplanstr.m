function histscplanstr(ChrReads, chr)
% histogram and scatter plot analysis for a chromosome data;
ChrTextTag = ['chromosome #',num2str(chr),'; '];
%% Histogram analysis for reads

rvect = single(min(ChrReads(chr).r) - 2:5:max(ChrReads(chr).r)-2);
fvect = 0.05:0.1:0.95;
figure('name',[ChrTextTag, 'Histogram analysis for reads'])
subplot(2,1,1)
hist(single(ChrReads(chr).r), rvect);
xlim([0, max(ChrReads(chr).r)]);
xlabel('total number of reads')
% title('reads')
subplot(2,1,2)
ChrReads(chr).f = single(ChrReads(chr).q)./single(ChrReads(chr).r);
binCtrs = min(ChrReads(chr).f):0.025:max(ChrReads(chr).f);
histf = hist(ChrReads(chr).f, binCtrs);
bar(binCtrs,histf,'hist');
h = get(gca,'child');
set(h,'FaceColor',[.9 .9 .9]);
%xlabel('Time to Failure'); ylabel('Probability Density'); ylim([0 0.1]);
% xgrid = linspace(0,20,100);
hold all;
r0 = 40;
theorbinom = binopdf(1:1:r0,r0,1/4);
plot((1:1:r0)./r0,sum(histf)./sum(theorbinom).*theorbinom,'r-','linewidth',2)
xlim([min(ChrReads(chr).f)-0.025, max(ChrReads(chr).f)]);
legend({'empiric','theoretical'})
xlabel('mutant SNP frequencies in reads')
%% 3d histogram
%     H3 = zeros(size(rvect,2),size(fvect,2));
%     for ri = 1:length(rvect)
%     H3(ri,:) = hist(ChrReads(chr).f(ChrReads(chr).r == ri), fvect);
%     end
%     [RRf,FFr] = meshgrid(rvect,fvect);
%
%     figure('name',[ChrTextTag, '3D Histogram: P(r,f)'])
%     pcolor(RRf,FFr,H3');
%     view(0,90); xlabel('r');ylabel('f'); axis square
%
% figure;bar3(H3','detached');
%% Scatter plot for the reads
figure('name',[ChrTextTag, 'Scatter plot analysis for reads'])
% scatter3(log10(single(ChrReads(chr).r)), ChrReads(chr).f, -ones(size(ChrReads(chr).f)), 5,'kx')
cloudPlot(log10(single(ChrReads(chr).r)), ChrReads(chr).f,[],...
    1,[30,30],'AlphaData', 0.8)
hold on
scatter(log10(single(ChrReads(chr).r)), ChrReads(chr).f, 5,'ko')
% scattercloud(single(ChrReads(chr).r), ChrReads(chr).f, 25, 0.7, 'kx',jet)
xlabel('r'); ylabel('f'); ylim([0,1])

figure('name',[ChrTextTag, '3D Scatter plot analysis for reads'])
% scatter3(log10(single(ChrReads(chr).r)), ChrReads(chr).f, -ones(size(ChrReads(chr).f)), 5,'kx')

scatter3(single(ChrReads(chr).x)*1e-6, ChrReads(chr).f, log10(single(ChrReads(chr).r)), 5, log10(single(ChrReads(chr).r)),'o')
% scattercloud(single(ChrReads(chr).r), ChrReads(chr).f, 25, 0.7, 'kx',jet)
xlabel('x'); ylabel('f'); zlabel('r'); ylim([0,1])
