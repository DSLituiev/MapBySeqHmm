close all; clear all
addpath('C:\MATLABuserfunctions');
addpath('C:\MATLABuserfunctions\cm_and_cb_utilities');
% q = fliplr([3 6 12 30 300]);
% r = fliplr([4 8 16 40 400]);
q = ([1 2 4  10 100]);
r = ([4 8 16 40 400]);
N=50;
k = linspace(0,N,N+1);
[K,Q] = ndgrid(k,q);
[~,R] = ndgrid(k,r);

L = binopdf(Q,R,K./N);
L = bsxfun(@rdivide,L,sum(L,1));

%for ii = 1:length(q)
KK = zeros(size(K,1)*2,size(K,2));
KK(1:2:end) = K-0.5;
KK(2:2:end) = K+0.5;
LL = zeros(size(L,1)*2,size(L,2));
LL(1:2:end) = L;
LL(2:2:end) = L;
figf=figure;
% plot(F,L, 'linewidth',3)
plot(KK,LL, 'linewidth',3)
set(gca,'tickdir','out', 'xlim',[0,N])
% stairs(N*F,L, 'linewidth',3);
Leg = {'{\itq} = 1 / {\itr} = 4','     2 / 8','     4 / 16','    10 / 40','   100 / 400'};
% legend(Leg(end:-1:1))
cm = hsv;
cmlines(gca,cm(1:round(end/2),:))
legend(gca,Leg, 'FontSize',16, 'FontName','Courier' )
fig(figf,'width',20,'height',12)
legend('boxoff')
export_fig LikelihoodBin -eps -painters 