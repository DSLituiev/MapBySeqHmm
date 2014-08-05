% f = figure('name','Sample path');
% stairs(timevals(1:t0),sumStates(1:t0), 'r-','linewidth',3)
% hold all
% arbtp = [1; sort(1+round((t0-1)*rand(12,1))); t0];
% plot(timevals(arbtp),sumStates(arbtp), 'bo-','linewidth',2)
% set(gca,'xtick',[])
% 
% 
% N= 50;
% for kk = 0:N
% Pst(kk+1) = (2^(-N))*nchoosek(N,kk);
% end
% 
% f= figure;
% bar(0:N, Pst)
% 
% export_fig P50  -r600 -eps  -painters
%% General params
Npl = 10;
Copies = 20;
Perc = 0.4;
%% N plants

[X Y] = meshgrid(1:2,1:(Npl/2));
CU = logical(Y<=Perc*Npl).*logical(X<2);
f = figure;
% scatter(X(:),Y(:),200,Cv,'filled')
scatter(X(:),Y(:),200,CU(:),'filled')
colormap(hsv)
caxis([-0.5 1])
set(gca,'xlim',[0 2*Copies+1],'ylim',[0 Npl+1], 'xtick',[], 'ytick',[],...
    'xcolor','w','ycolor','w')
fig(f,'width',2*Copies,'height',Npl)
export_fig PlantsFr -eps -painters 

%% ReadsOrd
[X Y] = meshgrid(1:Copies,1:Npl);
CU = logical(Y<=Perc*Npl);
% Cv = bsxfun(@times,CU(:),[0.498 1 0]) + bsxfun(@times,~CU(:),[1 0 0]);
f = figure;
% scatter(X(:),Y(:),200,Cv,'filled')
scatter(X(:),Y(:),200,CU(:),'filled')
colormap(hsv)
caxis([-0.5 1])
set(gca,'xlim',[0 Copies+1],'ylim',[0 Npl+1], 'xtick',[], 'ytick',[],...
    'xcolor','w','ycolor','w')
fig(f,'width',Copies,'height',Npl)
export_fig ReadsOrd   -eps -painters 
% exportfig(f,'ReadsOrd.eps','Format','eps','Color' ,'rgb','height',10);

%% ReadsMix
randseq = randperm((Copies*Npl));
% C = logical(round(rand(Copies,Npl)./Perc));
C = reshape(CU(randseq),Npl,Copies);

%Cv = bsxfun(@times,C(:),[0.498 1 0]) + bsxfun(@times,~C(:),[1 0 0]);

%% ReadsMixSub
f = figure;
scatter(X(:),Y(:),200,C(:),'filled')
%scatter(X(:),Y(:),200,C(:),'filled')
set(gca,'xlim',[0 Copies+1],'ylim',[0 Npl+1], 'xtick',[], 'ytick',[],...
    'xcolor','w','ycolor','w')
colormap(hsv)
caxis([-0.5 1])
fig(f,'width',Copies,'height',Npl)
export_fig ReadsMix  -r600 -eps  -painters

S = logical((X(:)-Copies/2).^2+(Y(:)-Npl/2).^2<=3.6);
f = figure;
scatter(X(:),Y(:),200*(S+0.001),C(:),'filled')
set(gca,'xlim',[0 Copies+1],'ylim',[0 Npl+1], 'xtick',[], 'ytick',[],...
    'xcolor','w','ycolor','w')
colormap(hsv)
caxis([-0.5 1])
fig(f,'width',Copies,'height',Npl)
export_fig ReadsMixSub  -r600 -eps  -painters


sum(C(:))./length(C(:))
sum(S)
sum(C(:).*S)


