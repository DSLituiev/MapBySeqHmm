function TransFreqMatr(sumStates, t0, Zmatr)
%% calculate matrix of the transition frequencies based on the simulation

TransitionFreqMatr = sum(uint8(bsxfun(@and,permute(bsxfun(@eq,sumStates(t0:end-1),0:size(Zmatr,1)),[3,2,1]),...
    permute(bsxfun(@eq,sumStates(t0+1:end),0:size(Zmatr,1)),[2,3,1]))) ,3);

filt = fspecial('gaussian',5,2);
PFlsmooth = conv2(TransitionFreqMatr, filt,'same');

figure('name','Joint')
surf(TransitionFreqMatr,'linestyle','none'); view(0,90)
axis tight square
sum(sum(PFlsmooth-PFlsmooth'))
% 
% figure('name','Joint')
% surf(log10(PFl),'linestyle','none'); view(0,90)
% axis tight square
% sum(sum(PFlsmooth-PFlsmooth'))


figure('name','Joint, smooth')
surf(PFlsmooth,'linestyle','none'); view(0,90)
axis tight square

Pst = sum(uint8(bsxfun(@eq,sumStates(t0:end-1),0:size(Zmatr,1))),1);
Pstsmooth = smooth(Pst,5);
T = bsxfun(@rdivide,TransitionFreqMatr,Pst');
T(isnan(T)) = 0;
Tsmooth = bsxfun(@rdivide,PFlsmooth,Pstsmooth);
TsmoothT = bsxfun(@rdivide,PFlsmooth,Pstsmooth');

Tsmooth(isnan(Tsmooth)) = 0;
TsmoothT(isnan(TsmoothT)) = 0;


figure('name','T')
surf( conv2(T,filt,'same'),'linestyle','none'); view(0,90)
axis tight square

figure('name','Tsmooth')
surf(Tsmooth,'linestyle','none'); view(0,90)
axis tight square
hold on
plot(0:1:AveragingTimes,0:1:AveragingTimes,'k')

figure('name','Stationary distribution')
plot(0:1:AveragingTimes,Pst)
xlim([0,AveragingTimes])
% colormap(jet)