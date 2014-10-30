function histscplan(Dataxrq,MinM,MaxM,chr)
% histogram and scatter plot analysis for a chromosome data;
 rx = single(Dataxrq(:,2)); 
 qx = single(Dataxrq(:,3));
  ChrTextTag = ['chromosome #',num2str(chr),'; '];
%% Histogram analysis for reads
    rvect = single(MinM(chr,2)-2:5:MaxM(chr,2)-2);
    fvect = 0.05:0.1:0.95;
    figure('name',[ChrTextTag, 'Histogram analysis for reads'])
       subplot(2,1,1)
    hist(rx,rvect);
    xlim([0, MaxM(chr,2)]);
    xlabel('total number of reads')
    % title('reads')
    subplot(2,1,2)
    fx = qx./rx;
    binCtrs = min(fx):0.025:max(fx);
    myhi = hist(fx,binCtrs);
    bar(binCtrs,myhi,'hist');
    h = get(gca,'child');
    set(h,'FaceColor',[.9 .9 .9]);
    %xlabel('Time to Failure'); ylabel('Probability Density'); ylim([0 0.1]);
    % xgrid = linspace(0,20,100);
    hold all;
    r0 = 40;
    theorbinom = binopdf(1:1:r0,r0,1/4);
    plot((1:1:r0)./r0,sum(myhi)./sum(theorbinom).*theorbinom,'r-')
    xlim([min(fx)-0.025, max(fx)]);
    xlabel('mutant SNP frequencies in reads')
    %% 3d histogram
    H3 = zeros(size(rvect,2),size(fvect,2));
    for ri = 1:length(rvect)
    H3(ri,:) = hist(fx(rx == ri),fvect);
    end
    [RRf,FFr] = meshgrid(rvect,fvect);
    figure('name',[ChrTextTag, '3D Histogram: P(r,f)'])
    pcolor(RRf,FFr,H3');% ,'LineStyle', 'none')
   % set(gca,'XTick',MbLine,'XTickLabel',ticklab,...
   %     'YTick',MbLine,'YTickLabel',ticklab,'TickDir','out')
    view(0,90);% xlabel('x');ylabel('y'); axis square
    % figure;bar3(H3','detached');
 %% Scatter plot for the reads   
     figure('name',[ChrTextTag, 'Scatter plot analysis for reads'])        
     scatter(rx,fx)  
     xlabel('r'); ylabel('f')
