
clear Lst Lfl Lsl
Alpha = 10.^(-2:0.5:3); % 2.^(-4:0.5:7)';
% Alpha = 1:5:101;
for cc = 1:5
    obj = AR.HMM{cc};
    obj.resetFlag = true;
    
    for aa = numel(Alpha):-1:1
        Lst(aa,cc) = obj.statLikelihood(Alpha(aa));
        Lsl(aa,cc) = obj.selLikelihood( Alpha(aa));
        Lfl(aa,cc) = obj.totLikelihood( Alpha(aa));
    end
end

for cc = 1:5
    figure
    plot(Alpha, [Lst(:,cc)- Lfl(:,cc), Lsl(:,cc)- Lfl(:,cc)])
    set(gca, 'xscale', 'log')
end


for cc = 1:5
    figure
    plot(Alpha(2:end), [diff(Lst(:,cc)), diff(Lsl(:,cc)), diff(Lfl(:,cc))])
    set(gca, 'xscale', 'log')
end


[mL, ind] = max(Lsl);
Alpha(ind)
%
% medL = median(Lsl)
% log10(sum( Alpha.*10.^(Lsl' - medL) )) - medL
