function CalcSetLegpos(spl, leg)
% calculates and sets the legend position
splpos = cell2mat(get(spl, 'position'));
[~, chrm] = min(splpos(:,3));
%= left bottom width height]
legpos = [ splpos(chrm,1) + splpos(chrm,3)+0.05, splpos(chrm,2),...
    1- splpos(chrm, 3) - 2*splpos(chrm,1),  splpos(chrm,4) ];
set(leg, 'position', legpos )