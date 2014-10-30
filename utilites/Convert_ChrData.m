% convert ChrData file to a structure
load ChrData
for ii = 1:length(ChrData)
    ChrMap(ii).nt = ChrData{ii}(:,1);
    ChrMap(ii).cM = ChrData{ii}(:,2);
end

save ChrMap