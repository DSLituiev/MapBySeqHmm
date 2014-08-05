% convert ChrData file to a structure
clear all
ftag = 'HL7';
load(['..\datafiles\Data',...
    ftag, 'uint32cell']);

for ii = 1:length(Data_xrqExNsSt)
    ChrReads(ii).x = Data_xrqExNsSt{1,ii}(:,1);
    ChrReads(ii).r = Data_xrqExNsSt{1,ii}(:,2);
    ChrReads(ii).q = Data_xrqExNsSt{1,ii}(:,3);
 %   ChrReads(ii).PriorExNsSt = Data_xrqExNsSt{2,ii};
    ChrReads(ii).Pr.snpExon = Data_xrqExNsSt{2,ii}(:,1);
    ChrReads(ii).Pr.snpNonsyn = Data_xrqExNsSt{2,ii}(:,2);
    ChrReads(ii).Pr.snpStop = Data_xrqExNsSt{2,ii}(:,3);
end
% clear Data_xrqExNsSt  ii ans chr ChrNumber 
save(['..\datafiles\ChrReads_',...
    ftag], 'ChrReads');
