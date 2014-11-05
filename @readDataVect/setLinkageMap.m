function setLinkageMap(obj,  mapPath)
%=       load the recombination map (contains positions of the markers
%= and genetic distance in cM between them)
load(mapPath)
obj.chrMap = ChrMap;
clear ChrMap;
