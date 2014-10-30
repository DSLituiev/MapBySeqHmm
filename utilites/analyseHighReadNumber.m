function analyseHighReadNumber( RR )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

nTop = 20;
[~, tt] = maxk(RR.r, nTop);


nonRepeat = strcmpi('NO', RR.repeat(tt) );

for ii = 1:numel(nTop)    
    fprintf('chr.\t%u\t position\t%u\t',  RR.maxHitGene{ii}, RR.chromosome(ii) )
    
    if nonRepeat(ii)
        fprintf('gene: \t%s\n',  RR.maxHitGene{  ii } )
    else
        fprintf('type:\t%s\t name:\t%s\n', RR.repeat{ii}, RR.repeatName{ii})
    end
    
end

