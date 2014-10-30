function outStr = initializeSO(terms, outStr, xind)

for ii = 1:numel(terms)
    outStr.(terms{ii}){xind,1} = [];
end
