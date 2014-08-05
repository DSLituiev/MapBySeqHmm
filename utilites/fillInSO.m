function outStr = fillInSO(terms, inStr, outStr, xind)

for ii = 1:numel(terms)
    effects = regexp(inStr.effect,'&', 'split');
    for jj = 1:numel(effects)
        inds = strcmp(terms{ii}(4:end), effects{jj});
        if any(inds)
            outStr.(terms{ii}){xind,1} = inStr.gene;
        end
    end
    
    if ~isempty(inStr.CDS)
         outStr.(terms{ii}){xind,2} = uint32(sscanf(inStr.CDS, '%u'));
         outStr.(terms{ii}){xind,3} = inStr.AA;
         outStr.(terms{ii}){xind,4} = inStr.Codone;
    end
end
