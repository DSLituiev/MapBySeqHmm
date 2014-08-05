function outStr = cleanSO(terms, outStr)

for ii = 1:numel(terms)
    if all( cellfun(@isempty, outStr.(terms{ii})) )
        outStr = rmfield(outStr, terms{ii} );
    end
end
