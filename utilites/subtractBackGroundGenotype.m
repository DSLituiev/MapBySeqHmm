function [testReads, annotation] = subtractBackGroundGenotype(testPath, varargin)
%% check the input parameters
p = inputParser;
addRequired(p, 'testPath', @(x)(ischar(x) & exist(x, 'file')) );
addOptional(p, 'bkgrPath', '', @(x)(ischar(x) & exist(x, 'file')));
 addParamValue(p,     'version',          2, @isscalar);
parse(p, testPath, varargin{:});

%%
if nargin>1 && ~isempty(p.Results.bkgrPath)
    refPath = p.Results.bkgrPath;
    
    if p.Results.version == 1
    [testReads, ~,  annotation] = readSequencingDataCsv(testPath);
    
    [refReads, ~] = readSequencingDataCsv(refPath, 'noannotation');
    elseif p.Results.version == 2
        
    [testReads, ~,  annotation] = readSequencingDataCsv2(testPath);
    
    [refReads, ~] = readSequencingDataCsv2(refPath, 'noannotation');
    end
    
    testPositions = [testReads.chromosome,  testReads.x];
    
    refPositions = [refReads.chromosome,  refReads.x];
    
    
    [~, sCommon]  = intersect(testPositions, refPositions, 'rows');
    iUnique = true( size(testPositions,1), 1);
    
    iUnique(sCommon) = false;
    
    fprintf('%u out of %u reads left\n', sum(iUnique),  size(testPositions,1) );
    
    testReads = selectFieldIndices( testReads , iUnique);
    annotation =  selectFieldIndices( annotation, iUnique);
else
    if p.Results.version == 1
        [testReads, ~,  annotation] = readSequencingDataCsv(testPath);
    elseif p.Results.version == 2
        [testReads, ~,  annotation] = readSequencingDataCsv2(testPath);
    end
end

end