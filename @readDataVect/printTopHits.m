function printTopHits(obj, filename, varargin)

%% check the input parameters
p = inputParser;
addRequired(p, 'obj', @isobject);
addRequired(p, 'filename', @ischar);
%
addOptional(p,'numPerChr', 100, @(x)isnumeric(x)&isscalar(x)  );
addOptional(p,'rankField', 'xPosteriorNorm', @(x)isprop(obj, x)  );
%
addParamValue(p, 'cutoffValue',           [],  @(x)(isnumeric(x) & (isscalar(x)|isempty(x)))  )
addParamValue(p, 'cutoffField', 'xPselNorm', @(x)isprop(obj, x) );

parse(p, obj, filename, varargin{:});

numPerChr = p.Results.numPerChr;
rankField = p.Results.rankField;
%%
soDict = soDictionary();

if ~isempty(p.Results.cutoffValue)
    hitInds = (p.Results.cutoffValue < obj.(rankField) );
else
    hitInds = (obj.(p.Results.cutoffField) < obj.(rankField) );    
end

for  chr = obj.chrNumber:-1:1
    a =  max(obj.(rankField)(hitInds & obj.chromosome == chr));
    if ~isempty(a)
        maxChrHit(chr,1) = a;
    else
        maxChrHit(chr,1) = NaN;
    end
end

[~, chrSortOrder] = sort( maxChrHit, 1, 'descend');

chrRanking = uint32(tiedrank(-maxChrHit));
xHitRanking = uint32(tiedrank(- obj.(rankField)(hitInds) ));

xRanks = zeros(size(hitInds), 'double');
xRanks(hitInds) = xHitRanking;
fID = fopen(filename, 'wt' );
if fID == -1
    error('printTopHits:wrongFileID','the file cannot be created')
end

fprintf( fID, ['chr;pos;rank w/i chr;chr rank;rank;',...
    'LogLHOdds;LikeliHood;Posterior;',...
    'geneID;geneSO;mut.pos in prot;mut in prot;mut in CDS; TAIR link; eFP browser link;',...
    '[mt] tot reads;[mt] alt reads;[mt] SNP ratio;']);
if obj.flagWT
    fprintf( fID, '[wt] tot reads;[wt] alt reads;[wt] SNP ratio;');
end
if ~isempty(obj.xRnaPrior)
    fprintf( fID, 'Array Prior / mean:%1.4f;', nanmean(obj.xArrayPrior));
    fprintf( fID, 'RNA Prior / mean:%1.4f;', nanmean(obj.xRnaPrior));
end


fprintf(fID, '\n');

for  cc = 1:obj.chrNumber% :-1:1
    chr = chrSortOrder(cc);
    cHitSubs = find(hitInds & obj.chromosome == chr);
    [~, cxHitRanking] = sort( obj.(rankField)(cHitSubs), 1, 'descend');
    
    [~, cxHitRankingSubs] = sort( cxHitRanking, 1, 'ascend');
    cxTopHitSubs = cHitSubs(cxHitRanking);
    for ii = 1: min(numel(cxTopHitSubs), numPerChr)
        jj = cxTopHitSubs(ii);
        currEffect = soDict.lookup(obj.geneSO(jj));
        %                      currEffect = terms{1}{so==obj.geneSO(jj)};
        fprintf(fID, '%u;%u;%u;%u;%u;', obj.chromosome(jj), obj.x(jj), cxHitRanking(cxHitRankingSubs(ii)), chrRanking(chr), xRanks(jj) );
        fprintf(fID, '%g;%g;%g;', obj.xLogOdds(jj), obj.xPsel(jj),  obj.xPosteriorNorm(jj));
        fprintf(fID, '%s;%s;%u;', obj.geneID{jj}, currEffect, obj.mutPosProt(jj));
        fprintf(fID, '%s;%s;', obj.mutProt{jj}, obj.mutCDS{jj});
        fprintf(fID, '=HYPERLINK("http://arabidopsis.org/servlets/TairObject?name=%s&type=locus", "go to TAIR");', obj.geneID{jj});
        fprintf(fID, ['=HYPERLINK("', ...
            'http://www.bar.utoronto.ca/efp/cgi-bin/efpWeb.cgi?dataSource=Developmental_Map&modeInput=Absolute&primaryGene=%s', ...
            '&secondaryGene=None&modeMask_low=None&modeMask_stddev=None','", "go to eFP");'], obj.geneID{jj});
        
               fprintf(fID, '%u;%u;%g;', obj.r(jj), obj.q(jj), obj.f(jj));
        if obj.flagWT
               fprintf(fID, '%u;%u;%g;', obj.rw(jj), obj.qw(jj), obj.fw(jj));
        end        
        if ~isempty(obj.xRnaPrior)
            fprintf( fID, '%g;', obj.xArrayPrior(jj) );
            fprintf( fID, '%g;', obj.xRnaPrior(jj) );
        end
        fprintf(fID, '\n');
    end
end
fclose(fID);
end