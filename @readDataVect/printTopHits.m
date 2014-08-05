function printTopHits(obj, filename, numPerChr)
soDict = soDictionary();

hitInds = (obj.xPselNorm < obj.xPosteriorNorm);

for  chr = obj.chrNumber:-1:1
    maxChrHit(chr,1) = max(obj.xPosteriorNorm(hitInds & obj.chromosome == chr));
end

[~, chrSortOrder] = sort( maxChrHit, 1, 'descend');

chrRanking = uint32(tiedrank(-maxChrHit));
xHitRanking = uint32(tiedrank(- obj.xPosteriorNorm(hitInds) ));

xRanks = zeros(size(hitInds), 'double');
xRanks(hitInds) = xHitRanking;
fID = fopen(filename, 'wt' );
if fID == -1
    error('printTopHits:wrongFileID','the file cannot be created')
end

fprintf( fID, ['chr;pos;rank w/i chr;chr rank;rank;',...
    'LogLHOdds;LikeliHood;Posterior;',...
    'geneID;geneSO;mut.pos in prot;mut in prot;mut in CDS; TAIR link;',...
    '[mt] tot reads;[mt] alt reads;[mt] SNP ratio;']);
if obj.flagWT
    fprintf( fID, '[wt] tot reads;[wt] alt reads;[wt] SNP ratio;');
end

fprintf(fID, '\n');

for  cc = 1:obj.chrNumber% :-1:1
    chr = chrSortOrder(cc);
    cHitSubs = find(hitInds & obj.chromosome == chr);
    [~, cxHitRanking] = sort( obj.xPosteriorNorm(cHitSubs), 1, 'descend');
    
    [~, cxHitRankingSubs] = sort( cxHitRanking, 1, 'ascend');
    cxTopHitSubs = cHitSubs(cxHitRanking);
    for ii = 1: min(numel(cxTopHitSubs), numPerChr)
        currEffect = soDict.lookup(obj.geneSO(cxTopHitSubs(ii)));
        %                      currEffect = terms{1}{so==obj.geneSO(cxTopHitSubs(ii))};
        fprintf(fID, '%u;%u;%u;%u;%u;', obj.chromosome(cxTopHitSubs(ii)), obj.x(cxTopHitSubs(ii)), cxHitRanking(cxHitRankingSubs(ii)), chrRanking(chr), xRanks(cxTopHitSubs(ii)) );
        fprintf(fID, '%g;%g;%g;', obj.xLogOdds(cxTopHitSubs(ii)), obj.xPsel(cxTopHitSubs(ii)),  obj.xPosteriorNorm(cxTopHitSubs(ii)));
        fprintf(fID, '%s;%s;%u;', obj.geneID{cxTopHitSubs(ii)}, currEffect, obj.mutPosProt(cxTopHitSubs(ii)));
        fprintf(fID, '%s;%s;', obj.mutProt{cxTopHitSubs(ii)}, obj.mutCDS{cxTopHitSubs(ii)});
        fprintf(fID, '=HYPERLINK("http://arabidopsis.org/servlets/TairObject?name=%s&type=locus", "click here");', obj.geneID{cxTopHitSubs(ii)});
        
               fprintf(fID, '%u;%u;%g;', obj.r(cxTopHitSubs(ii)), obj.q(cxTopHitSubs(ii)), obj.f(cxTopHitSubs(ii)));
        if obj.flagWT
               fprintf(fID, '%u;%u;%g;', obj.rw(cxTopHitSubs(ii)), obj.qw(cxTopHitSubs(ii)), obj.fw(cxTopHitSubs(ii)));
        end        
        fprintf(fID, '\n');
    end
end
fclose(fID);
end