function [xPrior , geneID, maxHitEffect,...
    maxHitCDS, maxHitAA, maxHitCodone ] = constructPriorStr(annotation)
logCS.SO_stop_gained = 12;
logCS.SO_missense_variant = 10;
logCS.SO_splice_region_variant = 4;
logCS.SO_splice_donor_variant = 5;
logCS.SO_splice_acceptor_variant = 5;
logCS.SO_synonymous_variant = 3;
logCS.SO_5_prime_UTR_variant = 2;
logCS.SO_3_prime_UTR_variant = 2;
logCS.SO_non_coding_exon_variant = 1;
logCS.SO_intron_variant= 1;
logCS.SO_upstream_gene_variant = 0;
logCS.SO_downstream_gene_variant = 0;

fieldsConst = fieldnames(logCS);
fieldsAnnot = fieldnames(annotation);
for ii = 1:numel(fieldsAnnot)
    clear annDouble
    temp = ~cellfun(@isempty, annotation.(fieldsAnnot{ii})(:, 1));
    annDouble = -Inf(size(temp));
    annDouble(temp) = 1;
    annBin.(fieldsAnnot{ii}) = annDouble(:);
end

[~, sC,~] = intersect(fieldsConst, fieldsAnnot);

iC = true(size(fieldsConst));
iC(sC) = false;
logCS = rmfield(logCS, fieldsConst(iC));


strout = bsxfields(@mtimes, logCS, annBin);
fieldsNew = fieldnames(strout);

for ii = numel(fieldsNew):-1:1
    fieldsNewShort{ii,1} = fieldsNew{ii}(4:end);
end

Base = 2;

catExp = catfields(strout);
catExp(isnan(catExp)) = 0;
noAnn = all((catExp==0), 2);

mTot = size(annBin.(fieldsConst{1}),1);

geneID = cell(mTot,1);
maxHitCDS = cell(mTot,1);
maxHitAA = cell(mTot,1);
maxHitCodone = cell(mTot,1);
[~, inds] = max(catExp, [], 2);
for ii = mTot:-1:1
    geneID(ii)   = annotation.(fieldsNew{inds(ii)})(ii,1);
    maxHitCDS(ii)    = annotation.(fieldsNew{inds(ii)})(ii,2);
    maxHitAA(ii)     = annotation.(fieldsNew{inds(ii)})(ii,3);
    maxHitCodone(ii) = annotation.(fieldsNew{inds(ii)})(ii,4);
end
maxHitEffect = fieldsNewShort(inds);
maxHitEffect(noAnn) = {[]};
geneID(noAnn) = {[]};
maxHitCDS(noAnn) = {[]};
maxHitAA(noAnn) = {[]};
maxHitCodone(noAnn) = {[]};

Prior = 1 + sum(Base.^catExp, 2);
% PrDenum = sum(catfields(logCS));

% fnames = fieldnames(logCS);

% logCoeff = [logCS.(fnames{:})]

% Prior = 1 + Base.^(ExNsSt*logCoeff);
% Prior = Prior./ (1 + Base.^(logCoeff));
Prior = Prior./sum(Prior);
xPrior = log10(Prior);
% Px0 = 0.01 + 0.99.*ExNsSt(:,1).*(0.1+0.3.*ExNsSt(:,2)+0.6.*ExNsSt(:,3));
