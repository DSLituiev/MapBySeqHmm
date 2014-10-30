% DescrVars = unique(AB1_cell(:,7));
%
% IDs = char(AB1_cell(:,1));
% IDs = IDs(:,2:end);
%
% AB1_cell(:,1) =  cellstr(IDs);
%
clc; clear all
addpath('C:\Users\Dmitri\Documents\MATLAB\Information_and_Statistics\H_M');
load('C:\Users\Dmitri\Documents\MATLAB\Information_and_Statistics\H_M\datafiles\AB1_annotation_cell.mat')

clear IDs DescrVars
% [Positions, UniqInd, ~] = unique(AB1_cell(:,1));
% Positions = char(Positions);
% Pos = cell2mat(AB1_cell(UniqInd,2));
% DescrVars = char(DescrVars);

VerbosePositions = char(AB1_cell(:,1));

[ PositionID, iaf, ~ ] = unique(VerbosePositions,'rows','first');
[ ~ , ial , ~ ] = unique(VerbosePositions,'rows','last');
Pos = uint32(cell2mat(AB1_cell(iaf,2)));

if logical(length(VerbosePositions) - length(PositionID) - sum(ial - iaf))
    error('the positions are not sorted!')
end

clear VerbosePositions

Chr = AB1_cell(iaf,3);
Chr(strcmp(Chr,'Mt')) = {10};
Chr(strcmp(Chr,'Pt')) = {20};
Chr = cell2mat(Chr);


Ups     = strcmp(AB1_cell(:,7),'UPSTREAM');
Downs   = strcmp(AB1_cell(:,7),'DOWNSTREAM');
Intergs = strcmp(AB1_cell(:,7),'INTERGENIC');

% inds    = cell(size(Positions,1),1);
Annot   = cell(size(PositionID,1),1);
ATnum   = cell(size(PositionID,1),1);
VerboseIndAnnot = cell(size(PositionID,1),1);

for ii = 1:length(PositionID)
    [Annot{ii}, iaAnnot, ~] = unique(AB1_cell(iaf(ii):ial(ii),7));
    VerboseIndAnnot{ii} = iaf(ii) + iaAnnot - 1;
    ATnum{ii} = AB1_cell(VerboseIndAnnot{ii}, 6);
end



clear iaAnnot iaf ial
for ii = 1:length(PositionID)
    if length(Annot{ii})>1
        if any(Ups(VerboseIndAnnot{ii})) &&( any(Downs(VerboseIndAnnot{ii})) || any(Intergs(VerboseIndAnnot{ii})) )
            %== non-informative = upstream OR downstream OR intergenic
            NonInformInds = Ups(VerboseIndAnnot{ii})| Downs(VerboseIndAnnot{ii})|Intergs(VerboseIndAnnot{ii});
            if sum(NonInformInds) == size(Annot{ii},1)
                if any(Ups(VerboseIndAnnot{ii}))
                    Annot{ii} = {'UPSTREAM'};
                else
                    Annot{ii} = {'DOWNSTREAM'};
                end
                ATnum{ii} = {''};
            else
                Annot{ii} = Annot{ii}(~NonInformInds);
                ATnum{ii} = ATnum{ii}(~NonInformInds);
            end
            clear NonInformInds
        end
        
    end
end

clear VerboseIndAnnot Ups Downs Intergs

% nn = 0;
% att = 0;
% for ii = 1:length(Positions)
%     if length(Annot{ii})>1
%         att = max(att, length(unique(ATnum{ii})));
%         nn = nn+1;
%         Annot{ii}
%         % nn = max(nn,length(Annot{ii}));
%     end
% end

Prior.STOP_GAINED = false(length(PositionID),1);
Prior.FRAMESHIFT_CODING = false(length(PositionID),1);
Prior.NON_SYNONYMOUS_CODING = false(length(PositionID),1);
Prior.ESSENTIAL_SPLICE_SITE = false(length(PositionID),1);
Prior.SYNONYMOUS_CODING = false(length(PositionID),1);
Prior.STOP_LOST = false(length(PositionID),1);
Prior.COMPLEX_INDEL = false(length(PositionID),1);
Prior.P5P_UTR = false(length(PositionID),1);
Prior.P3P_UTR = false(length(PositionID),1);
Prior.INTRONIC = false(length(PositionID),1);
Prior.NON_CODING_GENE = false(length(PositionID),1);



for ii = 1:length(PositionID)
    
    if any(strcmpi(Annot{ii},'STOP_GAINED'))
        Prior.STOP_GAINED(ii) = true;
        ATnum{ii} = ATnum{ii}(strcmpi(Annot{ii},'STOP_GAINED') );
        continue
    end
    
    if any(strfindx(Annot{ii},'ESSENTIAL_SPLICE_SITE')) % !!!!!
        Prior.ESSENTIAL_SPLICE_SITE(ii) = true;
        continue
    end
    
    if any(strcmpi(Annot{ii},'NON_SYNONYMOUS_CODING')) ||...
            any(strcmpi(Annot{ii},'NON_SYNONYMOUS_CODING,SPLICE_SITE'))
        Prior.NON_SYNONYMOUS_CODING(ii) = true;
        ATnum{ii} = ATnum{ii}( strcmpi(Annot{ii},'NON_SYNONYMOUS_CODING') |...
            strcmpi(Annot{ii},'NON_SYNONYMOUS_CODING,SPLICE_SITE') );
        continue
    end
    
    if any(strfindx(Annot{ii},'3PRIME_UTR'))
        Prior.P3P_UTR(ii) = true;
        ATnum{ii} =   ATnum{ii}( cell2mat( regexp(Annot{ii}, '.*3PRIME_UTR.*') ) );
        continue
    end
    
    if any(strfindx(Annot{ii},'5PRIME_UTR'))
        Prior.P5P_UTR(ii) = true;
        ATnum{ii} =   ATnum{ii}( cell2mat(regexp(Annot{ii}, '.*5PRIME_UTR.*')) );
        continue
    end
    
    if any(strfindx(Annot{ii},'INTRONIC'))
        Prior.INTRONIC(ii) = true;   % cellfun(@(x)x(1),
        ATnum{ii} =  ATnum{ii}( cell2mat(regexp(Annot{ii}, '.*INTRONIC.*')) );
        continue
    end
    
    if any(strcmpi(Annot{ii},'STOP_LOST'))
        Prior.STOP_LOST(ii) = true;
        ATnum{ii} = ATnum{ii}( strcmpi(Annot{ii},'STOP_LOST') );
        continue
    end
    
    if any(strcmpi(Annot{ii},'SYNONYMOUS_CODING')) ||...
            any(strcmpi(Annot{ii},'SYNONYMOUS_CODING,SPLICE_SITE'))
        Prior.SYNONYMOUS_CODING(ii) = true;
        ATnum{ii} =  ATnum{ii}( strcmpi(Annot{ii},'SYNONYMOUS_CODING') |...
            strcmpi(Annot{ii},'SYNONYMOUS_CODING,SPLICE_SITE') );
        continue
    end
    
    if any(strfindx(Annot{ii},'FRAMESHIFT_CODING'))
        Prior.FRAMESHIFT_CODING(ii) = true;
        continue
    end
    
    if any(strfindx(Annot{ii},'NON_CODING_GENE'))
        Prior.NON_CODING_GENE(ii) = true;
        continue
    end
    
    if any(strfindx(Annot{ii},'COMPLEX_INDEL'))
        Prior.COMPLEX_INDEL(ii) = true;
        continue
    end
    
end

% ind1 = find(strcmpi(cellstr(PositionID),'4#13260324'))
% Annot{ind1}
% Prior.NON_SYNONYMOUS_CODING(ind1)
% ATnum{ind1}

clear ii

save('C:\Users\Dmitri\Documents\MATLAB\Information_and_Statistics\H_M\datafiles\AB1_annotation_ShortCell.mat')
