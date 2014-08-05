%== matches the Chromosome structur
clear all
ftag = 'AB1';

load(['C:\Users\Dmitri\Documents\MATLAB\Information_and_Statistics\H_M\datafiles\',...
    'AB1_annotation_ShortCell']);

load(['C:\Users\Dmitri\Documents\MATLAB\Information_and_Statistics\H_M\datafiles\',...
    'ChrReads_AB1']);

ChrReads = rmfield(ChrReads, {'AT','Pr'});
PrFNames = fieldnames(Prior);

ATrep = cellfun(@(x)size(x,1), ATnum);

fprintf('%u non-unique ATs\n', sum(logical(ATrep )) )

ATInds = false(size(ATrep));
ATInds(ATrep ==1) = true;
%== collect repeated entries
ATRepeat = strcmp( cellfun(@(x)x{1,:} ,   ATnum(ATrep==2), 'un', 0), cellfun(@(x)x{2,:} ,   ATnum(ATrep == 2), 'un', 0));
ATInds(ATrep==2) =  ATRepeat;
ATRepeat = strcmp( cellfun(@(x)x{1,:} ,   ATnum(ATrep==3), 'un', 0), cellfun(@(x)x{2,:} ,   ATnum(ATrep == 3), 'un', 0)) &...
    strcmp( cellfun(@(x)x{1,:} ,   ATnum(ATrep==3), 'un', 0), cellfun(@(x)x{3,:} ,   ATnum(ATrep == 3), 'un', 0));
ATInds(ATrep==3) =  ATRepeat;
ATnum(~ATInds) = {{' '}};
ATchar = char([cellfun(@(x)x{1,:} ,   ATnum, 'un', 0)]);

for ii = 1:5    
    % yy = Pos(Chr == ii);  
    xshift = int32(cellfun(@(x)size(x,2), cellstr(ChrReads(ii).Alleles.Alt)) - cellfun(@(x)size(x,2), cellstr(ChrReads(ii).Alleles.Ref)));
    x = int32(ChrReads(ii).x) - xshift;
    y = int32(Pos(Chr == ii));
     [~, xInds, yInds] = intersect(x,y);
    %==   x    ~  y(Inds)
    
    
%     [Val, Inds] = min( abs(bsxfun(@minus, x, y') ) , [], 2);
%     %==   x    ~  y(Inds)
%     %==   Val  =  | x - y(Inds) |
%     
%     [IndsSort si]= sort(Inds);
%     z = find(~logical(diff(IndsSort)));
%     %== z -- indices of the repeated Inds
%     %[ IndsSort(z), IndsSort(z+1) ]
%       RepDist = subsrefx(Val(si), {[z , z+1]},'()' );
%      subsrefx(y(IndsSort(z)), {xor(RepDist(:,1), RepDist(:,2))},'()' )
%      subsrefx(y(IndsSort(z)), {~xor(RepDist(:,1), RepDist(:,2))},'()' )
%        
%     fprintf( 'chr.# %u\tmaximal distance mismatch %u \n', ii, max(Val))
%     fprintf('\t\tredundant entries:\t%u\n', numel(unique(ChrReads(ii).x)) - numel(unique(Inds)) )
%     
    for ff = 1:length(PrFNames)
        CurrPr.(PrFNames{ff}) =  Prior.(PrFNames{ff})(Chr == ii);
        ChrReads(ii).Pr.(PrFNames{ff}) = false(size(ChrReads(ii).x));
        ChrReads(ii).Pr.(PrFNames{ff})(xInds) = CurrPr.(PrFNames{ff})(yInds);
    end
    CurrAT =  ATchar(Chr == ii,:);
    ChrReads(ii).AT(xInds,:) = CurrAT(yInds,:);
end

ind0 = find(strcmpi(cellstr(PositionID),'4#13260324'))
Prior.NON_SYNONYMOUS_CODING(ind0)
ATnum{ind0}
ATchar(ind0,:)
find(strcmp(cellstr(ATchar), 'AT4G26180') )

 ind1 = find(ChrReads(4).x == 13260324 )
ChrReads(4).Pr.NON_SYNONYMOUS_CODING(ind1)
ChrReads(4).AT(ind1,:)
% find(strcmp(cellstr(ChrReads(4).AT),'AT4G20550') )
% repinds = find( strcmp( [cellfun(@(x)x{1}, ATnum, 'un',0)] ,'AT4G20550') )
% Pos(repinds)

ii =4;  CurrAT =  ATchar(Chr == ii,:); CurrAT(1162,:)

find(strcmp(cellstr(CurrAT), 'AT4G26180') )
 
% find(strcmp(cellstr(ChrReads(4).AT), 'AT4G26180') )

clear ii ff PrFNames Val Inds XX YY Chr Annot AB1_cell ATnum Positions Prior Pos

clear  ans chr ChrNumber  PositionsID
 
save(['C:\Users\Dmitri\Documents\MATLAB\Information_and_Statistics\H_M\datafiles\ChrReads_',...
    ftag], 'ChrReads');
