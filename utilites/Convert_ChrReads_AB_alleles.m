%== Extracts the data from the csv file
ChrNumber = 5;

ftag = 'AB1';

filename = 'C:\Users\Dmitri\Documents\MATLAB\Information_and_Statistics\H_M\raw_tables\AB1.csv';
fid = fopen(filename);
CellAllT = textscan(fid, '%u32 %u32 %s %s %u32 %s','delimiter',';','HeaderLines',1);
% CellAll = textscan(fid, '%u32 %u32 %s %s %u32 %u32','delimiter',';','HeaderLines',1);
fclose(fid);

%=== problem: the annotation contains comma-separated variant alleles!
CommaInds =  cellfun( @(x)~isempty(x) , strfind(CellAllT{end}, ',') );

%% the double entries are droped out
% DoubleQ = arrayfun(@(x)textscan( x{:}, '%u32 %u32','delimiter',','),...
%     CellAllT{end}( CommaInds ),'un',0);
% DoubleQ = cell2mat(vertcat( DoubleQ{:}));
% 
% DoubleAlt = arrayfun(@(x)textscan( x{:}, '%s %s','delimiter',','),...
%     CellAllT{4}( CommaInds ),'un',0);   % CommaInds & mutSize<40 
% 
% DoubleAlt = vertcat( DoubleAlt{:});
% DoubleAlt = cat(2, vertcat( DoubleAlt{:,1} ), vertcat( DoubleAlt{:,2} ) );


[chromosome, position, RefC, AltC, r, q] = CellAllT{:};
clear CellAllT

chromosome = chromosome(~ CommaInds);
position = position(~ CommaInds);
RefC = RefC(~ CommaInds);
AltC = AltC(~ CommaInds);
r = r(~ CommaInds);
q = uint32(str2double(q(~ CommaInds)));

TextValid = cat(2, RefC, AltC);




SNP = cellfun(@(x)size(x,2)==1, TextValid);
SNPInd = all(SNP,2);
ShortIns = SNP(:,1)& ~SNP(:,2) & cellfun(@(x)size(x,2)<=2, TextValid(:,2));
ShortDel = ~ SNP(:,1)& SNP(:,2)& cellfun(@(x)size(x,2)<=2, TextValid(:,1));
sum(ShortIns) + sum(ShortDel) + sum(SNPInd)
Shorts = SNPInd|ShortIns|ShortDel;
% CellAll(ShortDel,3)
Ref = char(zeros( size(TextValid,1),2));
Alt = Ref;
% 
 Ref(Shorts,:) = char(RefC(Shorts));
 Alt(Shorts,:) = char(AltC(Shorts));

%% short insertions
AltIns = char(TextValid(ShortIns,2));
RefIns = char(TextValid(ShortIns,1));

% Ref(ShortIns) = ' ';
if all(AltIns(:,1)  == RefIns)
    Alt(ShortIns) = AltIns(:,2);
else
    warning('delmismatch')
end
%% short dels
AltDel = char(TextValid(ShortDel,2));
RefDel = char(TextValid(ShortDel,1));

% Ref(ShortIns) = ' ';
if all(AltDel  == RefDel(:,1))
    Alt(ShortDel) = RefDel(:,2);
else
    warning('delmismatch')
end

% 
% ATchar = char([ATnum{:}]');
%  ChrNumber = 5;
for chr = 1:ChrNumber
    inds = logical(chromosome == chr);
    ChrReads(chr).r = uint32(r(inds));
    ChrReads(chr).q = uint32(q(inds));
    ChrReads(chr).x = uint32(position(inds));
    ChrReads(chr).Alleles.Ref = Ref(inds,:);
    ChrReads(chr).Alleles.Alt = Alt(inds,:);    
end

save(['C:\Users\Dmitri\Documents\MATLAB\Information_and_Statistics\H_M\datafiles\ChrReads_',...
    ftag], 'ChrReads');
