function [Reads, ChrNumber, annotation] = readSequencingDataCsv(dataPath, varargin)
%% check the input parameters
prsr = inputParser;
addRequired(prsr, 'dataPath', @(x)(ischar(x) & exist(x, 'file')) );
%
addOptional(prsr,'AnnotFlag', '', @ischar);
% addParamValue(p,     'modeNum',             4, @isscalar);
%
parse(prsr, dataPath, varargin{:});
%%
fid = fopen(dataPath);

flagSkipTextChromosomes = true;
% 1 chr	
% 2 pos
% 3 rating
% 4 refAllele
% 5 altAllele
% 6 geneID	geneSO	
% 8 mutCDS	mutProt	
% 10 mutPosCDS	mutPosProt
% 12 distance
% 13 wt_totCount	
% 14 wt_refCount		
% 15 wt_altCount		
% mt_totCount
% mt_refCount
% mt_altCount

header = textscan(fid, repmat('%s ',[1, 20]), 1, 'delimiter', ';');
fclose(fid);

if strcmpi( 'wt_altCount', header)
    flagWT = true;
else
    flagWT = false;
end

indEcotype = strcmpi(header, 'snpEcotypesInfo');
flagEcotype = any(indEcotype);

dataPattern = cell(20,1);
dataPattern(1:20) = {'%s'};
dataPattern(3) = {'%f'}; % quality
dataPattern([2,7,11:15]) = {'%u'};

% dataPattern(19) = {'%f'};
% dataPattern([20]) = {'%u8'};

if flagWT && ~isempty(header{16}{:}) % 17, 18
    dataPattern([2,7,11:18]) = {'%u'};
end

indRepeatType = strcmpi('repeatType', header);
flagRepeats = any(indRepeatType);
if flagRepeats
    dataPattern{indRepeatType} = '%s';
    dataPattern{strcmpi('repeatName', header)} = '%s';
end

fid = fopen(dataPath);
data = textscan(fid, strjoin(dataPattern, ' '), 'delimiter', ';', 'HeaderLines', 1);
fclose(fid);

fprintf('data file has been read\n')

chromosomeStr = data{1};

chrNames = unique(chromosomeStr);
chrNumbers = str2double(chrNames );

literalChrNames = find(isnan(chrNumbers));
    
if flagSkipTextChromosomes
    textChrInds = false(size(chromosomeStr));    
    for jj = 1:numel(literalChrNames)
        textChrInds = textChrInds | strcmpi(chromosomeStr, chrNames(literalChrNames(jj)));
    end
    for jj = 1:size(data, 2)
       data{jj}(textChrInds) = [];
    end
    chromosome = uint8( str2double( data{1} ) );
else
    maxChrNum = nanmax(chrNumbers);
    for jj = 1:numel(literalChrNames)
        maxChrNum = maxChrNum +1;
        chrNumbers(literalChrNames) = maxChrNum;
    end 
    chromosome = zeros(size(chromosomeStr), 'uint8');
    for ii = 1:numel(chrNames)
        chromosome(strcmpi(chromosomeStr, chrNames{ii}), 1) = chrNumbers(ii);
    end
end

position = data{2};
r = data{13};
q = data{15};

ChrNumber = max(chromosome);

Reads = readDataVect(chromosome, double(position), double(q), double(r), [], false);

mutationRating = data{3};
BASE = 2;
Reads.xPrior = BASE.^(mutationRating);
Reads.xPrior = Reads.xPrior./sum(Reads.xPrior);
Reads.xPrior = log10(Reads.xPrior);

if flagRepeats
    Reads.notaRepeat = strcmp(data{indRepeatType}, 'NO');
end

if flagEcotype
    % Reads.xPrior(~logical(data{19})) = -Inf;
    Reads.snpFrequencyInEcotypes = double(data{19});
    Reads.snpEcotypesInfo = logical(data{indEcotype}>0);
end

Reads.geneID = data{6};
Reads.geneSO = data{7};
Reads.mutCDS = data{8};
Reads.mutProt = data{9};
Reads.mutPosProt = data{11};

if flagWT
    Reads.flagWT = true;
    Reads.rw = double(data{16});
    Reads.qw = double(data{18});
end

fprintf('%u reads extracted\n', Reads.Mtot)

% 1...5: chr; pos; rating; refAllele; altAllele; 
% 6...9: wt_totCount;mt_totCount;    wt_refCount; mt_refCount;  
% 10...11: wt_altCount;mt_altCount;
% 12..17: geneID;geneSO; mutCDS;mutProt; mutPosCDS;mutPosProt;
% 18: distance;
% for ii = 6:11
%     annotation.(header{ii}{1}) = data{ii};
% end

%     fid = fopen('SO_terms.csv');
%     % terms = textscan(fid, '%s %s %s %s', 'delimiter', ',', 'HeaderLines', 0);
%     [terms] = textscan(fid, '%s %s %s %s', 'delimiter', ';', 'HeaderLines', 0);
%     fclose(fid);
%     termsID = strcat('SO_',terms{1});
%     % termsID = terms{3};
%     % termsID = char(termsID);
%     % termsID(:,3) = '_';
%     
%     clear terms
annotation = [];
% annotation = initializeSO(termsID, annotation, numel(position));
  