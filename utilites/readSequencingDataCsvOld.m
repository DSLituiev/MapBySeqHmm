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
data = textscan(fid, ['%u %u %s %s %n',...
                      ' %u %u %u %u %u',...
                      ' %s %s %s'], 'delimiter', ';', 'HeaderLines', 1);
fclose(fid);

fprintf('data file has been read\n')

chromosome = data{1};
position = data{2};
qual = data{5};
p = data{7};
q = data{8};
r = p + q;

% pRefQual = data{9};
% qAltQual = data{10};

ChrNumber = max(chromosome);

Reads = readDataVect(chromosome, double(position), double(q), double(r));
% qual, strcmp(data{11}, 'NO')
    fprintf('%u reads extracted\n', Reads.Mtot)

    
    
% Reads.pRefQual = pRefQual;
% Reads.qAltQual = qAltQual ;

% Reads.repeat = data{11};
% Reads.repeatName = data{12};
% Reads.notaRepeat = strcmp(data{11}, 'NO');

if isempty(regexp( prsr.Results.AnnotFlag,'no.*','once'))
    fprintf('reading the annotation...\n')
    rawAnnotation = data{13};
    clear data;
    
    fid = fopen('SO_terms.csv');
    % terms = textscan(fid, '%s %s %s %s', 'delimiter', ',', 'HeaderLines', 0);
    [terms] = textscan(fid, '%s %s %s %s', 'delimiter', ';', 'HeaderLines', 0);
    fclose(fid);
    termsID = strcat('SO_',terms{1});
    % termsID = terms{3};
    % termsID = char(termsID);
    % termsID(:,3) = '_';
    
    clear terms
            
    p1 = '(?<gene>\w*)';
    p2 = '(?<effect>\w*)';
    p3 = '(?<CDS>\d*)';
    p4 = '(?<AA>\S*)';
    p5 = '(?<Codone>\S*)';
    expr = [p1 ':' p2 ':' p3 ':' p4 ':' p5 ];

    annotation = [];
    annotation = initializeSO(termsID, annotation, numel(position));
    
    for ii = numel(position):-1:1
        annCell{ii} = regexp( rawAnnotation{ii}, '\|' , 'split' );
        for jj = numel(annCell{ii}):-1:1
            inStr = regexp(annCell{ii}{jj}, expr, 'names');
            if ~isempty(inStr)
                annotation = fillInSO(termsID, inStr, annotation, ii);
                %                 Reads.ann(ii).eff(jj) = eff0;
                %                 effects(xx) = eff0(:).effect;
                %                 xx = xx+1;
            end
        end
    end
    
    annotation = cleanSO(termsID, annotation );    
end