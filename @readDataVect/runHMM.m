function [obj, varargout] = runHMM(obj, varargin)
% RUNHMM  --  creates `hmm_cont` objects for each chromosome and runs them
% to calculate likelihoods P[ model_m | x_ ]
% the obtained likelihood values are passed into arrays of the parent object
%
% ## INPUT (optional)
% chromosomes to process

obj.lastWarn = lastwarn;
p = inputParser;
addRequired(p, 'obj', @isobject);
addOptional(p, 'chr',  0,  @(x)(isscalar(x) && x<= obj.chrNumber));
addOptional(p, 'plFlag',  '',  @ischar );
%           addParamValue(p,     'exp10',             false, @isscalar);
parse(p, obj, varargin{:});

if isempty(obj.Alpha)
    obj.Alpha = 1;
    fprintf('setting Alpha to 1\n')
end
%%
if isempty(obj.E) || ~all(size(obj.E) == [ obj.Mtot, obj.Np ])
    obj.calcEmission
else
    fprintf('keeping the emission matrix from a previous run...\n')
end
%% which chromosome(s) to process?
if  p.Results.chr>0
    chrV = p.Results.chr;
else
    chrV = 1:obj.chrNumber;
end

%% initialize empty vectors for output probabilities
obj.cPstat = zeros(obj.chrNumber,1);
obj.cPsel = zeros(obj.chrNumber,1);
obj.cNormConst = zeros(obj.chrNumber,1);

%% set the P[z] models
kPstat = obj.pop.Pstat;
kPflat = ones(1, obj.Np)./obj.Np;
kPsel = zeros(1, obj.Np);
if obj.selType
    kPsel(obj.Np) = 1;
else
    kPsel(1) = 1;
end
%% loop through the chromosomes
for chr = chrV
    ticInit = tic;
    fprintf('processing the chromosome\t%u\t...', chr)
    msgLength = 0;
    
    %% initialize the `hmm_cont` model for this chromosome
    obj.HMM{chr} = hmm_cont(obj.pop, obj.E(obj.cSta(chr):obj.cEnd(chr), :));
    
    for linkageLooseningCurrent = obj.Alpha
        if numel(obj.x(obj.chromosome==chr))<2
            continue
        end
        obj.setTransitionMatrix(chr)
        %% 'flat' model
        [xPflat, xkPflat] = obj.HMM{chr}.getLikelihoodOfAModel( kPflat );
        %% static model
        [xPstat, ~] = obj.HMM{chr}.getLikelihoodOfAModel( kPstat);
        %% selection model
        [xPsel, ~] = obj.HMM{chr}.getLikelihoodOfAModel( kPsel );
        xPsel(isnan(obj.xPsel)) = -Inf; %%%%%%%%%
        %% pass to the outer object fields
        obj.xkPflat(obj.ci{chr},:) = xkPflat;
        obj.xPflat(obj.ci{chr}) = xPflat;
        obj.xPstat(obj.ci{chr}) = xPstat;
        obj.xPsel(obj.ci{chr}) = xPsel;
        %% normalize
        obj.cPstat(chr) = calcMarginal(obj.xPstat(obj.ci{chr}));
        obj.cPsel(chr)  = calcMarginal(obj.xPsel( obj.ci{chr}));
        obj.cNormConst(chr) = - ( obj.pop.Np + calcMarginal(xPflat) ) ;
        obj.xLogOdds(obj.ci{chr}) =  xPsel - xPstat;
        
        if ~isempty(obj.xPrior) &&  numel(obj.xPrior) == numel(obj.xPsel)
            obj.xPosterior(obj.ci{chr},1) = obj.xPsel(chr) + obj.xPrior(obj.ci{chr});
            obj.cPosterior(chr) = calcMarginal(obj.xPosterior(obj.ci{chr}));
        else
            warning('readDataVect:run:noPrior', 'Prior is not defined!\n')
        end
        
        if strcmpi('plot', p.Results.plFlag)
            f_out = plotProbabilitiesOnAChromosome(obj, chr);
            varargout{1} = f_out;
        elseif nargout>1
            varargout{1} = [];
        end
        %% log printing block:
        if strcmp(obj.lastWarn, lastwarn)
            obj.lastWarn = '';
        elseif ~isempty(lastwarn)
            obj.lastWarn = lastwarn;
        end
        
        if isempty(obj.lastWarn) || strcmp(obj.lastWarn, 'Directory already exists.')
            fprintf(repmat('\b',1, msgLength));
        else
            fprintf('\n')
        end
        
        textA = sprintf('Alpha =\t%4.2f\t...', linkageLooseningCurrent);
        fprintf(textA)
        msgLength = numel(textA);
    end
    
    fprintf(repmat('\b',1, msgLength));
    fprintf('\b\b\b took\t%4.2f\t s \t SNPs:\t%u \n',  toc(ticInit), obj.M(chr) )
end % chr
end