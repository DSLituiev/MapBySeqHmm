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
addParamValue(p, 'silent',  false,  @islogical );
addParamValue(p, 'keepTransitionMatrix',  false,  @islogical );
%           addParamValue(p,     'exp10',             false, @isscalar);
parse(p, obj, varargin{:});

if isempty(obj.Alpha)
    obj.Alpha = 1;
    fprintf('setting Alpha to 1\n')
end
%%
if isempty(obj.E) || ~all(size(obj.E) == [ obj.Mtot, obj.Np ])
    obj.calcEmission;
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

log10N_LinkLoosening = log10(numel(obj.Alpha));
%% loop through the chromosomes
for chr = chrV
    ticInit = tic;
    if obj.cSta(chr)<obj.cEnd(chr) && size(obj.E,1)>=obj.cEnd(chr)
        if ~p.Results.silent
            fprintf('processing the chromosome\t%u\t...', chr)
        end
    else
        fprintf('skipping the chromosome\t%u\t... \n', chr)
        continue
    end
    msgLength = 0;
    
    %% initialize the `hmm_cont` model for this chromosome
    
    if ~p.Results.keepTransitionMatrix || isempty(obj.HMM) || numel(obj.HMM)< chr || isempty(obj.HMM{chr}.T)
        obj.HMM{chr} = hmm_cont(obj.pop, obj.E(obj.cSta(chr):obj.cEnd(chr), :));
    else
        obj.HMM{chr}.E = obj.E(obj.cSta(chr):obj.cEnd(chr), :);
        obj.HMM{chr}.resetFlag = true;
    end
    
    xPflat = zeros( obj.M(chr), numel(obj.Alpha) );
    xPstat = zeros( obj.M(chr), numel(obj.Alpha) );
    xPsel  = zeros( obj.M(chr), numel(obj.Alpha) );
    xkPflat  = zeros( obj.M(chr), obj.Np, numel(obj.Alpha) );
    
    for ii = 1:numel(obj.Alpha)
        if numel(obj.x(obj.chromosome==chr))<2
            continue
        end
        if ~p.Results.keepTransitionMatrix || isempty(obj.HMM{chr}.T)
            obj.setTransitionMatrix(chr, obj.Alpha(ii))
        end
        %% 'flat' model
        [xPflat(:, ii) , xkPflat(:, :, ii)] = obj.HMM{chr}.getLikelihoodOfAModel( kPflat );
        %% static model
        [xPstat(:, ii), ~] = obj.HMM{chr}.getLikelihoodOfAModel( kPstat);
        %% selection model
        [xPsel(:, ii), ~] = obj.HMM{chr}.getLikelihoodOfAModel( kPsel );
        xPsel(isnan(xPsel(:, ii)), ii) = -Inf; %%%%%%%%%
        %% plot
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
        
        if ~p.Results.silent
            if isempty(obj.lastWarn) || ...
                    strcmp(obj.lastWarn, 'Directory already exists.') || ...
                    strcmp(obj.lastWarn(1:8), 'T matrix')
                fprintf(repmat('\b',1, msgLength));
            else
                fprintf('\n')
            end
        end
        if ~p.Results.silent 
            textA = sprintf('Alpha =\t%4.2f\t...', obj.Alpha(ii));
            fprintf(textA)
            msgLength = numel(textA);
        end
    end
        %% pass to the outer object fields
        obj.xkPflat(obj.ci{chr},:) = calcMarginal(xkPflat, 3) - log10N_LinkLoosening;
        obj.xPflat(obj.ci{chr},1) =  calcMarginal(xPflat, 2) - log10N_LinkLoosening;
        obj.xPstat(obj.ci{chr},1) =  calcMarginal(xPstat, 2) - log10N_LinkLoosening;
        obj.xPsel(obj.ci{chr},1)  =  calcMarginal(xPsel,  2) - log10N_LinkLoosening;
        %% normalize
        obj.cPstat(chr) = calcMarginal(obj.xPstat(obj.ci{chr}));
        obj.cPsel(chr)  = calcMarginal(obj.xPsel( obj.ci{chr}));
        obj.cNormConst(chr) = - ( obj.pop.Np + calcMarginal(obj.xPflat(obj.ci{chr})) ) ;
        obj.xLogOdds(obj.ci{chr}) =  obj.xPsel(obj.ci{chr}) - obj.xPstat(obj.ci{chr});
        
        if ~isempty(obj.xPrior) &&  numel(obj.xPrior) == numel(obj.xPsel)
            obj.xPosterior(obj.ci{chr},1) = obj.xPsel(obj.ci{chr}) + obj.xPrior(obj.ci{chr});
            obj.cPosterior(chr) = calcMarginal(obj.xPosterior(obj.ci{chr}));
        elseif ~p.Results.silent
            warning('readDataVect:run:noPrior', 'Prior is not defined!\n')
        end
        
     if ~p.Results.silent 
        fprintf(repmat('\b',1, msgLength));
        fprintf('\b\b\b took\t%4.2f\t s \t SNPs:\t%u \n',  toc(ticInit), obj.M(chr) )
     end
end % chr
end