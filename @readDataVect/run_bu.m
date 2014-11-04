function [obj, varargout] = run(obj, varargin)
% RUN

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

if  p.Results.chr>0
    chrV = p.Results.chr;
else
    chrV = 1:obj.chrNumber;
end

for chr = chrV
    ticInit = tic;
    fprintf('processing the chromosome\t%u\t...', chr)
    msgLength = 0;
    for Alpha0 = obj.Alpha
        if numel(obj.x(obj.chromosome==chr))<2
            continue
        end
        obj = obj.calcT(chr, Alpha0);
        obj = obj.crossMatr(chr);
        obj = obj.cumMatr(chr);
        
        obj = obj.runFBflat(chr);
        obj = obj.runFBstat(chr);
        obj = obj.runFBselection(chr);
        
        obj.xPsel(isnan(obj.xPsel)) = -Inf; %%%%%%%%%
        obj.cPstat(chr) = calcMarginal(obj.xPstat(obj.ci{chr}));
        obj.cPsel(chr)  = calcMarginal(obj.xPsel( obj.ci{chr}));
        
        obj.cNormConst(chr) = - ( obj.pop.Np + calcMarginal(obj.xPflat(obj.ci{chr}) ) ) ;
        obj.xLogOdds(obj.ci{chr}) =  obj.xPsel(obj.ci{chr}) -  obj.xPstat(obj.ci{chr});
        
        if ~isempty(obj.xPrior) &&  numel(obj.xPrior) == numel(obj.xPsel)
            obj.xPosterior(obj.ci{chr},1) = obj.xPsel(chr) + obj.xPrior(obj.ci{chr});
            obj.cPosterior(chr) = calcMarginal(obj.xPosterior(obj.ci{chr}));
        else
            warning('readDataVect:run:noPrior', 'Prior is not defined!\n')
        end
        
        if strcmpi('plot', p.Results.plFlag)
            varargout{1} = figure;
            subplot(2,1,1)
            plot(obj.x(obj.ci{chr}), obj.xPstat(obj.ci{chr}), 'g-')
            hold all
            plot(obj.x(obj.ci{chr}), obj.xPout(obj.ci{chr}), 'r-')
            plot(obj.x(obj.ci{chr}), obj.xPflat(obj.ci{chr}), 'b-')
            legend({'stationary', 'selection', 'flat'})
            
            subplot(2,1,2)
            plot(obj.x(obj.ci{chr}), obj.xPstat(obj.ci{chr})  + obj.cNormConst(chr) , 'g-')
            hold all
            plot(obj.x(obj.ci{chr}), obj.xPsel(obj.ci{chr}) + obj.cNormConst(chr) , 'r-')
            legend({'stationary norm',  'selection norm'})
            
        elseif nargout>1
            varargout{1} = [];
        end
        obj.T{chr} = [];
        obj.A{chr} = [];
        obj.logAlpha{chr} = [];
        obj.logBeta{chr} = [];
        
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
        
        textA = sprintf('Alpha =\t%4.2f\t...', Alpha0);
        fprintf(textA)
        msgLength = numel(textA);
    end
    
    
    fprintf(repmat('\b',1, msgLength));
    fprintf('\b\b\b took\t%4.2f\t s \t SNPs:\t%u \n',  toc(ticInit), obj.M(chr) )
end % chr
end