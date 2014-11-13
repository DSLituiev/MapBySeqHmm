classdef readDataVect < dynamicprops % < handle
    properties
        chromosome
        x
        q
        r
        f
        qw
        rw
        fw
        fmMeanF
        fwMeanF        
        fmMedianF
        fwMedianF
        flagWT = false;
        preSelectionGen = 0;
        qual % quality
        notaRepeat;
        contrib;
        c;
        E;
        dx;        
        dxFiltered;
        emissionHandle;
        Alpha;
        ci = {}; % chromosome indices (cell)
        cSta; % chromosome start (vector)
        cEnd; % chromosome end   (vector)
        cNormConst;
        M; % number of SNPs on each chromosome (vector)
        xLL; % log-likelihood of observed data
        xLPost % log-posterior
        chrNumber;
        cMaxX % chromosome length
        Mtot;
        annotation;
        pop;
        chrMap;
        xPflat;
        xPstat;
        xPsel;
        xLogOdds;
        xPselNorm;
        xkPflat;
        cPsel;
        cPstat;
        cPtot;
        Pz;
        prevFig;
        prevAxes;
        prevLine;
        prevYLims;
        prevYNames;
        xTrue = 0;
        cTrue = 0;
        selType = true; % 'true' for positive selection, 'false' for negative.
        xPosterior
        cPosterior
        xPosteriorNorm
        cPostTot
        xPy
        f0
        z
        HMM
    end
    
    properties (SetAccess = private, GetAccess = private)
        lastWarn = ''
    end
    
    properties
        xPriorRaw;
        xPrior;
        geneID;
        geneSO;
        mutCDS;
        mutProt;
        mutPosProt;
        xRnaPresence
        xRnaPrior
        xArrayPrior
        snpFrequencyInEcotypes;        
        snpEcotypesInfo;
        xTemp
    end
    
    properties  (Dependent = true)
        Np; % N + 1 , where N is the number of individuals
    end
    
    methods
        
        function obj = readDataVect(chromosome, x,q,r, varargin)
            if nargin>0 && numel(q) == numel(x) && numel(r) == numel(x)
                obj.chromosome = chromosome;
                obj.x = x;
                obj.q = q;
                obj.r = r;
                obj.f = q./r;
                obj.chrNumber = max(chromosome(:));
                obj.Mtot = numel(x);
                obj.applyFunctionChromosomeWise(@max, 'cMaxX', 'x');
                if (nargin>5)
                    %                     obj.qual = varargin{1};
                    %                     obj.notaRepeat = varargin{2};
                    %                     if (nargin>6)
                    %                         obj.annotation = varargin{3};
                    %                     end
                end
                obj = chrInds(obj);
            end
        end
        
        %% linkage map
        setLinkageMap(obj,  mapPath)
        
        function applyFunctionChromosomeWise(obj, fhandle, fieldOut, fieldIn, varargin)
            assert(isprop(obj, fieldIn))
            assert(isprop(obj, fieldOut))
            obj.(fieldOut) = zeros(obj.chrNumber, 1);
            for cc = obj.chrNumber:-1:1
                obj.(fieldOut)(cc) = feval(fhandle, obj.(fieldIn)(obj.chromosome == cc) );
            end
        end
       %% visualize statistics
        f_out = visualizeStat(obj)
        %% set population object
        function set.pop(obj, N)
            if isnumeric(N)  && isscalar(N)
                obj.pop = population(N);
            elseif isobject(N)
                obj.pop = N;
            else
                warning('readDataVect:pop:unknownType', 'unknown type of input')
            end
            
        end
        
        %% unmixing
        [ obj, mixtObj, mu, iTrue, varargout ] = unmixRepeats( obj, firstParam, varargin );
        
        function [mixtObj] = unmix(obj, varargin)
            obj.calcDxMin;
            if nargin>1 && ~isempty(varargin{1})
                plotFl = varargin{1};
            else
                plotFl = '';
            end
            [obj, mixtObj, mu] = unmixRepeats( obj, 'dx', plotFl,'modeNum', 2);
            
            fprintf('___________________________________\n' )
            fprintf('mixture modes\t log10(dx)\t dx\n' )
            fprintf('lowest (spurious)\t%g\t%u\n', min(mu.dx) , round(10.^min(mu.dx)) )
            fprintf('highest    (true)\t%g\t%u\n', max(mu.dx) , round(10.^max(mu.dx)) )
            
            if any(isnan(obj.contrib))
                warning('unmix:NaNs', 'NaN values in the contribution vector! replacing with "1"')
                obj.contrib(isnan(obj.contrib)) = 1;
            end
        end
        
        function obj = calcDxMin( obj )
            obj = calcDx(obj,  @(x)nanmin(x, [], 2));
        end

        function obj = filterDx(obj, varargin )
            obj.dxFiltered = zeros(size(obj.dx));
            if nargin>2
                chromosomes = varargin{1};
            else
                chromosomes = 1:obj.chrNumber;
            end
            
            for chr = chromosomes                
%                 KERNEL = min(KERNEL, floor(numel(obj.dx(obj.chromosome == chr))/2) - 3 );
%                 fprintf('chromosome %u, kernel size: %u\n', chr, KERNEL)
                obj.dxFiltered(obj.chromosome == chr) = ...
                    smooth(obj.x(obj.chromosome == chr), obj.dx(obj.chromosome == chr), 'loess');
            end          
        end
        
        function obj = calcDx( obj, fh )
            obj.dx = zeros(numel(obj.x) , 1);
            for cc = obj.chrNumber: -1:1
                inds = logical(cc == obj.chromosome);
                dxRaw = [NaN; diff(double(obj.x(inds)))];
                dxRaw( dxRaw < 0) = NaN;
                dxPair = [dxRaw, circshift(dxRaw,-1)];
                obj.dx(inds,1) = feval( fh, dxPair);
            end
        end
       %% filtering
       filterFields(obj, field, fh, varargin)
        %% filter one field
        function filterField(obj, fieldName, inds)
            if ~isempty(obj.(fieldName))
                obj.(fieldName) = obj.(fieldName)(inds);
            end
        end
        %% chromosome retrieval
        function subobj = getChromosome(obj, cc, varargin)
            if nargin>2 && any(strcmpi('old', varargin))
                subobj = chromProbO();
            else
                subobj = chromProb();
            end
            propNames = properties(obj);
            inds = logical(cc == obj.chromosome);
            for ii = 1:numel(propNames)
                if any(strcmpi(propNames{ii}, properties(subobj)))
                    s = size(obj.(propNames{ii}),1);
                    switch s
                        case obj.Mtot
                            subobj.(propNames{ii}) = obj.(propNames{ii})(inds, :);
                        case obj.chrNumber
                            subobj.(propNames{ii}) = obj.(propNames{ii})(cc, :);
                        otherwise
                            subobj.(propNames{ii}) = obj.(propNames{ii});
                    end
                end
            end
            subobj.M = sum(inds);
        end
        %% set Np: number of individuals + 1
        function Np = get.Np(obj)
            Np = size(obj.E, 2);
            if ~(Np == obj.pop.Np)
                obj.calcEmission;
                Np = size(obj.E, 2);
                if ~(Np == obj.pop.Np)
                    dbstop if warning
                    warning('readDataVect:wrongNp', 'the Np is wrong')
                end
            end
        end
        
        %% normalization over the chromosomes
        function normalizeChromosomes(obj)
            obj.cPtot = calcMarginal( [obj.cPstat(:), obj.cPsel(:)], 2 );
            %             obj.cPostTot = calcMarginal( [obj.cPstat, obj.cPosterior], 2 );
            
            obj.cPtot(isnan(obj.cPtot)) = 0;
            for chr = obj.chrNumber:-1:1
                obj.xPselNorm(obj.ci{chr})      = obj.xPsel(obj.ci{chr}) - obj.cPtot(chr);
                %                 obj.xPosteriorNorm(obj.ci{chr}) = obj.xPosterior(obj.ci{chr}) - obj.cPostTot(chr);
            end
            
            normFactor = calcMarginal(obj.xPselNorm);
            if isnan(normFactor); normFactor = 0; end
            obj.xPselNorm  = obj.xPselNorm - normFactor;
            
            if ~isempty(obj.xPrior)
                obj.xPosteriorNorm = obj.xPselNorm + obj.xPrior;
                obj.xPosteriorNorm = obj.xPosteriorNorm - calcMarginal(obj.xPosteriorNorm);
            end
            %             normFactor = calcMarginal(obj.xPosteriorNorm);
            %             obj.xPosteriorNorm  = obj.xPosteriorNorm - normFactor;
            
        end
        
        f = plotChromosome2D(obj, chr);
        
        function varargout = calcMleZ(obj, varargin)
            [~, obj.z] =  max(obj.xkPflat,[], 2 );
            obj.z = obj.z -1;
            if nargin>1 && strcmpi(varargin{1}, 'plot')
                obj.plotChromosomes('z', 'yscale', 'lin', 'figure', 'new',...
                    'plotfun', @(x,y)stairs(x,y, 'linewidth', 2), 'ylim', [0, obj.pop.N] );
            end
            varargout{1} = obj.z;
        end
        
        [obj, varargout] = runHMM(obj, varargin)
        
        function includeRnaPresence(obj)
            obj.xPrior(~obj.xRnaPresence) = -Inf;
        end
        %% find chromosome indices and other chromosome-related constants
        function obj = chrInds(obj)
            obj.ci = cell(obj.chrNumber, 1);
            for chr = 1: obj.chrNumber
                obj.ci{ chr } = (obj.chromosome == chr);
                obj.M(chr) = sum(obj.ci{chr});
                                
                if obj.M(chr)>0; 
                    obj.cSta(chr) = find(obj.ci{ chr }, 1, 'first');
                    obj.cEnd(chr) = find(obj.ci{ chr }, 1, 'last');
                else
                    obj.cSta(chr) = obj.Mtot+1;
                    obj.cEnd(chr) = obj.Mtot+1;
                end                
                obj.cNormConst = zeros(obj.chrNumber, 1);
            end
        end
        %% emission
        obj = calcEmission(obj);
        %% transition
        function setTransitionMatrix(obj, chr, varargin)
            t =  0.01 * mapPhysicalToGeneticPositionCentiMorgan(obj, chr);
            
            if nargin>3 && isscalar(varargin{1})
                linkageLoosening = varargin{1};
            else
                linkageLoosening = 1;
            end
                 
            if nargin>3 && ischar(varargin{2})
                modelName = varargin{2};
            else
                modelName = 'HMM';
            end
            if ~isprop(obj, 'HMM')
                P = addprop(obj, modelName);
            end
            if iscell(obj.(modelName) ) && numel(obj.(modelName) )>=chr && isobject(obj.(modelName){chr}) 
                 obj.(modelName){chr}.t = t;
            else % initialise it
                 obj.(modelName){chr} = hmm_cont(obj.pop, [], t);
            end
            obj.(modelName){chr}.calcT(linkageLoosening);
        end
        
        function x_cM = mapPhysicalToGeneticPositionCentiMorgan(obj, chr)
            % interpolates position of the data points based on the
            % supplied genetic map
            assert( ~isempty(obj.chrMap(chr)) )
            if isscalar(obj.chrMap(chr)) && ~isstruct(obj.chrMap(chr))
                x_cM =  obj.chrMap(chr) * obj.x(obj.ci{chr});
            else
                x_cM = interp1(obj.chrMap(chr).nt, obj.chrMap(chr).cM, double(obj.x(obj.ci{chr})),'pchip','extrap');
            end
        end
                
        %% plotting
        varargout = plotChromosomes(obj, fields, varargin);
        
        varargout = plotOneChromosome(obj, fields, chr, varargin);
        
        plotSnpRatio(obj, KERNEL, varargin)
        
        %% plot likelihood and posterior
        function plotStems(obj, fields, varargin)
            
             if nargin>2 && isscalar(varargin{1}) && isnumeric(varargin{1})
                chr0 = varargin{1};
                 fieldPlotFun = @(x)obj.plotOneChromosome(chr0, x{:});
                 args = varargin(2:end);
            else                
                 fieldPlotFun = @(x)obj.plotChromosomes(x{:});
                 args = varargin;
             end
            
            markerSz = 4;
            fieldPlotFun({obj, fields, ...
                'plotfun', @(x,y)stem(x,y, 'MarkerSize', markerSz, 'linewidth', 0.8),... % 'MarkerEdgeColor', 'r',
                'exp10',false, 'yscale', 'lin', args{:}});
            
            set(obj.prevLine(:, end), 'BaseValue',( obj.prevYLims(1)-2) );
            set(obj.prevAxes(:, end), 'TickDir', 'out')
        end
        
        %% plot stems for Likeilihood and Posterior
        plotStemsLP(obj, varargin)
        
        %% scatter plot mutant vs wild type pool statistics
        f_out = plotScatterMW(obj, varargin);
        
        %% clear all the plots and delete them from the object register
        function clearPlots(obj)
            for ii = 1:numel(obj.prevFig);
                try
                    delete(obj.prevFig(ii));
                catch err
                    if  (strcmp(err.identifier, 'MATLAB:hg:udd_interface:CannotDelete'))
                        warning('MATLAB:readDataVect:deletedFig', 'the fig has been already deleted');
                    end
                end
            end
            obj.prevFig = [];
            obj.prevAxes = [];
            obj.prevLine = [];
            obj.prevYLims = [];
            obj.prevYNames = {};
        end
        
        printTopHits(obj, filename, varargin) ;
        
        function filterF(obj, KERNEL, varargin)
            %% requires `fastmedfilt1d` module !
            obj.fmMedianF = zeros(size(obj.f));
            obj.fmMeanF = zeros(size(obj.f));
            if ~isempty(obj.fw)
                obj.fwMeanF = zeros(size(obj.f));
                obj.fwMedianF = zeros(size(obj.f));
            end
            if nargin>2
                chromosomes = varargin{1};
            else
                chromosomes = 1:obj.chrNumber;
            end
            
            for chr = chromosomes
                if sum(obj.chromosome == chr) > 0
                    KERNEL = min(KERNEL, floor(numel(obj.f(obj.chromosome == chr))/2) );
                    fprintf('chromosome %u, kernel size: %u\n', chr, KERNEL)
                    obj.fmMedianF(obj.chromosome == chr) = fastmedfilt1d(obj.f(obj.chromosome == chr), KERNEL);
                    obj.fmMeanF(obj.chromosome == chr) = smooth(obj.q(obj.chromosome == chr), KERNEL)./smooth(obj.r(obj.chromosome == chr), KERNEL);
                    if ~isempty(obj.fw)
                        obj.fwMedianF = fastmedfilt1d(obj.fw(obj.chromosome == chr), KERNEL);
                        obj.fwMeanF(obj.chromosome == chr) = smooth(obj.qw(obj.chromosome == chr), KERNEL)./smooth(obj.rw(obj.chromosome == chr), KERNEL);
                    end
                end
            end
        end
        %% Baum-Welch : depricated
        obj =  AR.runBaumWelch(obj, chr);
        %%
        obj = set_cMaxX(obj, filePath);
%         function [cc, xIndOnChr] = findChrByIndex(obj, ind)
%             assert( all(ind <= obj.Mtot), 'index/ces out of bound!')
%             cInds = bsxfun(@le, obj.cSta(:)', ind(:)) & bsxfun(@le, ind(:), obj.cEnd(:)');
%             assert( all(sum(cInds,2) == 1), 'ambigous mapping!')
%             dict = double(1:obj.chrNumber)';
%             cc = double(cInds)*dict;
%             xIndOnChr = ind - obj.cSta(cc)' + 1;
%         end
    end
end