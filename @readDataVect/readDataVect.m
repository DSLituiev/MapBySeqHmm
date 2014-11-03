classdef readDataVect < handle
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
        T = {};
        A = {};
        logAlpha = {};
        logBeta = {};
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
        xkPout;
        xPout;
        xPstat;
        xkPstat;
        cPstat;
        xPsel;
        xkPsel;
        cPsel;
        cPtot;
        xPflat;
        xkPflat;
        xLogOdds;
        xPselNorm;
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
        nPlLines;
    end
    
    methods
        
        function nPlLines = get.nPlLines(obj)
            nPlLines = numel(obj.prevYNames);
        end
        
        function obj = readDataVect(chromosome, x,q,r, varargin)
            if nargin>0 && numel(q) == numel(x) && numel(r) == numel(x)
                obj.chromosome = chromosome;
                obj.x = x;
                obj.q = q;
                obj.r = r;
                obj.f = q./r;
                obj.chrNumber = max(chromosome);
                obj.Mtot = numel(x);
                obj.applyFunctionChromosomeWise(@max, 'cMaxX', 'x');
                obj.xPstat = -Inf(obj.Mtot , 1);
                obj.xPsel = -Inf(obj.Mtot , 1);
                obj.xPflat = -Inf(obj.Mtot , 1);
                obj.xPout = -Inf(obj.Mtot , 1);
                obj.xPosterior = zeros(obj.Mtot , 1);
                obj.xPosteriorNorm = zeros(obj.Mtot, 1);
                
                obj.xPselNorm = -Inf(obj.Mtot , 1);
                obj.cPstat = -Inf(obj.chrNumber , 1);
                obj.cPsel = -Inf(obj.chrNumber , 1);
                obj.cPosterior = -Inf(obj.chrNumber , 1);
                
                obj.xLogOdds =  zeros(obj.Mtot , 1);
                
                
                obj.T = cell(obj.chrNumber, 1);
                obj.A = cell(obj.chrNumber, 1);
                obj.logAlpha = cell(obj.chrNumber, 1);
                obj.logBeta = cell(obj.chrNumber, 1);
                if (nargin>5)
                    %                     obj.qual = varargin{1};
                    %                     obj.notaRepeat = varargin{2};
                    %                     if (nargin>6)
                    %                         obj.annotation = varargin{3};
                    %                     end
                end
                % uARnique(obj.chromosome)
                obj = chrInds(obj);
            end
        end
        
        function applyFunctionChromosomeWise(obj, fhandle, fieldOut, fieldIn, varargin)
            assert(isprop(obj, fieldIn))
            assert(isprop(obj, fieldOut))
            obj.(fieldOut) = zeros(obj.chrNumber, 1);
            for cc = obj.chrNumber:-1:1
                obj.(fieldOut)(cc) = feval(fhandle, obj.(fieldIn)(obj.chromosome == cc) );
            end
        end
        %% visualize statistics
        function visStat(obj)
            nBins = max( floor(numel(obj.dx)/25), 7 );
            %==
            df = 0.01*ceil(100/nBins);
            fx = df/2:df:(1-df/2);
            %==
            dr = max(1/nBins, 0.01);
            rx = df/2:dr:(log10(max(obj.r))-dr/2);
            %==
            figure
            
            subplot(2,2,1)
            myhist100(fx,  obj.f, 'r', 'edgecolor', 'none');
            set(gca, 'tickDir', 'out')
            title('SNP ratio [f]')
            xlabel('f');
            
            subplot(2,2,2)
            [~, modeX] = myhist100LogX(rx,  log10(obj.r), 'r', 'edgecolor', 'none');
            set(gca, 'tickDir', 'out')
            title( sprintf('coverage / read number per locus [r], \n median[r] = %u, mode[r] = %u', floor(median(obj.r)), floor(modeX)) );
            xlabel('r [log_{10}-scaled]');
            
            
            subplot(2,2,3)
            scatter(log10(obj.r), obj.f, 'r');
            xlabel('log_{10}(r)'); ylabel('f')
            set(gca, 'tickDir', 'out')
            title('read number [r] vs SNP ratio [f]')
            set(gca, 'ylim', [0,1])
            
            subplot(2,2,4)
            
            myhistLogX(nBins , log10(obj.dx), 'r')
            %             [hy, hx] = hist( log10(obj.dx), floor(numel(obj.dx)/25));
            %             bar(10.^hx, hy)
            set(gca, 'tickDir', 'out', 'xscale', 'log')
            %             xlim([0, 700])
            title('nt to closest neighbour SNP [${\Delta} x$]', 'interpreter', 'latex')
            xlabel('${\Delta} x$ [log$_{10}$-scaled]', 'interpreter', 'latex');
            
            %             hist(obj.qual, 0:10:700);
            %             xlim([0, 700])
            %             title('quality')
        end
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
        function obj = unmix(obj, varargin)
            obj.calcDxMin;
            if nargin>1 && ~isempty(varargin{1})
                plotFl = varargin{1};
            else
                plotFl = '';
            end
            [obj, mu] = unmixRepeats( obj, 'dx', plotFl,'modeNum', 2);
            
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
        function obj = filter(obj, field, fh, varargin)
            inds = fh(obj.(field));
            if ~isprop(obj, field)
                error('readDataVect:filter:noFieldFound', 'no field  "%s" found!', field)
            end
            if isempty(obj.(field) )
                error('readDataVect:filter:noFieldFound', 'the field  "%s" is empty! \n Cannot perform filtering', field)
            end
            if sum(inds) == 0
                error('readDataVect:filter:empty', 'no entries are left after filtering')
            end
            if nargin>3
                field2 = varargin{1};
                fh2 = varargin{2};
                inds = inds &  fh2(obj.(field2));
            end
            
            obj.chromosome = obj.chromosome(inds);
            obj.x = obj.x(inds);
            obj.q = obj.q(inds);
            obj.r = obj.r(inds);
            obj.f = obj.q./obj.r;
            if obj.flagWT
                obj.qw = obj.qw(inds);
                obj.rw = obj.rw(inds);
                obj.fw = obj.qw./obj.rw;
            end
            
            obj.chrNumber = max(obj.chromosome);
            obj.Mtot = numel(obj.x);
            obj.xPout = -Inf(obj.Mtot, 1);
            obj.xPstat = -Inf(obj.Mtot, 1);
            obj.xPsel  = -Inf(obj.Mtot, 1);
            obj.xPselNorm = -Inf(obj.Mtot, 1);
            obj.xPflat = -Inf(obj.Mtot, 1);
            obj.xLogOdds = zeros(obj.Mtot, 1);
            obj.E = [];
            obj.contrib = [];
            obj.calcDxMin();
            
            
            if ~isempty(obj.xPrior)
                obj.xPrior = obj.xPrior(inds);
                obj.xPosterior = -Inf(obj.Mtot, 1);
                obj.xPosteriorNorm = -Inf(obj.Mtot, 1);
            end
            
            optFields = {'qual', 'dxFiltered', 'notaRepeat', 'geneID', 'geneSO', 'mutCDS',...
                'mutProt', 'mutPosProt', 'xRnaPrior', 'xRnaPresence', 'xArrayPrior', 'snpFrequencyInEcotypes', 'snpEcotypesInfo'};
            
            for ii = numel(optFields):-1:1
                filterField(obj, optFields{ii} , inds);
            end
            
            obj.chrInds();
            
        end
        %%
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
        %% Np: number of individuals + 1
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
        
        %% operations on chromosomes
        function normalizeChromosomes(obj)
            obj.cPtot    = calcMarginal( [obj.cPstat, obj.cPsel], 2 );
            %             obj.cPostTot = calcMarginal( [obj.cPstat, obj.cPosterior], 2 );
            
            obj.cPtot(isnan(obj.cPtot)) = 0;
            for chr = obj.chrNumber:-1:1
                obj.xPselNorm(obj.ci{chr})      = obj.xPsel(obj.ci{chr})      - obj.cPtot(chr);
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
        
        function plotChromosome2D(obj, chr)
            Z =  bsxfun( @minus, obj.xkPflat , calcMarginal( obj.xkPflat , 2) );
            figure; 
            surf( obj.x(obj.chromosome == chr)*1e-6, obj.pop.kvect, Z(obj.chromosome == chr,:)', 'linestyle', 'none'); 
            view(0,90)
        end
        
        function varargout = calcMleZ(obj, varargin)
            [~, obj.z] =  max(obj.xkPflat,[], 2 );
            obj.z = obj.z -1;
            if nargin>1 && strcmpi(varargin{1}, 'plot')
                obj.plotChromosomes('z', 'yscale', 'lin', 'figure', 'new',...
                    'plotfun', @(x,y)stairs(x,y, 'linewidth', 2), 'ylim', [0, obj.pop.N] );
            end
            varargout{1} = obj.z;
        end
        
        function [obj, varargout] = run(obj, varargin)
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
                % [0, obj.cEnd(1:end-1)] + 1; 
                % cumsum(obj.M);                
                obj.cNormConst = zeros(obj.chrNumber, 1);
            end
        end
        %% emission
        obj = calcEmission(obj);
        %% transition
        function dx_cM = recdistances(obj, chr)
            % interpolates distances between the data points based on the
            % supplied genetic map
            x_cM = interp1(obj.chrMap(chr).nt, obj.chrMap(chr).cM, double(obj.x(obj.ci{chr})),'pchip','extrap');
            dx_cM = abs(diff(x_cM));
        end
        
        function obj = calcT(obj, chr, Alpha)
            if isempty(obj.chrMap)
                error('calcT:noChrMap', 'chromosome genetic map has not been supplied! \n please assign obj.chrMap = cm;\n cm.nt[M x 1];\n cm.cM[M x 1]' )
            end
            if ~isscalar(Alpha)
                error('calcT:vectorAlpha', 'Alpha must be a scalar')
            end
            %= calculate Transition matrices
            obj.T{chr} = zeros(obj.pop.Np, obj.pop.Np, obj.M(chr)-1);
            %             tvect = Alpha * 0.01* abs(diff(recdistances(obj.chrMap(chr), obj.x(obj.ci{chr})))) ;
            tvect = Alpha * 0.01* obj.recdistances(chr);
            %= calculate
            for ii = 1:(obj.M(chr)-1)
                if tvect(ii)~=0 && ~isinf(tvect(ii))
                    obj.T{chr}(:,:, ii) = expm(obj.pop.Q*tvect(ii));
                elseif  tvect(ii)==0
                    obj.T{chr}(:,:, ii) = eye(obj.pop.Np);
                elseif tvect(ii)== Inf
                    obj.T{chr}(:,:, ii) = repmat(obj.pop.Pstat, [ obj.pop.Np, 1]);
                end
            end
            %= sanity check
            if   any(obj.T{chr}(:)<0)
                if abs(obj.T{chr}(obj.T{chr}(:)<0)) > 1e-25
                    error('transition3D_expm:TNegInput', 'T matrix has negative values!')
                else
                    [~, msgid] = lastwarn;
                    if ~strcmpi(msgid, 'transition3D_expm:TNegInput')
                        warning('transition3D_expm:TNegInput', ['T matrix has very small negative values (possibly numeric error)! \n',...
                            ' Replacing negative numbers by zeros'])
                    end
                    obj.T{chr}(obj.T{chr}(:)<0) = 0;
                end
            end
            
        end
        %% cumulative matrices 'logAlpha' and 'logBeta'
        function obj = cumMatr(obj, chr)
            % function obj = cumMatrSafe(obj)
            if ~isempty(obj.Np)
                
                obj.logAlpha{chr} = -inf(obj.M(chr), obj.Np);
                obj.logBeta{chr}  = -inf(obj.M(chr), obj.Np);
                
                %% forward
                Ac = obj.E( obj.cEnd(chr), :)';
                obj.logAlpha{chr}(obj.M(chr), 1: obj.Np) = log10(Ac');
                
                scalePrev = 0;
                scaleA = zeros(obj.M(chr), 1);
                
                for m = (obj.M(chr) - 1):-1:1
                    Ac = obj.A{chr}(:,:, m) * (Ac * 10.^(scaleA(m+1) - scalePrev) );
                    obj.logAlpha{chr}(m, :)= log10(Ac) - scaleA(m+1);
                    scalePrev = scaleA(m + 1);
                    scaleA(m) = - max( obj.logAlpha{chr}( m, :) );
                end
                
                %% backward
                Bc = ones(1, obj.Np);
                obj.logBeta{chr}(1, 1: obj.Np) = 0;
                
                scalePrev = 0;
                scaleB = zeros(obj.M(chr)-1, 1);
                
                for m = 2:1:obj.M(chr)
                    Bc = Bc * obj.A{chr}(:,:, m-1)' * (10.^(scaleB(m-1) - scalePrev));
                    obj.logBeta{chr}(m, :) = log10(Bc) - scaleB(m-1);
                    scalePrev = scaleB(m-1);
                    scaleB(m) = - max(obj.logBeta{chr}(m, :));
                end
                
            else
                warning('cumMatrSafe:emptyNp', 'define T first!')
            end
        end
        %% crossMatrix A
        function obj = crossMatr(obj, chr)
            if isempty(obj.E)
                obj.calcEmission
            end
            
            if ~( numel(obj.T) < chr) && ~isempty(obj.T{chr})
                obj.A{chr} =  bsxfun(@times, ...
                    permute( obj.E(obj.cSta(chr):(obj.cEnd(chr)-1), :), [2, 3, 1] ), obj.T{chr} ) ;
            else
                error('crossMatr:emptyT', 'define transition matrix T first!');
            end
            
            %             fprintf('matrix A calculated for the chr %u\n', chr);
        end
        %% forward-backward wrapping functions
        function obj = runFB(obj, Pin, chr)
            if nargin<2
                obj.Pz = ones(1, obj.Np);
            else
                obj.Pz = Pin;
            end
            obj = runFBinternal(obj, chr);
        end
        
        function obj = runFBstat(obj, chr)
            obj.Pz = obj.pop.Pstat;
            obj = runFBinternal(obj, chr);
            obj.xkPstat(obj.ci{chr}, :)  = obj.xkPout(obj.ci{chr}, :) ;
            obj.xPstat(obj.ci{chr}) = obj.xPout(obj.ci{chr});% median(obj.xPout);
        end
        
        function obj = runFBflat(obj, chr)
            obj.Pz = ones(1, obj.Np)./obj.Np;
            obj = runFBinternal(obj, chr);
            obj.xkPflat(obj.ci{chr}, :)  = obj.xkPout(obj.ci{chr}, :) ;
            obj.xPflat(obj.ci{chr}) = obj.xPout(obj.ci{chr});
        end
        
        function PFF = getPoutFullFlat(obj, chr)
            if isempty(obj.xkPflat)
                obj = runFBflat(obj, chr);
            end
            PFF = obj.xkPout;
        end
        
        function obj = runFBselection(obj, chr)
            obj.Pz = zeros(1, obj.Np);
            if obj.selType
                obj.Pz(obj.Np) = 1;
            else
                obj.Pz(1) = 1;
            end

%             lambda = 1;
%             obj.Pz = lambda.^(obj.pop.N - obj.pop.kvect)./factorial(obj.pop.N - obj.pop.kvect) .* exp(-lambda);
            
%             obj.Pz = binopdf(obj.pop.kvect, 2*obj.pop.N, 1/2*(1 - 1/obj.pop.N) );
%             obj.Pz = obj.Pz./sum(obj.Pz);
%             figure; plot( obj.Pz )
            
            obj = runFBinternal(obj, chr);
            obj.xkPsel(obj.ci{chr}, :)  = obj.xkPout(obj.ci{chr}, :) ;
            obj.xPsel(obj.ci{chr}) = obj.xPout(obj.ci{chr});  % median(obj.xPout);
        end
        
        %% FB - final step
        function obj = runFBinternal(obj, chr)
            if isempty(obj.A)|| isempty(obj.logAlpha ) || isempty(obj.logAlpha )
                obj = obj.crossMatr(chr);
                obj = obj.cumMatr(chr);
            end
            if size(obj.Pz,2) ~= obj.pop.Np && ~isscalar(obj.Pz)
                error('runFBinternal:wrongPzSize','wrong prior distribution for z_m is submitted!')
            end
            
            obj.xkPout(obj.ci{chr}, :) = bsxfun(@plus, (obj.logAlpha{chr} + obj.logBeta{chr}), log10(obj.Pz(:)'));
            % obj.xkPout(obj.ci{chr}, :) = bsxfun(@plus, nansum([obj.logAlpha{chr}, obj.logBeta{chr}],2), log10(obj.Pz));
%             infInds = isinf(obj.logAlpha{chr}) & isinf(obj.logBeta{chr});
%             figure; p= pcolor(obj.logBeta{chr}); set(p, 'linestyle','none')
%             obj.xkPout(obj.Sta(chr), :)
            obj.xPout(obj.ci{chr}) = calcMarginal(obj.xkPout(obj.ci{chr}, :), 2);
            if any(isnan(obj.xPout(obj.ci{chr})))
                warning('runFBinternal:someEntriesAreNaN', 'some entries in the probability vector are NaN!')
            end
            if all(isnan(obj.xPout(obj.ci{chr})))
                error('runFBinternal:allEntriesAreNaN', 'all entries in the probability vector are NaN!')
            end
        end
        
        %% plotting
        varargout = plotChromosomes(obj, fields, varargin);
        
        varargout = plotOneChromosome(obj, fields, chr, varargin);
        
        plotSnpRatio(obj, KERNEL, varargin)
        
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
    
        function plotScatterMW(obj, varargin)
            if ~obj.flagWT
                return
            end
            if nargin>1 && isscalar(varargin{1}) && isnumeric(varargin{1})
                chr0 = varargin{1};
                inds = obj.chromosome == chr0 ;
            else
                inds = true(size(obj.x));
            end
            figure('name', 'mt vs wt SNP ratio');
            da(1) = scatter(obj.f(inds), obj.fw(inds), 'b.');
            
            hold on

            da(2) =  scatter(obj.fmMedianF(inds ), obj.fwMedianF(inds ), 4, [0.7,.3, .3], 'o');
            da(3) =  scatter(obj.fmMeanF(inds ), obj.fwMeanF(inds ), 4, [0.2,.8, .2], 'o');
            daInfo = {'raw [all]', 'median', 'mean'};
            
            dl(1) = plot([0,1], [0,1], 'm--');
            dl(2) = plot([0,.5], [0.5, 0], 'r-');
            dlInfo = {'m = w (errors)', 'm + w = 1/2 (expected)'};
            
            deInfo = {};
            if ~isempty(obj.snpEcotypesInfo)
                de(1) =  scatter(obj.f(inds & obj.snpEcotypesInfo), obj.fw(inds & obj.snpEcotypesInfo), 'g.');
                deInfo = 'BG ecotype';                
                uistack(da(2:3), 'top')
                uistack(dl, 'top')
            end
            
            axis equal
            xlabel('f_{mt}'); ylabel('f_{wt}')
            xlim([0,1]); ylim([0,1])
            legend([da, dl, de], {daInfo{:}, dlInfo{:}, deInfo}, 'location', 'northwest')
        end
        
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
            % addpath(fullfile(USERFNCT_PATH, 'fastmedfilt1d'));
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
        %% Baum-Welch
        obj =  AR.runBaumWelch(obj, chr);
        %%
        obj = set_cMaxX(obj, filePath);
    end
end