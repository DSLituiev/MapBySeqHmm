classdef hmm_cont < handle
%HMM_CONT -- a class for handling of continous space/time Markov chains
%   
% Inputs:
% =======
% `pop` --  a `population` object, 
%             describing the hidden state space
% `t`   -- time/space coordinate 
%             (without differencing)
% `E`   -- emission matrix
%

    properties
        t
        Pz
        M
        Np
        E
        T
        A
        linkageLoosening = 1;
        logAlpha
        logBeta
        pop
        xkPplain
        xPplain
        selType = true; % 'true' for positive selection, 'false' for negative.
    end
    
    methods
        %% initialization / construction
        function obj = hmm_cont(population_obj, emission_matrix, dist_genetic)
            obj.pop = population_obj;
            
            if nargin<2
                return
            end
            
            if ~isempty(emission_matrix)          
                obj.E = emission_matrix;                
                obj.M = size(obj.E, 1);
                obj.Np = size(obj.E, 2);
                assert(obj.pop.Np == obj.Np)
            end
            
            if nargin>2
                obj.t = dist_genetic(:);
                if ~isempty(emission_matrix)
                    assert(numel(obj.t) == obj.M, 'dimension mismatch: size(t,1) ~= M ' )
                else
                    obj.M = numel(obj.t);
                    obj.Np = obj.pop.Np;
                end
            else
                obj.t = [];
            end
        end
        
        %% transition
        function obj = calcT(obj, varargin)
            if nargin>1
                obj.linkageLoosening = varargin{1};
            end
            
            assert(~isempty(obj.t), 'the node (recombination) distances "t" have not been set!')
            % assert(~isempty(obj.chrMap), 'chromosome genetic map has not been supplied! \n please assign obj.chrMap = cm;\n cm.nt[M x 1];\n cm.cM[M x 1]' )
            assert( isscalar(obj.linkageLoosening), 'linkageLoosening factor must be a scalar')
            %= calculate Transition matrices
            obj.T = zeros(obj.pop.Np, obj.pop.Np, obj.M-1);
            %             t = Alpha * 0.01* abs(diff(recdistances(obj.chrMap(chr), obj.x(obj.ci{chr})))) ;
            delta_t = diff(obj.t, 1);
            %             obj.t = linkageLoosening * 0.01* obj.recdistances(chr);
            %= calculate
            for ii = 1:(obj.M-1)
                if delta_t(ii)~=0 && ~isinf(delta_t(ii))
                    obj.T(:,:, ii) = expm(obj.pop.Q* delta_t(ii));
                elseif  delta_t(ii)==0
                    obj.T(:,:, ii) = eye(obj.pop.Np);
                elseif delta_t(ii)== Inf
                    obj.T(:,:, ii) = repmat(obj.pop.Pstat, [ obj.pop.Np, 1]);
                end
            end
            %= sanity check
            if   any(obj.T(:)<0)
                if abs(obj.T(obj.T(:)<0)) > 1e-25
                    error('transition3D_expm:TNegInput', 'T matrix has negative values!')
                else
                    [~, msgid] = lastwarn;
                    if ~strcmpi(msgid, 'transition3D_expm:TNegInput')
                        warning('transition3D_expm:TNegInput', ['T matrix has very small negative values (possibly numeric error)! \n',...
                            ' Replacing negative numbers by zeros'])
                    end
                    obj.T(obj.T(:)<0) = 0;
                end
            end
            
        end
        %% crossMatrix A
        function obj = crossMatr(obj)
            assert(~isempty(obj.E), 'provide the emission matrix!')            
%             assert(~isempty(obj.T), 'define transition matrix T first!');
            if isempty(obj.T)
                obj.calcT();
            end
            
            obj.A =  bsxfun(@times, ...
                permute( obj.E(1:end-1, :), [2, 3, 1] ), obj.T ) ;
        end
        %% cumulative matrices 'logAlpha' and 'logBeta'
        function obj = cumMatr(obj)
            if isempty(obj.A)
                obj.crossMatr();
            end
            
            if isempty(obj.Np)
                warning('cumMatrSafe:emptyNp', 'define T first!')
                return
            end
            
            obj.logAlpha = -inf(obj.M, obj.Np);
            obj.logBeta  = -inf(obj.M, obj.Np);
            
            %% forward
            Ac = obj.E( end, :)';
            obj.logAlpha(obj.M, 1: obj.Np) = log10(Ac');
            
            scalePrev = 0;
            scaleA = zeros(obj.M, 1);
            
            for m = (obj.M - 1):-1:1
                Ac = obj.A(:,:, m) * (Ac * 10.^(scaleA(m+1) - scalePrev) );
                obj.logAlpha(m, :)= log10(Ac) - scaleA(m+1);
                scalePrev = scaleA(m + 1);
                scaleA(m) = - max( obj.logAlpha( m, :) );
            end
            
            %% backward
            Bc = ones(1, obj.Np);
            obj.logBeta(1, 1: obj.Np) = 0;
            
            scalePrev = 0;
            scaleB = zeros(obj.M, 1);
            
            for m = 2:1:obj.M
                Bc = Bc * obj.A(:,:, m-1)' * (10.^(scaleB(m-1) - scalePrev));
                obj.logBeta(m, :) = log10(Bc) - scaleB(m-1);
                scalePrev = scaleB(m-1);
                scaleB(m) = - max(obj.logBeta(m, :));
            end
        end
        %% Forward-Backward with no assumption about the model
        function runFBinternal(obj)
            if isempty(obj.logAlpha ) || isempty(obj.logAlpha )
                obj.cumMatr();
            end
            
            obj.xkPplain = obj.logAlpha + obj.logBeta;
            obj.xPplain = calcMarginal(obj.xkPplain, 2);
            if any(isnan(obj.xPplain ))
                warning('runFBinternal:someEntriesAreNaN', 'some entries in the probability vector are NaN!')
            end
            if all(isnan(obj.xPplain))
                error('runFBinternal:allEntriesAreNaN', 'all entries in the probability vector are NaN!')
            end
        end
        %% Apply a provided assumption about the model and return the results
        function [xPout, xkPout] = getLikelihoodOfAModel(obj, model_P_z)
            if isempty(obj.xkPplain)
                obj.runFBinternal;
            end
            xkPout = bsxfun(@plus, obj.xkPplain, log10(model_P_z(:)'));
            xPout = calcMarginal(xkPout, 2);
        end
    end
end

