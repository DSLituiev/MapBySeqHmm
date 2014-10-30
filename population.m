classdef population < hgsetget
    properties
        parent = [];
        N
        Np
        kvect
        f
        Pstat
        Pflat
        Q
    end
    
    methods
        function obj = population(Nin, varargin)
            if nargin>0
                obj.N = Nin;
                obj.Np = Nin + 1;
                obj.kvect = 0:1:Nin;
                if (nargin>1) && any(strcmpi(varargin{1}, {'homozygous','selfed','self', 's'}))
                    obj.f = obj.kvect./obj.N;
                else % heterozygous population from a back-cross
                    obj.f = obj.kvect./(2*obj.N);
                end
                obj = obj.stationaryDistribution;
                obj = obj.calcQmatrix;                
                obj.Pflat = ones(1, obj.Np)./ obj.Np;
            end
            
        end
        
        function obj = changePopulation(obj, Nin, varargin)
            if nargin>1
                obj.N = Nin;
                obj.Np = Nin + 1;
                obj.kvect = 0:1:Nin;
                if (nargin>2) && any(strcmpi(varargin{1}, {'homozygous','selfed','self', 's'}))
                    obj.f = obj.kvect./obj.N;
                else % heterozygous population from a back-cross
                    obj.f = obj.kvect./(2*obj.N);
                end
                obj = obj.stationaryDistribution;
                obj = obj.calcQmatrix;
                obj.Pflat = ones(1, obj.Np)./ obj.Np;
            end
            
        end
        
        function obj = stationaryDistribution(obj)
            obj.Pstat = zeros(1, obj.Np);
            
            if obj.N <= 20;
                ii = 0;
                obj.Pstat(ii+1) = 2^(-obj.N).* nchoosek(obj.N, ii);
                for ii = 1:obj.N
                    obj.Pstat(ii+1) = 2^(-obj.N).*nchoosek(obj.N,ii);
                end
            else % for large numbers
                ii = 0:1:obj.N;
                obj.Pstat(ii+1) = 2^(-obj.N).* round(exp(gammaln(obj.N + 1)-gammaln(ii + 1)-gammaln(obj.N - ii + 1)));
            end
            
        end
        
        function obj = calcQmatrix(obj)
            % infinitisemal transition matrix
            % Calculates infinitesimal transition matrix Q
            % for continuous Ehrenfest process
            %
            % ===== Input: =====
            % N - number of plants/molecules
            
            % short k-vector
            kvectsh = 0:(obj.N-1);
            % initialize
            obj.Q = zeros(obj.Np);
            % increment
            ind = sub2ind([obj.Np obj.Np], 1:obj.N, 2:(obj.N+1));
            obj.Q(ind) = obj.N - kvectsh;
            % decrement
            ind = sub2ind([obj.Np obj.Np], 2:(obj.N+1), 1:obj.N);
            obj.Q(ind) = kvectsh+1;
            % staying same
            ind = sub2ind([obj.Np obj.Np], 1:obj.N+1, 1:obj.N+1);
            obj.Q(ind) = - sum(obj.Q,2);
            % Q is done!
        end
        %         function Qout = get.Q(obj)
        %             if ~isempty(obj.Q)
        %                 Qmatrix(obj);
        %             end
        %             Qout = pbj.Q;
        %         end
    end
end