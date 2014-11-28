classdef MapBySeqGeneric < handle
    %MAPBYSEQGENERIC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        t
        r
        q
        f
        emissionHandle
        E
        preSelectionGen = false;
        contrib = [];
        xPsel
        xPstat
        xLogOdds
        pop
        flagWT = false;
        HMM
    end
    
    methods
        function obj = MapBySeqGeneric(varargin)
            if ~isempty(varargin)
                if nargin>=1
                    N = varargin{1};
                    if ~isobject(N)
                        obj.pop = population(N);
                    else
                        obj.pop = N;
                    end
                end
                argList = {'q', 'r', 't', 'x'};
                
                for ii = 1:(numel(argList) )
                    if numel(varargin)>=ii+1
                        obj.(argList{ii}) = varargin{ii+1};
                        
                    end
                end
            end
        end
        %         function initPop(obj, Nplants)
        %             obj.pop = population(Nplants);
        %         end
        
        calcEmission(obj);
        
        function initHMM(obj, varargin)
            assert(~isempty(obj.t), 'supply the genetic distances first!')
            obj.calcEmission();
            
            if nargin>1 && isnumeric(varargin{1})
                chromosome = varargin{1};
                for ii = numel(chromosome)
                    obj.HMM{chromosome(ii)} = hmm_cont(obj.pop, obj.E, obj.t);
                end
            else
                obj.HMM = hmm_cont(obj.pop, obj.E, obj.t);
            end
            
        end
    end
    
end

