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
    end
    
    methods
%         function obj = MapBySeqGeneric(N)
%             obj.pop = N;
%         end
%         function initPop(obj, Nplants)
%             obj.pop = population(Nplants);
%         end
        calcEmission(obj);
    end
    
end

