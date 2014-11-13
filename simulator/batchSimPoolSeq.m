classdef batchSimPoolSeq < simPooledSeq
    %BATCHSIMPOOLSEQ Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        numRepeatIter;
        xnLogOdds
        xnPsel
        xnPstat
        xn_k
        n_t
        n_t0_predicted
    end
    
    methods
        function obj = batchSimPoolSeq(varargin)
            obj = obj@simPooledSeq(varargin{:});
        end
        
        function iterate(obj, numRepeatIter, varargin)
            if nargin>2 && ischar(varargin{1})
                mode = varargin{1};
            else
                mode = [];
            end
            msgLength = 0;
            obj.numRepeatIter = numRepeatIter;
            
            for nn = obj.numRepeatIter:-1:1
                [obj.xnLogOdds(:,nn), obj.xnPsel(:,nn), obj.xnPstat(:,nn), obj.xn_k(:,nn), obj.n_t0_predicted(nn), obj.n_t(:,nn)] = ...
                    runSilently(obj, mode);                
                msgLength = updateLog( msgLength, nn );
            end
            fprintf('\n')            
        end
    end
    
end

