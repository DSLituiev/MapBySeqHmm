classdef simPooledSeq < MapBySeqGeneric
    %SIMPOOLEDSEQ Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        r_expected
        k_t
        FalsePositives = 0;
        SamplingPoints
        xCausativeSNP
        Selection
        T_max
        Rec_max
        tChr
        kChr
        indCausativeSNP
        t0
        est_waiting_time
        HMM
        maxLogOdds
        t0_predicted 
    end
    
    methods
        function obj = simPooledSeq(xCausativeSNP, Nplants, r_expected, SamplingPoints, varargin)
            
            obj.pop = population( Nplants);
            obj.r_expected = r_expected; % expected number of reads
            obj.FalsePositives = 0;
            obj.SamplingPoints = SamplingPoints;
            obj.xCausativeSNP = xCausativeSNP;
            obj.Selection = true;
            if nargin>4 && isscalar(varargin{1})
                obj.T_max  = varargin{1};
            else
                obj.T_max = 1; %<= maximal length of the path
            end
            
            obj.Rec_max = obj.T_max * 10; %<= maximal number of recombinations per one path
        end
        
        [tChr, fChr, t0, dummy , T_max] = generatePath(obj);
        [t,q,r, fAtSamplingPoints] = sequencePath(obj, varargin);
        
        
        [ f ] = plotSeqSim( obj );
        
        function [xLogOdds, xkPsel] = runHMM(obj)
            obj.emissionHandle = @(qq, rr, ff)emissionk0(qq, rr, obj.pop);
            obj.calcEmission();
            obj.HMM = hmm_cont(obj.pop.N, obj.E, obj.t);
            model_P_z_sel = zeros(obj.pop.Np, 1);
            model_P_z_sel(obj.pop.Np - obj.FalsePositives) = 1;
            [obj.xPsel, xkPsel] = getLikelihoodOfAModel(obj.HMM, model_P_z_sel);
            obj.xPstat = getLikelihoodOfAModel(obj.HMM, obj.pop.Pstat);
            obj.xLogOdds = obj.xPsel - obj.xPstat;
            xLogOdds = obj.xLogOdds;
            
            [obj.maxLogOdds, maxind] = max(obj.xLogOdds);
            obj.t0_predicted = obj.t(maxind);
        end
        
        function [xLogOdds, varargout] = runSilently(obj, varargin)
            %% == generate a path
            obj.generatePath();
            %% simulate sequencing
            obj.sequencePath(varargin{:});
            %% == estimate mean recombination waiting time
            obj.est_waiting_time =  max(obj.tChr)/size(obj.tChr,1)*obj.pop. N;
            %% Forward-backward
            xLogOdds = obj.runHMM();
            
            numRecPoints = numel(obj.tChr);
            varargout = {obj.xPsel, obj.xPstat, obj.k_t, obj.t0_predicted, obj.t, numRecPoints };
        end
    end
    
end

