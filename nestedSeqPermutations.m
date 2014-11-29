classdef nestedSeqPermutations < readDataVect
    %NESTEDSEQPERMUTATIONS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mixture
        bin_num = 8;
        bins
        bin_membership
        ind_vect
        new_ind_vect
        bin_size
    end
    
    methods
        function obj = nestedSeqPermutations(init_obj)
            obj = obj@readDataVect(init_obj.chromosome, init_obj.x,init_obj.q, init_obj.r);
            obj.chrMap = init_obj.chrMap;
            obj.emissionHandle = init_obj.emissionHandle;
            obj.cMaxX = init_obj.cMaxX;
            obj.pop = init_obj.pop;
            obj.contrib = init_obj.contrib;
            obj.filterFields('r', @(x)(x<quantile(obj.r, 0.98))); % mutant reads
            obj.filterFields('q', @(x)(x>7)); % mutant reads
            obj.filterFields('f', @(x)(x<.8)); % SNP ratio
            
            obj.calcDxMin;
            obj.mixture = obj.unmix();
        end
        function split(obj, varargin)
            if nargin > 1 && isscalar(varargin{1})
                obj.bin_num = varargin{1};
            else
                obj.bin_num = 8;
            end
            
            %   obj.bin_membership = fix(2*log10(obj.dx));
            if obj.bin_num>1
                [~,inds] = sort(obj.dx);
                minBin = numel(obj.dx)/obj.bin_num;                
                obj.bin_membership = halve(inds, minBin);
                obj.bin_num = max(obj.bin_membership);
            else
                obj.bin_membership = ones(size(obj.dx));
            end
            %   figure; plot(obj.bin_membership); xlim([0, numel(obj.bin_membership)])
            
            
            obj.ind_vect = (1:obj.Mtot)';
            for ii = obj.bin_num:-1:1
                obj.bins{ii} = obj.ind_vect(obj.bin_membership == ii);
            end
        end
        
        function permute(obj, varargin)
            if isempty(obj.bin_membership)
                obj.split(varargin{:});
            end
            obj.new_ind_vect = zeros(size(obj.bin_membership));
            rng('shuffle')
            for ii = obj.bin_num:-1:1
                p = randperm(sum(obj.bin_membership == ii));
                temp = obj.ind_vect(obj.bin_membership == ii);
                obj.new_ind_vect(obj.bin_membership==ii) = temp(p);
            end
            
            key_fields = {'q', 'r'};
            for ii = 1:numel(key_fields)
                obj.(key_fields{ii}) = obj.(key_fields{ii})(obj.new_ind_vect);
            end
        end
        
        function [xnPsel, xnLogOdds] = iterate(obj, n_iter, varargin)            
            obj.split(varargin{:});
            warning off;
            msgLength = 0;
%             f1 = figure('name', 'f');hold on
%             f2 = figure('name', 'p_sel'); hold on
%%
                nn = n_iter;
                obj.runHMM('silent', true);
                % obj.normalizeChromosomes;
                xnPsel(:, nn) = obj.xPsel;                
%                 figure(f1); plot(obj.q(obj.chromosome==1)./obj.r(obj.chromosome==1), 'x');
%                 figure(f2) ; plot(obj.xPsel(obj.chromosome==1))
                xnLogOdds(:, nn) = obj.xLogOdds;
                msgLength = updateLog( msgLength, nn );
                
            for nn = n_iter-1:-1:1
                obj.permute();
                obj.E = [];
                obj.runHMM('silent', true, 'keepTransitionMatrix', true);
                % obj.normalizeChromosomes;
                xnPsel(:, nn) = obj.xPsel;                
%                 figure(f1); plot(obj.q(obj.chromosome==1)./obj.r(obj.chromosome==1), 'x');
%                 figure(f2) ; plot(obj.xPsel(obj.chromosome==1))
                xnLogOdds(:, nn) = obj.xLogOdds;
                msgLength = updateLog( msgLength, nn );
            end
            %%
            warning on;
            fprintf('\n')            
        end
    end
    
end

