%==== A simulator script for the mapping process

close all; clear all; clc; %<= magic spells

%% add useful functions
addpath('D:\MATLABuserfunctions\mtimesx');
addpath('D:\MATLABuserfunctions\MinMaxSelection');
addpath('D:\MATLABuserfunctions');
%= add the parent path
addpath(cd(cd('..')))

%% == Initialize parameters:

NumRepeatIter = 1000; %= number of iterations
sweepSamplingPoints = (25:5:250)';
NumSweepIter = numel(sweepSamplingPoints);
par.Nplants = 50;
% par.SNPperPlant = 10;
par.r_expected = 12;
par.FalsePositives = 1;
par.cSNPpos = .5;
par.Selection = true;
par.Rec_max = 10; %<= maximal number of recombinations per one path
par.T_max = 1; %<= maximal length of the path

for jj = NumSweepIter:-1:1
    tjj = tic;
    par.SamplingPoints = sweepSamplingPoints(jj);
    
    %% initialize
    estWT = NaN(NumRepeatIter, NumSweepIter);
    
    P_Gt = NaN(NumRepeatIter, NumSweepIter);
    P_MaxP = NaN(NumRepeatIter, NumSweepIter);
    
    t_Gt = NaN(NumRepeatIter, NumSweepIter);
    t_MaxP = NaN(NumRepeatIter, NumSweepIter);
    
    tind_Gt = NaN(NumRepeatIter, NumSweepIter);
    tind_MaxP = NaN(NumRepeatIter, NumSweepIter);
    
    for ii = 1:NumRepeatIter
        %% == generate a path
        %= tChr : vector of mutation locations
        %= fChr : vector of the number of mutant plants
        %= t_Gt : location of the mutation under selection
       %  [tChr, fChr, t_Gt(ii,jj), tind_Gt(ii,jj), ~] = generatePath(par);
         [tChr, fChr, ~, ~, ~] = generatePath(par);
         
        %% simulate sequencing
        [t, q, r, ~, t_Gt(ii,jj), tind_Gt(ii,jj)] = sequencePath(fChr, tChr, par);
        
        %% == estimate mean recombination waiting time
        estWT(ii,jj) =  max(tChr)/size(tChr,1)*par.Nplants;
        
        %% Forward-backward
        % ti1 = tic;
        logPobsDy = runFB(t, q, r,  par.Nplants);
        % fprintf('iter.#\t%3u\tF-B algorithm took\t%3.2g\ts\n', ii, toc(ti1))
        
        %% = find the ML point
        [P_MaxP(ii,jj), tind_MaxP(ii,jj)] = max(logPobsDy);
        t_MaxP(ii,jj) = t(tind_MaxP(ii,jj));
        P_Gt(ii,jj) = logPobsDy( tind_Gt(ii,jj) );
    end
    
    fprintf('sweep.#\t%3u\t took\t%3.2g\ts\n', jj, toc(tjj))    
    
end

save('batchHTSMsimulation', 'par', 'sweepSamplingPoints', 'NumSweepIter', 'NumRepeatIter', 'P_Gt' , 'P_MaxP', 't_Gt', 't_MaxP', 'tind_Gt', 'tind_MaxP')

% load('batchHTSMsimulation')
%% == plot
tHist = -.6:.1:.6;
for ii = 1:NumSweepIter   
   f_tHist(:,ii) = hist(t_MaxP(:,ii) - t_Gt(:,ii), tHist)';   
end

figure('name','time')
surf(f_tHist)

figure('name','time index')
hist(tind_MaxP(:) - tind_Gt(:), floor(NumRepeatIter/12) )

figure('name', 'logL diff')
hist(P_MaxP(:) - P_Gt(:), floor(NumRepeatIter/12) )

figure('name', 'logL Gt')
hist(P_MaxP , floor(NumRepeatIter/12) )

figure('name', 'max logL')
hist(P_MaxP , floor(NumRepeatIter/12) )

figure('name', 'estimated waiting time')
hist(estWT , floor(NumRepeatIter/12) )

%
% xlim([0, par.T_max])
% ylim([ floor(min(logPobsDy)), ceil(max(logPobsDy))] )

