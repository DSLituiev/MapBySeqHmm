% Error using hmm_cont
% Too many input arguments.
close all
clear all
clc
%%
USERFNCT_PATH = './dependencies';
addpath(genpath(USERFNCT_PATH));
addpath(genpath('.'));
%%
N = 5;
q = [1, 2, 2, 1, 3, 4, 1, 3]';
r = [4, 6, 3, 8, 6, 7, 5, 9]';
f = q./r;
chromosome = ones(1, numel(r));
% x = ceil(sort(rand(numel(r),1))* 7* numel(r));
x = [ 4     5     8    17    30    32    44    53]';
xPselNormGroundTruth = - 0.01*[89.4467   90.8446   90.4565   95.1174   91.6904   87.1670   85.8210   92.6367]';
%%
pop = population(N);

theta = 0; lambda1 = 1;
emissionHandle = @(q, r)emissionMixBetaBinomial(q, r, pop, theta, lambda1);
emission_matrix = feval(emissionHandle, q,r);

HMM = hmm_cont(pop, emission_matrix, 0.01*x);
HMM.calcT();
HMM.runFBstat();
HMM.xPstat
HMM.runFBflat();

assert( all(all(bsxfun(@minus, mtimesx(HMM.pop.Pstat, HMM.T),  HMM.pop.Pstat) < 1e-12 ) ), 'stationarity for Kolmogorov fwd equation is violated')
HMM.runFBselection()

xPselNorm = HMM.xPsel - calcMarginal(HMM.xPsel);
%%
[ AR ] = run_readDataVect_simple( x, q, r, chromosome, N, 'plot' );

xPselNorm = AR.xPsel - calcMarginal(AR.xPsel);
%%
figure
plot(x, xPselNorm )

assert( norm(xPselNorm -xPselNormGroundTruth ) <= 7.8e-07)
