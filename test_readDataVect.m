close all
clear all
clc
%%
USERFNCT_PATH = './dependencies';
addpath(genpath(USERFNCT_PATH));
addpath(genpath('.'));
%%
N = 20;
q = [1, 2, 2, 1, 3, 4, 1, 3]';
r = [4, 6, 3, 8, 6, 7, 5, 9]';
f = q./r;
chromosome = ones(1, numel(r));
% x = ceil(sort(rand(numel(r),1))* 7* numel(r));
x = [ 4     5     8    17    30    32    44    53]';
xPselNormGroundTruth = - 0.01*[89.4467   90.8446   90.4565   95.1174   91.6904   87.1670   85.8210   92.6367]';
%%

[ AR ] = run_readDataVect_simple( x,q,r,chromosome, N, 'plot' );

xPselNorm = AR.xPsel - calcMarginal(AR.xPsel);
figure
plot(x, xPselNorm )

assert(norm(xPselNorm -xPselNormGroundTruth ) <= 7.8e-07)
