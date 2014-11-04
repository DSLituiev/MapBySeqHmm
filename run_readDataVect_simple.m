function [ AR ] = run_readDataVect_simple( x,q,r, chromosome, N, varargin )
%RUN_READDATAVECT_SIMPLE Summary of this function goes here
%   Detailed explanation goes here

AR = readDataVect(chromosome, x,q,r);
AR.pop = N;
AR.chrMap = 1;

theta = 0;
lambda1 = 1;

AR.emissionHandle = @(q, r, study)emissionMixBetaBinomial(q, r, AR.pop, theta, lambda1);

AR.Alpha =  1;
AR.calcEmission;
AR.run();

if nargin>4 && strcmpi(varargin, 'plot')
    AR.normalizeChromosomes;
    
    AR.plotChromosomes('xPsel', 'yscale', 'lin', 'norm', true, 'figure', 'new');
    
%     AR.plotStemsLP( 'ylim', [-.0002,0])
    
%     AR.plotChromosomes('f', 'yscale', 'lin', 'norm', true, 'figure', 'new');
end
end

