function [logBB, logUBB] = logMixBetaBinomUniform(k,n,theta, N)



f_full = 0:1/(2*N):1;
    
LL = logBetaBinomialThetaMu0(k, n, f_full, theta);

logBB = LL(:, f<1/2);


U = 1/(2*N);
% sum_{z=0}^{2 N} ( B(k|n, f_z) * U ) = 1/(2 N) * sum_{z=0}^{2 N} B(k|n, f_z)

logUBB = log(sum(10.^(LL + log(U)) , 2));