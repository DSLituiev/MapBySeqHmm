function [BBL, UBB] = mixBetaBinomUniform(q,r, N, theta)
% mixture components for a heterozygous population
%
% inputs:
% q     - number of mutant reads per locus
% r     - total number of reads per locus
% N     - population size
% theta - overdispersion parameter of beta-binomial distribution 
%         (theta = 0 => pure binomial)
%
% outputs:
% BBL   - beta-binomial pdf (lower/heterozygous-only part)
% UBB   - convolved uniform with beta-binomial distribution
%         sum_{z=0}^{2 N} ( BB(k|n, f_z) * U ) =
%         = 1/(2 N) * sum_{z=0}^{2 N} BB(k|n, f_z)

% note that f_full goes from 0 to 1:

if  (N - floor(N)) > 0
    error('mixBetaBinomUniform:nonIntegerN', 'non-integer population size N!')
end

if ~isreal(theta) || theta<0
    error('mixBetaBinomUniform:complexTheta', 'complex or negative theta!')
end

f_full = 0:(1/(2*N)):1;

% split computation into two parts to overcome memory problems;
% precompute the sum of the upper half of the b-b matrix
betaBinomU = 10.^logBetaBinomialThetaMu0(q, r, f_full(f_full > 0.5), theta);
sumBBU = sum(betaBinomU , 2);
clear betaBinomU

%== beta-binomial distribution
% take only lower half (as the population is heterozygous
BBL = 10.^logBetaBinomialThetaMu0(q, r, f_full(f_full <= 0.5), theta);

%== convolved uniform with beta-binomial distribution
% take both halves of the b-b
UBB = 1/(2*N) * ( sum( BBL , 2) + sumBBU );
