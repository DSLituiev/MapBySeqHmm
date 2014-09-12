function [E, c] = emissionMixBetaBinomial(q, r, pop, theta, lambda1)
% emissionMixBetaBinomial returns emission BetaBinomial-Uniform mixture probabilities
%
%

if ~isreal(theta)|| ~isreal(lambda1)
    error('emissionMixBetaBinomial:complexParams', 'complex theta or lambda!')
end
%% P[q|r,k]
[BB, UBB] = mixBetaBinomUniform(q,r, pop.N, theta);

E = bsxfun(@plus, (1-lambda1) .* UBB , bsxfun(@times, lambda1, BB) );
c = 1;
