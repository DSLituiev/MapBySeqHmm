function totLh = totLikelihood(obj, Alhpa)
%STATLIKELIHOOD Summary of this function goes here
%   Detailed explanation goes here
kPflat = ones(1, obj.Np)./obj.Np;

if Alhpa > 0
    obj.calcT(Alhpa);
    [xPsel, ~] = obj.getLikelihoodOfAModel(kPflat);
    totLh = calcMarginal(xPsel);
else
    totLh = -Inf;
end
end

