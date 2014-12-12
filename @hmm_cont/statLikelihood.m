function totLh = statLikelihood(obj, Alhpa)
%STATLIKELIHOOD Summary of this function goes here
%   Detailed explanation goes here
if Alhpa > 0
    obj.calcT(Alhpa);
    [xPstat, ~] = obj.getLikelihoodOfAModel(obj.pop.Pstat);
    totLh = calcMarginal(xPstat);
else
    totLh = -Inf;
end
end

