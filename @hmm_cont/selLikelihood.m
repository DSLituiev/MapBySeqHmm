function totLh = selLikelihood(obj, Alhpa)
%STATLIKELIHOOD Summary of this function goes here
%   Detailed explanation goes here

kPsel = zeros(1, obj.Np);
if obj.selType
    kPsel(obj.Np) = 1;
else
    kPsel(1) = 1;
end


if Alhpa > 0
    obj.calcT(Alhpa);
    [xPsel, ~] = obj.getLikelihoodOfAModel(kPsel);
    totLh = calcMarginal(xPsel);
else
    totLh = -Inf;
end
end

