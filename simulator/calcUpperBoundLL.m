function [ UB ] = calcUpperBoundLL( N, r, t )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here



r = 20*ones(30,1);
mbs = MapBySeqGeneric(N, 0.25*r, r);
mbs.t = ones(numel(r),1)*0.01;

Eub = emission_ub_r([], r, mbs.pop);

% figure; pcolor(mbs.pop.kvect, r, ub)
figure; plot(mbs.pop.kvect, Eub)

mbs.emissionHandle = @emission_ub_r;

mbs.calcEmission();

mbs.initHMM();
[xPout, xkPout] = mbs.HMM.getLikelihoodOfAModel(mbs.HMM.getHiddenStateModel('fla'));


end

