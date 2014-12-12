function [ Alhpa_opt, Lh ] = optimizeTransition( obj )
%optimizeTransition Summary of this function goes here
%   Detailed explanation goes here
obj.resetFlag = true;

Alpha0 = 16;

fun = @(x)statLikelihood(obj, x);

options = optimset('TolFun', 1e-3, 'Display', 'iter');

[Alhpa_opt, Lh] = fminsearch(fun, Alpha0);

%% static model

end

