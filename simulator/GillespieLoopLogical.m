%%% Gillespie's Loop Function for a 1D Poisson elementary process

function [state tvector] =  GillespieLoopLogical(Stmax)
%% generate copy number values
state = false(Stmax,1);
state(1:2:end) =  1;
%% generate random variables
tvector = cumsum(log(1./rand(1, Stmax,'single')));
end