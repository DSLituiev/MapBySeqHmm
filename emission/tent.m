function [ P ] = tent( t, T )
%TENT Summary of this function goes here
% 
% 1/T*(1-)

t_rat = bsxfun(@rdivide, t, T);

P = 0.5 - abs(t_rat - 0.5);

P(abs(t_rat-0.5)>0.5) = 0;

P = bsxfun(@times, P, 2./T);
end

