function [ M_noise, ind_noise ] = addSequencingNoise(obj, varargin )
%ADDSEQUENCINGNOISE Summary of this function goes here
%   Detailed explanation goes here

assert(obj.SamplingPoints == numel(obj.q))

if isscalar(varargin{1})
    perc_noise = varargin{1};
    M_noise = obj.SamplingPoints*perc_noise;
    ind_noise = randperm(obj.SamplingPoints, M_noise);
elseif isnumeric(varargin{1}) && numel(varargin{1}) < obj.SamplingPoints
    M_noise =  numel(varargin{1});
    ind_noise = varargin{1};
end



fk_noise = rand(M_noise,1);

obj.q(ind_noise) = bsxfun(@(x,y)random('binom',x,y), obj.r(ind_noise), fk_noise);


end

