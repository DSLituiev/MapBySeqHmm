function [t, q, r, varargout] = sequencePath(obj, varargin)
% simulates sequencing

%% == sample the path
%== generate the sampling points
% obj.t = sort([obj.T_max * rand(obj.SamplingPoints-3, 1); 0; t0; obj.T_max]);
if nargin>1 && ~isempty(varargin{1}) 
    if ~ischar(varargin{1})
        obj.t = varargin{1};
    elseif ~isempty(obj.t) && ischar(varargin{1}) 
    % nothing
    elseif isempty(obj.t) && ischar(varargin{1}) && strcmpi('grid', varargin{1}) 
        obj.t = linspace(0, obj.T_max, obj.SamplingPoints)';
    else
        obj.t = sort(obj.T_max * rand(obj.SamplingPoints, 1));
    end
else
    obj.t = sort(obj.T_max * rand(obj.SamplingPoints, 1));
end

%== sample
Z = bsxfun(@le, obj.tChr, obj.t');
[~, idx1] = max(flipud(Z),[],1); % find the indices (over-use)

% Zs = xv((size(Z,1) - idx1'+1)).*any(Z,1)' + ceil(obj.Nplants./2)* (~any(Z,1)');
%= sampled points:
obj.k_t = obj.kChr((size(Z,1) - idx1'+1)).*any(Z,1)' + obj.kChr(1)* (~any(Z,1)');
% == study the distribution of the intervals between the sampling points
% plotWaitingTimes(SamplingT)

%% == emission
obj.r = random('Poisson', obj.r_expected, [obj.SamplingPoints, 1] );
obj.q = bsxfun(@(x,y)random('binom',x,y), obj.r, obj.k_t/obj.pop.N/2);

[~, t0ind] = min( abs(obj.t -  obj.xCausativeSNP) );
t0 = obj.t(t0ind);

t = obj.t;
q = obj.q;
r = obj.r;
obj.f = obj.q./obj.r;
if nargout>0
    varargout{1} = obj.k_t;
    varargout{2} = t0;
    varargout{3} = t0ind;    
end