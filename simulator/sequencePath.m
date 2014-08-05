function [t,q,r, varargout] = sequencePath(fChr, tChr, par)
% simulates sequencing

%% == sample the path
%== generate the sampling points
% t = sort([par.T_max * rand(par.SamplingPoints-3, 1); 0; t0; par.T_max]);
t = sort(par.T_max * rand(par.SamplingPoints, 1));
%== sample
Z = bsxfun(@le, tChr, t');
[~, idx1] = max(flipud(Z),[],1); % find the indices (over-use)

% Zs = xv((size(Z,1) - idx1'+1)).*any(Z,1)' + ceil(par.Nplants./2)* (~any(Z,1)');
%= sampled points:
fAtSamplingPoints = fChr((size(Z,1) - idx1'+1)).*any(Z,1)' + fChr(1)* (~any(Z,1)');
% == study the distribution of the intervals between the sampling points
% plotWaitingTimes(SamplingT)

%% == emission
r = random('Poisson', par.r_expected, [par.SamplingPoints, 1] );
q = bsxfun(@(x,y)random('binom',x,y), r, fAtSamplingPoints/par.Nplants/2);

[~, t0ind] = min( abs(t -  par.cSNPpos) );
t0 = t(t0ind);

if nargout>0
    varargout{1} = fAtSamplingPoints;
    varargout{2} = t0;
    varargout{3} = t0ind;    
end