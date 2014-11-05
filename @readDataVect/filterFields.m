function filterFields(obj, field, fh, varargin)
%FILTER Summary of this function goes here
%   Detailed explanation goes here
%% filtering

inds = fh(obj.(field));

assert(isprop(obj, field), sprintf('no field  "%s" found!\n', field))

if isempty(obj.(field) )
    error('readDataVect:filter:noFieldFound', 'the field  "%s" is empty! \n Cannot perform filtering', field)
end

assert(sum(inds) ~= 0,  'no entries are left after filtering')

if nargin>3
    field2 = varargin{1};
    fh2 = varargin{2};
    inds = inds &  fh2(obj.(field2));
end

obj.chromosome = obj.chromosome(inds);
obj.x = obj.x(inds);
obj.q = obj.q(inds);
obj.r = obj.r(inds);
obj.f = obj.q./obj.r;
if obj.flagWT
    obj.qw = obj.qw(inds);
    obj.rw = obj.rw(inds);
    obj.fw = obj.qw./obj.rw;
end

obj.chrNumber = max(obj.chromosome);
obj.Mtot = numel(obj.x);
obj.xPstat = -Inf(obj.Mtot, 1);
obj.xPsel  = -Inf(obj.Mtot, 1);
obj.xPselNorm = -Inf(obj.Mtot, 1);
obj.xPflat = -Inf(obj.Mtot, 1);
obj.xLogOdds = zeros(obj.Mtot, 1);
obj.E = [];
obj.contrib = [];
obj.calcDxMin();


if ~isempty(obj.xPrior)
    obj.xPrior = obj.xPrior(inds);
    obj.xPosterior = -Inf(obj.Mtot, 1);
    obj.xPosteriorNorm = -Inf(obj.Mtot, 1);
end

optFields = {'qual', 'dxFiltered', 'notaRepeat', 'geneID', 'geneSO', 'mutCDS',...
    'mutProt', 'mutPosProt', 'xRnaPrior', 'xRnaPresence', 'xArrayPrior', 'snpFrequencyInEcotypes', 'snpEcotypesInfo'};

for ii = numel(optFields):-1:1
    filterField(obj, optFields{ii} , inds);
end

obj.chrInds();

end
