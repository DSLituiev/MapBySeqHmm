function struout = bsxfields(fun, stru1, stru2, varargin)
% multiplies corresponding fields of two structures
%

p = inputParser;
addRequired(  p,        'fun',         @(x)isa(x,'function_handle'));
addRequired(  p,        'stru1',       @isstruct);
addRequired(  p,        'stru2',       @isstruct);
addOptional(  p,  'inputfields', {[]}, @iscellstr);
% addOptional(  p, 'inputfields', {[]}, @iscell);
% addParamValue(p,         'dim',    1, @isscalar);
parse(p, fun, stru1, stru2, varargin{:});

if any(cellfun(@isempty,p.Results.inputfields))
    fnames1 = fieldnames(stru1);
    fnames2 = fieldnames(stru2);
    commonfields = intersect(fnames1, fnames2);
    % clear fnames1 fnames2
else
    fnames1 = p.Results.inputfields;
    fnames2 = p.Results.inputfields;
    commonfields = p.Results.inputfields;
end

if length(fnames1)<length(fnames2)
    struout = struct(stru1);
else
    struout = struct(stru2);
end
    
for ii = 1:length(commonfields)
    struout.(commonfields{ii}) = fun(stru1.(commonfields{ii}), stru2.(commonfields{ii}));
end
