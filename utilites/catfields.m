function cf = catfields(stru, varargin)
% catenates fields of a structure
% if the second input is not specified, operates on all fields
%
% ====== Input  ====== 
%
% - stru                   -- structure
% - inputfields (optional) -- field names to catenate
%
% ====== Output ======
%
% - cf -- an matrix of catenated fields
%

p = inputParser;

addRequired(  p,        'stru',       @isstruct);
addOptional(  p, 'inputfields', {[]}, @iscell);
% addParamValue(p,         'dim',    1, @isscalar);
parse(p, stru, varargin{:});

if ~isempty(p.Results.inputfields)
    fnames = fieldnames(stru);
else
    fnames = p.Results.inputfields;
end

cellc = cell(1,length(fnames));
for ii = 1:length(fnames)
    cellc{ii} = stru.(fnames{ii})(:,1);
end

cf = cell2mat(cellc);