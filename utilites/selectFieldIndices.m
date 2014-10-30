function [ x ] = selectFieldIndices( x , indices)
%selects sub-elements of a structure vector elements
%   for  x.A [ N x 1 ] , x.B [ N x 1 ] etc and  indices [ N x 1 ]  
%   returns y so that: y.A = x.A(indices), y.B = x.B(indices) etc.

fieldNames = fieldnames(x);

for ii = 1:numel(fieldNames )
    if ~isscalar(x.(fieldNames{ii}))
            x.(fieldNames{ii}) = x.(fieldNames{ii})(indices,:);
    end    
end

