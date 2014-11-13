function [ msgLength ] = updateLog( msgLength, nn, varargin )
%UPDATELOG Summary of this function goes here
%   Detailed explanation goes here
if nargin > 2  && ischar(varargin{1})
    message = varargin{1};
else
    message = 'iterations left :\t%u\t...';
end

fprintf(repmat('\b', 1, msgLength));
textA = sprintf(message, nn);
fprintf(textA)
msgLength = numel(textA);

end

