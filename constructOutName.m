function [ dataID ] = constructOutName( dataID0, varargin)
%CONSTRUCTOUTNAME Summary of this function goes here
%   Detailed explanation goes here

if nargin>1
    linkageLoosening = varargin{1};
    if isscalar(linkageLoosening)
        dataID = sprintf(strcat(dataID0, '_%u'), linkageLoosening);
    else
        dataID = sprintf(strcat(dataID0, '_%3f_%3f'), linkageLoosening(1), linkageLoosening(end));
    end
else
    dataID = dataID0;
end

end

