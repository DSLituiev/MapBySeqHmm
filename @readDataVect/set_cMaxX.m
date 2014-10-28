function [ obj ] = set_cMaxX(obj, filePath)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

assert(isprop(obj, 'cMaxX'))

readOut = dlmread(filePath, '\t', 0, 1);
obj.cMaxX = uint32(readOut(:,1));

end

