function f = plotChromosome2D(obj, chr)
%PLOTCHROMOSOME2D Summary of this function goes here
%   Detailed explanation goes here

Z = bsxfun( @minus, obj.xkPflat , calcMarginal( obj.xkPflat , 2) );

f = figure;
surf( obj.x(obj.chromosome == chr)*1e-6, obj.pop.kvect, Z(obj.chromosome == chr,:)', 'linestyle', 'none');
view(0,90)

end

