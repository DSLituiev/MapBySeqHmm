function [ f ] = plotSeqSim( obj )
%PLOTSEQSIM Summary of this function goes here
%   Detailed explanation goes here

if obj.Selection
    fname = 'chromosome under selection';
else
    fname = 'stationary chromosome';
end

f = figure('name', fname);

stairs(obj.tChr, obj.kChr, 'b', 'linewidth', 2.5)
hold on
plot(obj.t, obj.k_t, 'rs-','markersize',8,'MarkerFaceColor','r')
plot(obj.t, 2*obj.pop.N * obj.q./obj.r, 'go','markersize',8,'MarkerFaceColor','g')
xlim([0, obj.T_max])
ylim([0, 10*round(.15 * obj.pop.N)])
xlabel('linkage, morgans')
legend({'k: chromosome chunks'; 'k: SNP genotyping'; 'f: mutant read frequency x 2\itN'})

end

