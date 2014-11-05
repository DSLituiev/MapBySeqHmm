function f_out = plotProbabilitiesOnAChromosome(obj, chr)

f_out = figure;

subplot(2,1,1)
plot(obj.x(obj.ci{chr}), obj.xPstat(obj.ci{chr}), 'g-')
hold all
plot(obj.x(obj.ci{chr}), obj.xPout(obj.ci{chr}), 'r-')
plot(obj.x(obj.ci{chr}), obj.xPflat(obj.ci{chr}), 'b-')
legend({'stationary', 'selection', 'flat'})

subplot(2,1,2)
plot(obj.x(obj.ci{chr}), obj.xPstat(obj.ci{chr})  + obj.cNormConst(chr) , 'g-')
hold all
plot(obj.x(obj.ci{chr}), obj.xPsel(obj.ci{chr}) + obj.cNormConst(chr) , 'r-')
legend({'stationary norm',  'selection norm'})

