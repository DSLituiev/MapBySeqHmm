function [ f ] = plotHmmInferenceRes( obj )
%PLOTHMMINFERENCERES Summary of this function goes here
%   Detailed explanation goes here

fprintf('t0 predicted:\t%g\n', obj.t0_predicted)
fprintf('t0 set      :\t%g\n', obj.t0)

f = figure('name', 'inference results: log-odds');
stem(obj.t, obj.xLogOdds,  'go', 'BaseValue', floor(min(obj.xLogOdds))-2,...
    'markersize',3,'MarkerFaceColor','g')
hold on
plot(obj.t, obj.xLogOdds, 'y-')
plot(obj.t0_predicted, obj.maxLogOdds, 'bv')
plot(obj.t0, max(obj.xLogOdds), 'r^')

if any(obj.xLogOdds<0)
    plot( obj.t([1,end]) , [0,0] , 'k--' )
end

xlim([0, obj.T_max])
ylim([ floor(min(obj.xLogOdds(~isinf(obj.xLogOdds)))), ceil(max(obj.xLogOdds))] )


end

