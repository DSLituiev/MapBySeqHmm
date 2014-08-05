function varargout = plotWaitingTimes(tv)
T_max = max(tv);

tau = T_max/size(tv,1);
dt = 0.25*tau;
ti = (0.5*dt):dt:(16*tau); %<= plotting times
f = figure('name', 'distribution of "waiting time"');
pEmp = hist(abs(diff(tv)), ti );

pTheor =  size(tv,1)*(exp( - (ti-.5*dt)./tau) - exp( - (ti + .5*dt)./tau) );

bar(ti, pEmp)
hold all
plot(tau, interp1(ti,pTheor,tau), 'rd', 'markersize', 5, 'MarkerFaceColor', 'r')
plot(ti, pTheor, 'r-')
xlim([0, (15*tau)])
xlabel('recombination distance'); ylabel('frequency');

if nargout>=1
    varargout = {tau, f, pEmp};
end