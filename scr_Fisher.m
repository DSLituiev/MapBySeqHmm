


addpath('D:\MATLABuserfunctions\mtimesx');

N = 10;
study = population(N);
kvect = 0:1:N;
tvect = 1./(0:0.01:1);

R = 10:10:100;
for rr= numel(R):-1:1
    T = transition3D_expm(N, tvect);
    
    emissionHandle = @(qq, rr, ff)emissionk0(qq, rr, study);
    
    AR.q = (1:1:R(rr))';
    AR.r = R(rr)*ones(R(rr),1);
    
    [AR.E, AR.c] = wrapEmissionMatrix(AR.q, AR.r, study, emissionHandle);
    
    AR.E = bsxfun(@rdivide, AR.E, sum(AR.E,2) );
    
    a = mtimesx( permute(AR.E, [4, 2, 3, 1]), T);
    % [1 x N+1, M, r1]
%     size(a)
    
    a = sum(a, 4)/5;
    
    size(a)
    
    
%     sum(a,2)
    
    a = sum( mtimesx( a, permute(AR.E, [2, 4, 3, 1])),4);
    
%     size(a)
    
    a = squeeze(a);
    
    p(:,rr) = a;
end


figure
plot(1./tvect,  p')

figure
pcolor(1./tvect, R, p')
% plot(1./tvect, (a- a(end))./(a(1)- a(end)) )
% % set(gca, 'yscale', 'log')
% set(gca, 'xscale', 'log')


figure
pcolor(1./tvect(2:end-1), R, (diff(p,2,1))')


figure
plot(1./tvect(2:end-1), sqrt(diff(a, 2)) )
set(gca, 'yscale', 'log')

e1 = ones(N+1,1)./(N+1);
eM = e1;% zeros(N+1,1);
% eM(end) = 1;

b = mtimesx( mtimesx(e1', T), eM);
b= squeeze(b);

figure
plot(1./tvect, b)
% plot(tvect, (b- b(end))./(b(1)- b(end)) )
set(gca, 'yscale', 'log')
set(gca, 'xscale', 'log')

figure
plot(1./tvect(2:end-1), sqrt(diff(b, 2)) )
set(gca, 'yscale', 'log')

% tvect