
N = 50;

kvect = 0:1:N;
kvectFull = 0:1:2*N;

% Pz = ones(1, N+1)./(N+1);

Pz = factorial(N)./factorial(N-kvect)./ factorial(kvect)*2^-N;

r = 30;
q = (0:1:r)';

B = bsxfun( @(x,y)binopdf(x,r,y), q, kvectFull/(2*N));

Px_y1 = sum(bsxfun(@times, B(:, kvect+1), Pz ), 2) ;
Px_y0 = sum(B , 2)./(2*N+1);

figure
plot(q, Px_y0 )
hold all
plot(q, Px_y1 )

lambda = 0.5;

figure
plot(q, lambda*Px_y1./((1-lambda)*Px_y0 + lambda*Px_y1) )
hold all



figure
pcolor(B)