clear all
close all

par.N = 20;
tau = 0.33;

Q = Qmatrix(par.N);
T = expm(Q*tau);
p = StationaryDistr(par.N);


norm(p)

norm(T*p - p)

norm( p'*T - p' )

figure
plot(0:1:par.N, p)


J = diag(p)*T;

norm(J - J')

% p0 = eig(T1);

figure
pcolor(J)
axis  equal tight

%% stationary

Ep = bsxfun(@rdivide, E, p) ;
Ep(:,[1, end] ) = E(:, [1, end] );

% J =   bsxfun(@mtimes, diag(p), T) ;


addpath('D:\MATLABuserfunctions\mtimesx');

J =  mtimesx( diag(p), T);
%% back joint
b = J(:,:,end) * Ep(:,end);

for m = size(T,3):-1:2
%     A1 = bsxfun(@times, E(:,m), J(:, :, m) );
%     A0 = diag( E(:,m) ) * J(:, :, m) ;
%     norm(A1-A0)
    b =  J(:, :, m-1) * diag( Ep(:,m) ) *  b;
end

pBJ = Ep(:,1)'*b

%% forth joint
f = ( Ep(:,1)' * J(:,:,1)) ;

for m = 2:1:size(T,3)
    f =  f* diag( Ep(:,m) ) * J(:, :, m) ;
end

pFJ = f* Ep(:,end)

%% back conditional
b = J(:,:,end) * E(:,end);

for m = size(T,3):-1:2
%     A1 = bsxfun(@times, E(:,m), J(:, :, m) );
%     A0 = diag( E(:,m) ) * J(:, :, m) ;
%     norm(A1-A0)
    b =  (diag( E(:,m) ) * T(:, :, m-1))' *  b;
end

pBC = E(:,1)'*b
%% forth conditional
f = ( J(:,:,1) * E(:,1) )';
% f = ( E(:,1)' * J(:,:,1)) ;

for m = 2:1:size(T,3)
    f =  f * diag( E(:,m) ) * T(:, :, m) ;
end

pFC = f* E(:,end)
%%
J0 = J(:,:,1);
norm(J0 - J0' )

figure
pcolor(J0)
axis  equal tight


% norm( J(:,:,end)* Ep(:,end) -  ( Ep(:,end)' * J(:,:,end))' )

% norm( J(:, :, m) * diag( E(:,m) ) - (J(:, :, m) * diag( E(:,m) ))' )

