% test
cd '/media/Processing/MATLAB_Dima/SeqMapping_20140519/c'

addpath('/media/Processing/MATLABuserfunctions/mtimesx');


E = [ 1 1 1; 2 2 2; 6 3 6];

T = cat(3, 0.1*ones(3,3), 0.5*ones(3,3));

% mtimesx(T)

mex -v -largeArrayDims mexCrossMatr.c CFLAGS="\$CFLAGS -std=c99"

A = mexCrossMatr(E, T);

A0 =  bsxfun(@times,  permute( E, [2, 3, 1] ), T ) ;

% 

mex -v mexCumMatrSafe.c CFLAGS="\$CFLAGS -std=c99"

E = y.E;
A = y.A;
T = y.T;

clear logAlpha logBeta

for tt = 1:20
tic

A =  bsxfun(@times,  permute( E(1:end-1, :), [2, 3, 1] ), T ) ;

M = size(E,1);
Ac = E( end, :)';
logAlpha(M, :) = log10(Ac');

scalePrev = 0;
scaleA = zeros(M, 1);
% for-loop frwd
for m = (M-1):-1: 1
    Ac = A(:,:, m) * (Ac * 10.^(scaleA(m+1) - scalePrev) );
    logAlpha(m, :)= log10(Ac) - scaleA(m+1);
    scalePrev = scaleA(m+1);
    scaleA(m) = - max( logAlpha( m, :) );
end

Np = size(E,2);
%% backward
Bc = ones(1, Np);
logBeta(1, 1: Np) = 0;

scalePrev = 0;
scaleB = zeros(M-1, 1);

for m = 2:1:M
    Bc = Bc * A(:,:, m-1)' * (10.^(scaleB(m-1) - scalePrev));
    logBeta(m, :) = log10(Bc) - scaleB(m-1);
    scalePrev = scaleB(m-1);
    scaleB(m) = - max(logBeta(m, :));
end


elt(tt) = toc;
end

mean(elt)

for tt = 1:20
   tic
   [logAlphaN, logBetaN, AN] = mexCumMatrSafe(E,T);
   eltN(tt) = toc;    
end
mean(eltN)


mex -v mexCumMatrSafe.c CFLAGS="\$CFLAGS -std=c99"