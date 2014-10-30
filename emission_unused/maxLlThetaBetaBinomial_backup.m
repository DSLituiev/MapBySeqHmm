function [ dLdTheta0, LL0,  dLdTheta_s, LL_s , L0] = maxLlThetaBetaBinomial(k, n, mu, theta, varargin)

N = numel(k);
% Res_n  = zeros(N,1); Res_n_k = zeros(N,1); Res_k = zeros(N,1);
% LL_n   = zeros(N,1); LL_n_k  = zeros(N,1);  LL_k = zeros(N,1);
% L_n    = zeros(N,1);  L_n_k  = zeros(N,1);   L_k = zeros(N,1);
dimK = find(size(k)>1);
if isempty(dimK)
    dimK = find(size(theta)==1);
end
    
dims = size(bsxfun(@plus, theta, k));
S.type = '()'; 
S.subs = cellstr(repmat(':', [numel(dims),1]));

Res_n = zeros(dims); LL_n = zeros(dims); L_n = zeros(dims);
Res_k = zeros(dims); LL_k = zeros(dims); L_k = zeros(dims);
Res_n_k = zeros(dims); LL_n_k = zeros(dims); L_n_k = zeros(dims);

for ii = numel(n):-1:1
    S.subs{dimK} = ii;
    
    r_n   = permute(0:(n(ii)-1), [1,3,2]);
    [Res_1, LL_1 ,  L_1] =  mySeries(r_n, 1, theta);
    
    Res_n = subsasgn(Res_n, S, Res_1);
    LL_n  = subsasgn(LL_n, S, LL_1);
    L_n   = subsasgn(L_n, S, L_1); 
     
    r_k   = permute( 0:(k(ii)-1), [1,3,2]);
    [Res_1, LL_1 ,  L_1] =  mySeries(r_k, mu, theta);
    Res_k = subsasgn(Res_k, S, Res_1);
    LL_k  = subsasgn(LL_k, S, LL_1);
    L_k   = subsasgn(L_k, S, L_1); 
     
    
    r_n_k =  permute( 0:(n(ii)-k(ii)-1), [1,3,2]);
    [Res_1, LL_1 ,  L_1] =  mySeries(r_n_k, 1 - mu, theta);
    Res_n_k = subsasgn(Res_n_k, S, Res_1);
    LL_n_k  = subsasgn(LL_n_k, S, LL_1);
    L_n_k   = subsasgn(L_n_k, S, L_1); 
%     Res_n(:, ii)   = sum( bsxfun(@rdivide, r_n, (1 + bsxfun(@times, r_n, theta))) , 2);
%     LL_n(:, ii)    = sum( log10( 1 + bsxfun(@times, r_n, theta) ) , 2);
%     L_n(:, ii)     = prod( ( 1 + bsxfun(@times, r_n, theta) ) , 2);
     
%     Res_k(:, ii)   = sum( bsxfun(@rdivide, r_k, ( mu + bsxfun(@times, r_k, theta) ) ), 2);
%     LL_k(:, ii)    = sum( log10( mu + bsxfun(@times, r_k, theta) ) , 2);
%     L_k(:, ii)     = prod( ( mu + bsxfun(@times, r_k, theta) ) , 2);

%     Res_n_k(:, ii)  = sum( bsxfun(@rdivide, r_n_k, bsxfun(@plus, (1-mu), bsxfun(@times, r_n_k, theta))), 2);
%     LL_n_k(:, ii)   = sum( log10( (1-mu) + bsxfun(@times, r_n_k, theta) ) , 2);
%     L_n_k(:, ii)    = prod( ( (1-mu) + bsxfun(@times, r_n_k, theta) ) , 2);
end
% residual:
dLdTheta0 = Res_k + Res_n_k - Res_n;
LL0  = bsxfun(@plus, LL_k   + LL_n_k   - LL_n , log10(binomial(n(:),k(:))) );
L0   = bsxfun(@times, binomial(n(:), k(:)), L_k.*L_n_k ./ L_n );

% Res0 = sum(Res0);
% L0 = sum(L0);


    
rr = permute(0:(max(n)-1), [1,3,2]);

S_k = cumS(k, rr, dimK);
S_n_k = cumS(n-k, rr,dimK);
S_n = cumS(n, rr, dimK);

rrTheta = bsxfun(@times, rr, theta);
% R0 = bsxfun(@times, rr, ...
%     bsxfun(@rdivide, (N - S_k),   (mu + rrTheta)     ) + ...
%     bsxfun(@rdivide, (N - S_n_k), (1-mu + rrTheta) ) - ...
%     bsxfun(@rdivide, (N - S_n ),  (1 + rrTheta)     )  );

% figure; surf(0:1:(n-1) ,theta, R0)
% figure; plot(0:1:(n-1) , R0)

% squeeze(S_k)

dLdTheta_s = sum( bsxfun(@times, rr, ...
    bsxfun(@rdivide, ( S_k ),   (mu + rrTheta)     ) + ...
    bsxfun(@rdivide, (S_n_k), (1-mu + rrTheta) ) - ...
    bsxfun(@rdivide, ( S_n ),  (1 + rrTheta)     )  ), 3 );

%  SS_k   =  sum( bsxfun(@times, (S_k)  , log10(    mu + rrTheta) ) , 3 );
%  SS_n_k =   sum( bsxfun(@times, (S_n_k)  , log10(  1-  mu + rrTheta) ) , 3 );
%  SS_n   =   sum( bsxfun(@times, (S_n)  , log10(    1 + rrTheta) ) , 3 );
% 
% norm( squeeze( SS_k ) - (LL_k), 2)
% norm( squeeze(SS_n_k ) - (LL_n_k) , 2)
% norm( squeeze(SS_n ) - (LL_n) , 2)

LL_s =  bsxfun(@plus ,...
    sum( ...
    bsxfun(@times, ( S_k)  , log10(    mu + rrTheta) ) + ...
    bsxfun(@times, ( S_n_k), log10(1 - mu + rrTheta) ) - ...
    bsxfun(@times, ( S_n ) , log10(1      + rrTheta) ) ,...
     3 ), sum( log10(binomial(n,k)), dimK ) ) ;


    function S_x = cumS(x, rr, dimK)
        F_x = uint32(bsxfun(@le, rr, x-1));
        S_x = sum(F_x, dimK);
%         S_x = (cumsum(g_x));
    end

    function [R, LL, L] = mySeries(x, m, t)
        R  = sum( bsxfun(@rdivide, x, bsxfun(@plus, m, bsxfun(@times, x, t))), 3);
        LL = sum( log10( m + bsxfun(@times, x, t) ) , 3);
        L  = prod( (  m + bsxfun(@times, x, t) ) , 3);
    end

end
