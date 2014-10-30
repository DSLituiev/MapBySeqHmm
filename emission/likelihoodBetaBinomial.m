function [ L, logLtot, Var, phi, L0   ] = likelihoodBetaBinomial(k, n, b, mu, varargin)

dimK = find(size(k)>1);
% b = b(:);
a = bsxfun(@times, b, mu./(1-mu));
b = bsxfun(@times, b, ones(size(mu)));
a_b = bsxfun(@times, b, 1./(1-mu));

k_a = bsxfun(@plus, k, a);
n_k_b = bsxfun(@plus, n - k, b);

L =  bsxfun(@beta, k_a, n_k_b).*  bsxfun(@rdivide, binomial(n,k), bsxfun(@beta, a, b) );

Var =  bsxfun(@rdivide, ...
    bsxfun(@times, n, a.*b).*(bsxfun(@plus, a_b, n) ),...
    ( (a_b).^2 .* (a_b+1) ) );

% Var =  bsxfun(@rdivide, ...
%     bsxfun(@times, bsxfun(@times, n, 1-t) ,  bsxfun(@plus, n, a_b)),...
%     a_b + 1 );

phi = 1./a_b;% bsxfun(@rdivide,  (1-t), b );

if nargin>4 && ~isempty(varargin{1})
    mult = varargin{1};
else
    mult = 1;
end

if ~isempty(dimK)
    valInds = ~isinf(log10(L)) ;
    for d = ndims(L):-1:1
        if d~= dimK
            valInds = all(valInds, d) ;
        end
    end
    S.type = '()'; 
    S.subs = cellstr(repmat(':', [ndims(L),1]));
    S.subs{dimK} = valInds;
%     logLtot = nansum(log10(subsref(L,S)), dimK);
    inds = ~isinf(log10(L))& ~isnan(log10(L));
    logLtot = nanmean(log10(bsxfun(@times, L, mult)).*inds, dimK);
    logLtot = logLtot*size(L, dimK);
else
    logLtot = log10(L.*mult);
end


% beta = mu./theta;
theta = mu./a;
% N = numel(k);
% L_n    = zeros(N,1);  L_n_k  = zeros(N,1);   L_k = zeros(N,1); 

for ii = numel(n):-1:1
    r_n   = 0:(n(ii)-1);
    L_n(:,ii)   = prod( ( 1 + bsxfun(@times, r_n, theta) ) , 2); 
%     L_n1(ii,1)   = prod( ( 1./theta + r_n)  , 2); 
    
    r_k   = 0:(k(ii)-1);
    L_k(:,ii)    = prod( ( mu + bsxfun(@times, r_k, theta) ) , 2);
%     L_k1(ii,1) = prod(  mu./theta + r_k  , 2);
    
    r_n_k = 0:(n(ii)-k(ii)-1);
    L_n_k(:,ii)   = prod( ( 1- mu + bsxfun(@times, r_n_k, theta) ) , 2);
%     L_n_k1(ii,1)  = prod(  (1-mu)./theta + r_n_k  , 2);
end

L0   = bsxfun(@times, binomial(n,k), L_k.*L_n_k ./ L_n );


% L0   = binomial(n,k).* L_k.*L_n_k ./ L_n ;
% 
% L_n0 = gamma(n + 1./theta)./gamma(1./theta);
% mean( (L_n0 - L_n.*theta.^-n).^2)
% 
% L_k0 = gamma(k + mu./theta)./gamma(mu./theta);
% mean( (L_k0 - L_k.*theta.^-k ).^2)
% 
% L_n_k0 = gamma(n-k + (1-mu)./theta)./gamma((1-mu)./theta);
% mean( (L_n_k0 - L_n_k.*theta.^-(n-k) ).^2)
% 
% %====
%  mean( (L_k0.*L_n_k0./L_n0  - ...
%     gamma(k+a).* gamma(n-k+b)./gamma(n+a+b) ./ (gamma(a).* gamma(b)./gamma(a+b))...
%     ).^2) 
% 
% 
% %====
% mean( (L_k0.*L_n_k0./L_n0 -  L_k.*L_n_k./L_n ).^2)
% 
% mean( ( bsxfun(@rdivide, bsxfun(@beta, k_a, n_k_b) , bsxfun(@beta, a, b) ) - ...
%      L_k.*L_n_k./L_n ).^2)
% 
% mean( ( bsxfun(@rdivide, bsxfun(@beta, k+a, n-k+b) , bsxfun(@beta, a, b) ) - ...
%      L_k.*L_n_k./L_n ).^2)
%   
%   
% mean( ( bsxfun(@beta, k+a, n-k+b) - gamma(k+a).* gamma(n-k+b)./gamma(n+a+b)  ).^2) 
%  
% mean( ( bsxfun(@beta, a, b) - gamma(a).* gamma(b)./gamma(a+b) ).^2) 
% 
% 
% 
% mean( ( bsxfun(@rdivide, bsxfun(@beta, k+a, n-k+b) , bsxfun(@beta, a, b) ) - ...
%      L_k.*L_n_k./L_n ).^2) 
%  
%  
% mean( ( bsxfun(@rdivide, bsxfun(@beta, k+a, n-k+b) , bsxfun(@beta, a, b) ) - ...
%     gamma(k+a).* gamma(n-k+b)./gamma(n+a+b) ./ (gamma(a).* gamma(b)./gamma(a+b))...
%     ).^2) 
% 
%  