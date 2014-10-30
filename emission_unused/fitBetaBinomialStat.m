function thetaEst = fitBetaBinomialStat(k, n, mu, varargin )

if nargin>3 && ~isempty(varargin{1}) && isscalar(varargin{1})
    theta0 = varargin{1};
else
    %% find maximum log-likelihood on a grid
    theta = [0, 10.^(-5:0.25:5) ];
    [~,  ~, LL ] = maxLlThetaBetaBinomial(k, n, mu, theta, true);
    [~, ii] = max( LL );
    
    %% if max LL (theta=0), put thetaEst = 0
    %= as the LL-derivative does not converge to zero from above
    if (ii == 1)
        thetaEst = 0;
        return
    end
    %% iterate second time
    theta = theta(ii).*10.^(-0.25:0.025:0.25);
    [~,  RR ] = maxLlThetaBetaBinomial(k, n, mu, theta, true);
    [~, ii] = min( abs(RR) );
    
    theta0 = theta(ii);
end
fprintf('theta_0 = \t%2.4g', theta0);

options = optimset('TolX',1e-8); % set TolX
[thetaEst, ~, ~] = newtonraphson(@(x) maxLlThetaBetaBinomial(k, n, mu, x), theta0, options);
