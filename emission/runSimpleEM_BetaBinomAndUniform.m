function [theta, lambda, varargout] = runSimpleEM_BetaBinomAndUniform(Pf, q, r, N, varargin)
% runs EM to unmix Beta-Binomial and Uniform Distributions treating each
% quasi-binomially sampled locus independently;
% assumes that the underlying state variable 'f = 0, 1/(2N), ... 1/2'
% has some distribution 'Pf'
% p(f_k) = 2^(-N).*nchoosek(N, k);
%

f = 0:1/(2*N):1/2;
%= P(f| stationary model)
% Pf = pop.Pstat;

%% check the input parameters
p = inputParser;
addRequired(p, 'q', @isnumeric);
addRequired(p, 'r', @isnumeric);
addRequired(p, 'N', @(x)isnumeric(x)&isscalar(x) );
%
addOptional(p,'verbose', '', @(x)(ischar(x)|isempty(x)) );
%
addParamValue(p,     'errTol',           1e-6, @isscalar);
addParamValue(p,     'theta',         1e-2, @isscalar);
addParamValue(p,     'lambda',        1-1e-2, @isscalar);
addParamValue(p,     'maxIter',        1e3, @isscalar);
addParamValue(p,     'contribution',     1, @isnumeric);

parse(p, q, r, N, varargin{:});
%%


if ~isempty(p.Results.verbose) && strcmpi(p.Results.verbose, 'v')
    verboseFlag = true;
else
    verboseFlag = false;
end

if ~isempty(p.Results.verbose) && strcmpi(p.Results.verbose, 'b')
    simpleFlag = true;
else
    simpleFlag = false;
end

theta = p.Results.theta;
lambda = p.Results.lambda;

if numel(p.Results.contribution) == numel(q)
   contr = p.Results.contribution;
else
   contr = 1;    
end
% %%
%     theta = [0, 10.^(-5:0.25:5) ];
%     [TH, F] = ndgrid(f, theta');
%     [LL ] = logBetaBinomialThetaMu0(q, r, F, TH);
%     Lf = squeeze(sum(sum(bsxfun(@plus, reshape(LL, [numel(q), numel(f), numel(theta)]), log10(Pf) ), 2),1));
%     [~, ii] = max( LL );
%     % any( imag(reshape(LL, [numel(q), numel(f), numel(theta)]))>0 )
%     
    %% if max LL (theta=0), put thetaEst = 0
    %= as the LL-derivative does not converge to zero from above
%     if (ii == 1)
%         thetaEst = 0;
%         return
%     end
%     %% iterate second time
%     theta = theta(ii).*10.^(-0.25:0.025:0.25);
%     [~,  RR ] = maxLlThetaBetaBinomial(k, n, mu, theta, true);
%     [~, ii] = min( abs(RR) );
%     
%     theta0 = theta(ii);
%% initialize
dLambda = 1;
ii = 0;
options = optimset('TolX', p.Results.errTol, 'FunValCheck', 'on', 'display', 'off'); % set TolX
% options = optimset('TolX', p.Results.errTol, 'FunValCheck', 'on'); % set TolX

if verboseFlag
    fprintf('iter.# |\ttheta |\t lambda |\n')
    msg = sprintf('%u\t%4.4e\t%4.3f\n', ii, theta, lambda);
    fprintf(msg)
    msgLength = numel(msg);
end

if ~simpleFlag
    jj = 0;
while (abs(dLambda) > p.Results.errTol) && ii< p.Results.maxIter
    ii = ii+1;
    [pBB, pU] = conditBetaBinomStat(q, r, theta, N, f, Pf);   
%   gi1 = contr.*lambda.*pBB./(contr.*lambda.*pBB + (1-contr.*lambda).*pU );
    gi1 = lambda.*pBB./(lambda.*pBB + (1 - lambda).*pU );
     
    clear pBB pU
    
    lambdaNew = nanmean(gi1);
    
    % [ dLdTheta] = ddTheta_LogBetaBinomialThetaMu0(q, r, f, theta);
    % sum( gi1.*sum(bsxfun(@times, dLdTheta, Pf), 2) , 1)
    
    fhRes = @(th)nansum( contr.*gi1.* ...
        nansum(bsxfun(@times,  ddTheta_LogBetaBinomialThetaMu0(q, r, f, th), Pf), 2) , 1);
    
%     nansum( contr.*gi1.*nansum(bsxfun(@times,  ddTheta_LogBetaBinomialThetaMu0(q, r, f, theta), Pf), 2) , 1)
    
    % figure
    % surf( ddTheta_LogBetaBinomialThetaMu0(q, r, f, theta) )
    
    % feval(fhRes, theta)
    % dLdT = ddTheta_LogBetaBinomialThetaMu0(q, r, f, theta);
    % nansum(bsxfun(@times,  dLdT, Pf), 2)
    % [~, sumdLdTheta] = maxLlThetaBetaBinomial(q, r, f, theta);
    % sum(dLdT,1) - sumdLdTheta
    
    [thetaEst, ~, ~] = newtonraphson(@(th)fhRes(th), theta, options);
    
    dLambda = lambdaNew - lambda;
    
    theta = abs(thetaEst);
    lambda = abs(lambdaNew);
    
%= if hits negative values of $theta$, increase theta and decrease lambda:
    if thetaEst<0 || ~isreal(thetaEst) 
         jj = jj+1;
         lambda = 1- (1-p.Results.lambda)*2^(jj/8);
         theta = p.Results.theta*2^(jj/8);
     end
     
    if verboseFlag
          fprintf(repmat('\b',1, msgLength));
          msg = sprintf('%u\t%4.4e\t%4.3f\n', ii, thetaEst, lambdaNew);
          fprintf(msg);
          msgLength = numel(msg);
%         fprintf('%u\t%4.4e\t%4.3fi\n', ii, theta*1i, lambda*1i)    
    end
end
else
     fhRes = @(th)nansum( contr.* ...
        nansum(bsxfun(@times, ddTheta_LogBetaBinomialThetaMu0(q, r, f, th), Pf), 2) , 1);
     [theta, ~, ~] = newtonraphson(@(th)fhRes(th), theta, options);
end


if ii == p.Results.maxIter
    warning('runSimpleEM_BetaBinomAndUniform:UnConverged',...
        'Convergence problem in lambda1! Cut on the 1000 iteration');
end

if nargout>2
    varargout{1} = gi1;
end
