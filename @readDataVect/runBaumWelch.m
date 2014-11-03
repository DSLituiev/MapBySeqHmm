function obj = runBaumWelch(obj, varargin)
% maxIter = 100;
% chr = 1;

warning('runBW:info', [' TO DO : the SNPs are selected mostly based on the distribution,', ...
    ' not much for their linkage'])
%% check the input parameters
p = inputParser;
addRequired(p, 'obj', @isobject);
%
addOptional(p, 'chr', 0, @(x)isnumeric(x)&isscalar(x) );
addOptional(p,'verbose', '', @(x)(ischar(x)|isempty(x)) );
%
addParamValue(p,     'errTol',           1e-6, @isscalar);
addParamValue(p,     'maxIter',           10, @isscalar);


parse(p, obj, varargin{:});
%% check pre-calculated matrices

if any(isnan(obj.contrib))
    warning('runBW:NaN_contrib', 'NaNs in the contributions');
    obj.contrib(isnan(obj.contrib)) = 1;
end
%% initial values:
% Pzm = obj.pop.Pstat;
% options = optimset('TolFun', p.Results.errTol.*obj.pop.Np, 'FunValCheck', 'on', 'display', 'off'); % set TolX

dLambda = NaN(p.Results.maxIter+1, 1);


Pzm = 10.^logBetaBinomialThetaMu0(obj.pop.kvect(:), obj.pop.N, 2/3, 2e-1)';
% binopdf(obj.pop.kvect, obj.pop.N, 0.33);
Pzm = Pzm./sum(Pzm);
figure('name', 'the initial distribution assumed for demixing')
plot(obj.pop.kvect, Pzm);

% ./obj.pop.Np;
lambda0 = 0.8 ;

% obj.contrib(:) = lambda0 + 0.2*(rand(obj.Mtot, 1) - 0.5);
chrV = p.Results.chr;
if ~chrV
    chrV = 1:obj.chrNumber;
end

for cc = chrV
    
    lambdaOld = obj.contrib(obj.ci{cc});
    Px_zm_y1 = zeros(obj.M(cc), obj.pop.Np);
    Px_zm_y0 = zeros(obj.M(cc), obj.pop.Np);
    
    obj.calcT(cc, obj.Alpha(1));
    %%
    
    % dLambda(1) = 1;
    lambda0EstLog = -Inf(obj.M(cc), 1);
    lambda1EstLog = -Inf(obj.M(cc), 1);
    
    for ii = 1 : p.Results.maxIter
        
        [mE1, mE0] = mixBetaBinomUniform(double(obj.q), double(obj.r),  obj.pop.N, 1e-2);
        mE0 = repmat(mE0, [1, obj.pop.Np]);

        %% P [x_m | z_m, lambda_m]
        % lambda_m = obj.contrib
        obj.E = bsxfun(@plus,...
            bsxfun(@times, mE1, obj.contrib), ...
            bsxfun(@times,mE0, (1-obj.contrib) ) ...
            );
        
        % obj.calcEmission;
        obj = obj.crossMatr(cc);
        obj = obj.cumMatr(cc);

        %% P[ x_ | z_m, lambda_m ]
        %= P[ x_  | z_m, y_m ]
        rawFB = obj.logAlpha{cc} + obj.logBeta{cc} - log10(obj.E(obj.ci{cc}, :)) ;
        % any(isnan(rawFB),2)
        rawFB(isnan(rawFB(:,1)), 1) = -Inf;
        Px_zm_y1 = rawFB + log10(mE1(obj.ci{cc}, :)) ;
        Px_zm_y0 = rawFB + log10(mE0(obj.ci{cc}, :)) ;
        notNan = ~any(isnan(rawFB),2);
        
        logQ_zm = bsxfun(@times,  lambdaOld(notNan) , Px_zm_y1(notNan, :)) + ...
            bsxfun( @times, ( 1 - lambdaOld(notNan) ) , Px_zm_y0(notNan, :));        
        logQ_zm = bsxfun(@minus, logQ_zm, calcMarginal(logQ_zm, 2));
        Q_zm = 10.^logQ_zm;
        Q_zm = bsxfun( @rdivide, Q_zm, sum(Q_zm, 2));
        
        %     figure; plot(obj.x(obj.ci{chr}), calcMarginal(Px_zm_y1(obj.ci{chr},  :) - Px_zm_y0(obj.ci{chr}, :), 2), 'x' )
        
        %% Q[y]
        % Ty = ones(2,2)./2;
        lambda1EstLog(notNan) =   log10(lambda0)     + sum(Q_zm(notNan, 2:end) .*  Px_zm_y1(notNan, 2:end), 2)   ;
        lambda0EstLog(notNan) =   log10(1 - lambda0) + sum(Q_zm(notNan, 2:end) .*  Px_zm_y0(notNan, 2:end), 2)   ;
        lambdaNew = 10.^(lambda1EstLog - calcMarginal([lambda0EstLog, lambda1EstLog],2));
        lambdaNew(isnan(lambdaNew)) = 0.5;
             
        %     figure; plot(obj.x(obj.ci{chr}), lambdaEst);
        %     figure; plot(lambdaEst);
        
        dLambda(ii) = norm(lambdaOld- lambdaNew)./numel(lambdaNew);
        
        
        if any(isnan(lambdaNew))
            warning('runBW:NaN_contrib', 'NaNs in the contributions. Replacing with 0.5');
            obj.contrib(isnan(lambdaNew)) = 1;
        end
        lambdaOld = lambdaNew;
        
        
        %       Pzm = calcMarginal( cat(3, bsxfun(@plus, Pz_y0, 1-obj.contrib(obj.ci{chr})) ,...
        %            bsxfun(@plus, Pz_y1, obj.contrib(obj.ci{chr}))  ), 3 );
        %
        %       calcMarginal(Pzm( :, 2:end),2)
        
        if dLambda(ii) < 1e-6 || all(isnan(Q_zm(:))) || ii == p.Results.maxIter
            fprintf('chromosome:\t%u\t# iterations:\t%u\tchange:\t%e \n', cc, ii, dLambda(ii) )
            break
        end
    end    
    obj.contrib(obj.ci{cc}) = lambdaOld;
end

%  figure
%  surf(obj.x(obj.ci{chr}), obj.pop.kvect, Pzm', 'linestyle', 'none'); view( 0 , 90 )

if p.Results.chr~=0
    figure('name', 'Q[z_m]')
    surf(obj.x(obj.ci{cc}), obj.pop.kvect, Q_zm', 'linestyle', 'none'); view( 0 , 90 )
    
    figure('name', 'b-b')
    surf(obj.x(obj.ci{cc}), obj.pop.kvect, Px_zm_y1', 'linestyle', 'none'); view( 0 , 90 )
    
    
    figure('name', 'uniform')
    surf(obj.x(obj.ci{cc}), obj.pop.kvect, Px_zm_y0', 'linestyle', 'none'); view( 0 , 90 )
    
    figure('name', 'uniform: indexwise')
    p = pcolor(1:obj.M(cc), obj.pop.kvect, Px_zm_y0'); view( 0 , 90 ); set(p,  'linestyle', 'none');
    
    
    % i0 = find(obj.ci{chr}, 1,'first');
    % ii= 30;
    % figure
    % plot(obj.pop.kvect, Px_zm_y0(ii, :)' - calcMarginal( Px_zm_y0(ii, :),2) );
    % hold all
    % plot(obj.pop.kvect, log10(obj.E(i0+ii, :))' );
    % legend({'P[x|z_m, y_m = 0]', 'E'})
    
    
    
    
    indsTrue = obj.contrib > 0.5;
    
    figure;
    plot(obj.x(indsTrue& obj.ci{cc}), obj.f(indsTrue& obj.ci{cc}), 'gx')
    hold all;
    plot(obj.x(~indsTrue& obj.ci{cc}), obj.f(~indsTrue& obj.ci{cc}), 'rx')
    
    
    figure('name' , 'membership on the chromosome');
    plot(obj.x( obj.ci{cc}), obj.contrib( obj.ci{cc}), 'bx-')
    
end