function obj = runBW(obj, chr)
maxIter = 100;
% chr = 1;
%% check pre-calculated matrices

%% initial values:
Pzm = repmat(log10(obj.pop.Pstat), [ obj.M(chr), 1]);
lambda = 0.99;
obj.contrib(:) = lambda;

Pxz_y1 = zeros(obj.Mtot, obj.pop.Np);
Pxz_y0 = zeros(obj.Mtot, obj.pop.Np);

obj.chrInds();
obj.calcT(chr, obj.Alpha(1));
%%
totL = NaN(maxIter+1, 1);
totL(1) = 1;

for ii = 1 : maxIter
    mE1 = feval(obj.emissionHandle, double(obj.q), double(obj.r), obj.pop);
    sumE = sum(mE1, 2)./obj.pop.Np;
    mE0 = repmat( sumE, [1, obj.pop.Np]);
    
    % Ty = ones(2,2)./2;
    obj.E = bsxfun(@plus,...
        bsxfun(@times, mE1, obj.contrib), ...
        bsxfun(@times,mE0, (1-obj.contrib) ) ...
        );
    
    % obj.calcEmission;
    
    obj = obj.crossMatr(chr);
    obj = obj.cumMatr(chr);
    
    %= P[ x_ , y_m, z_m | lambda]
    Pxz_y1(obj.ci{chr}, :) = bsxfun(@plus, (obj.logAlpha{chr} + obj.logBeta{chr}), Pzm);
    Pxz_y0(obj.ci{chr}, :) = Pxz_y1(obj.ci{chr}, :) - log10(mE1(obj.ci{chr}, :)) + log10(mE0(obj.ci{chr}, :));
    
    Px_y1 = calcMarginal(Pxz_y1(obj.ci{chr}, 2:end),2);
    Px_y0 = calcMarginal(Pxz_y0(obj.ci{chr}, 2:end),2);
    
    % P[ z_m | y_m ]
    Pz_y1 = bsxfun(@minus, Pxz_y1(obj.ci{chr}, :), Px_y1) ;
    Pz_y0 = bsxfun(@minus, Pxz_y0(obj.ci{chr}, :), Px_y0) ;
    
%     obj.contrib(obj.ci{chr}) = 10.^( xPzy1 -  calcMarginal([xPzy0,xPzy1] , 2));    
    obj.contrib(obj.ci{chr}) = lambda.*10.^Px_y1./(lambda.*10.^Px_y1 + (1 - lambda).*10.^Px_y0 );
    
    lambdaNew = nanmean(obj.contrib(obj.ci{chr}));
    dLambda = lambdaNew - lambda;
    lambda = abs(lambdaNew);
        
    
     Pzm = calcMarginal( cat(3, bsxfun(@plus, Pz_y0, 1-obj.contrib(obj.ci{chr})) ,...
          bsxfun(@plus, Pz_y1, obj.contrib(obj.ci{chr}))  ), 3 );
      
%       calcMarginal(Pzm( :, 2:end),2)
      
    
    totL(ii+1) = calcMarginal([Px_y1; Px_y0]);
    if abs(totL(ii+1) - totL(ii)) < 1e-6 || all(isnan(Pzm(:)))
        fprintf('iteration:\t%u\tchange:\t%e \n', ii, (totL(ii+1) - totL(ii)))
        break
    end
end

fprintf('iteration:\t%u\tchange:\t%e \n', ii, (totL(ii+1) - totL(ii)))


figure; plot(totL)

indsTrue = obj.contrib > 0.5;

 figure
 surf(obj.x(obj.ci{chr}), obj.pop.kvect, Pzm', 'linestyle', 'none'); view( 0 , 90 )
 
 figure
 surf(obj.x(obj.ci{chr}), obj.pop.kvect, Pxz_y0(obj.ci{chr},:)', 'linestyle', 'none'); view( 0 , 90 )
 
 figure
 surf(obj.x(obj.ci{chr}), obj.pop.kvect, Pxz_y1(obj.ci{chr},:)', 'linestyle', 'none'); view( 0 , 90 )

figure; plot(obj.x(obj.ci{chr}), [Px_y0, Px_y1])
figure; plot(obj.x(obj.ci{chr}), Px_y0 - Px_y1)

figure; plot(obj.x(indsTrue& obj.ci{chr}), obj.f(indsTrue& obj.ci{chr}), 'gx')
hold all; plot(obj.x(~indsTrue& obj.ci{chr}), obj.f(~indsTrue& obj.ci{chr}), 'rx')

% obj.xPout(obj.ci{chr}) = calcMarginal(obj.xkPout(obj.ci{chr}, :), 2);
% if any(isnan(obj.xPout(obj.ci{chr})))
%     warning('runFBinternal:someEntriesAreNaN', 'some entries in the probability vector are NaN!')
% end
% if all(isnan(obj.xPout(obj.ci{chr})))
%     error('runFBinternal:allEntriesAreNaN', 'all entries in the probability vector are NaN!')
% end
