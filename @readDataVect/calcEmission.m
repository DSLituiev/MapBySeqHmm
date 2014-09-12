function obj = calcEmission(obj)

%% initialize constants and compute Emission matrix and vectors
%= note that E is only normalized to the c: 
%=   c = sum(E, 1)

if obj.flagWT 
    mE = feval(obj.emissionHandle, double(obj.q), double(obj.r), obj.pop);
    wE = feval(obj.emissionHandle, double(obj.qw), double(obj.rw), obj.pop);
    obj.E = mE.* fliplr(wE);
else
    obj.E  = feval(obj.emissionHandle, double(obj.q), double(obj.r), obj.pop);
end

if obj.preSelectionGen
    % likelihood that the sequence has not recombined out
    obj.Ey(:,1) = mE(:,1).* wE(:,1);
    % obj.Ey(:,2) = sum(mE,2).* sum(wE,2);
end

%= if the contribution of the true (simple Poisson) emission
%= to the mixture is given
%= flatten the emission distribution where the contribution of the 'true'
%= emisson is low:
if ~isempty(obj.contrib)
   
%   Euniform = ones( 1, obj.pop.Np)./(obj.pop.Np);
    Euniform = ones( 1, obj.pop.Np );
    %= sanity check
    if any(obj.contrib > 1)
        %         warning('propagateProbAlongChr:ENegInput','contrib > 1 in
        %         %u cases !\n', sum(contrib > 1) )
        fprintf('\t\t\tcontrib > 1 in %u cases !\n', sum(obj.contrib > 1) )
        obj.contrib(obj.contrib > 1) = 1;
    end
    %
    %     obj.E = bsxfun(@times, E, 1-(1-contrib).*(1-c)) + bsxfun(@times, Euniform, (1-contrib).*(1-c));
    
   sumE = sum(obj.E, 2);
    
   obj.E = bsxfun(@plus,...
           bsxfun(@times, obj.E, obj.contrib), ...
           bsxfun(@times, Euniform, sumE.*(1-obj.contrib) ) ...
           );
    %         obj.E = bsxfun(@times, obj.E, contrib) + bsxfun(@times, Euniform, 1-contrib);
end

% if any( abs(sum(E,2)-1)> 1e-12 )
%     warning('propagateProbAlongChr:wrongEmissionVectors', ...
%         'emission probabilities do not sum to 1 !')
%     [col] = find(abs(sum(E,2)'-1)> 1e-12);
%     fprintf('Error \t@col(x ind) = \t%6u\t |r|=\t%1.3g\n',...
%         [col', sum(E(col,:),2)-1]' )
% end

%% sanity check
if   any(obj.E(:)<0)
    warning('propagateProbAlongChr:ENegInput', 'E matrix has negative values!')
    [row, col] = find(E<0);
    fprintf('Negative value @row(k) = \t%u\t@col(x ind) = \t%u\t E =\t%1.3g\n',...
        row, col, E(row,col))
end

if   any(~isreal(obj.E(:)))
    warning('propagateProbAlongChr:Imag', 'E matrix has complex values!')
    [row, col] = find(E<0);
    fprintf('Complex value @row(k) = \t%u\t@col(x ind) = \t%u\t E =\t%1.3g\n',...
        row, col, obj.E(row,col))
end
