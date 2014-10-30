function [E, c] = wrapEmissionMatrix(q, r, study, emissionHandle, varargin)

%% process the input
contribInds = (cellfun(@numel, varargin) == numel(r));

%% initialize constants and compute Emission matrix and vectors
%= note that E is only normalized to the c: 
%=   c = sum(E, 1)
N = study.N;
[E, c] = feval(emissionHandle, double(q), double(r), study);

%= if the contribution of the true (simple Poisson) emission
%= to the mixture is given
%= flatten the emission distribution where the contribution of the 'true'
%= emisson is low:
if (nargin > 4) && any(contribInds)
    contrib = varargin{contribInds};
    Euniform = ones( 1, N+1)./(N+1);
    %= sanity check
    if any(contrib > 1)
        %         warning('propagateProbAlongChr:ENegInput','contrib > 1 in
        %         %u cases !\n', sum(contrib > 1) )
        fprintf('\t\t\tcontrib > 1 in %u cases !\n', sum(contrib > 1) )
        contrib(contrib > 1) = 1;
    end
    %
    %     E = bsxfun(@times, E, 1-(1-contrib).*(1-c)) + bsxfun(@times, Euniform, (1-contrib).*(1-c));
    
    E = bsxfun(@times, E, contrib) + bsxfun(@times, Euniform, (1-contrib.*c) );
    %         E = bsxfun(@times, E, contrib) + bsxfun(@times, Euniform, 1-contrib);
end

% if any( abs(sum(E,2)-1)> 1e-12 )
%     warning('propagateProbAlongChr:wrongEmissionVectors', ...
%         'emission probabilities do not sum to 1 !')
%     [col] = find(abs(sum(E,2)'-1)> 1e-12);
%     fprintf('Error \t@col(x ind) = \t%6u\t |r|=\t%1.3g\n',...
%         [col', sum(E(col,:),2)-1]' )
% end

%% sanity check
if   any(E(:)<0)
    warning('propagateProbAlongChr:ENegInput', 'E matrix has negative values!')
    [row, col] = find(E<0);
    fprintf('Negative value @row(k) = \t%u\t@col(x ind) = \t%u\t E =\t%1.3g\n',...
        row, col, E(row,col))
end

if   any(~isreal(E(:)))
    warning('propagateProbAlongChr:Imag', 'E matrix has complex values!')
    [row, col] = find(E<0);
    fprintf('Complex value @row(k) = \t%u\t@col(x ind) = \t%u\t E =\t%1.3g\n',...
        row, col, E(row,col))
end
