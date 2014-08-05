function varargout = generatePath(par)


%% === Generate the path recombination times:

k_tau = 1; % propensity

tmatr = 1/k_tau * cumsum( - log( rand(par.Rec_max, par.Nplants, 'single') ) );

%== find the cutoff time (the time when the last recombination in any path occured)
% [~, lastMinInd] = min(tmatr,[],2);
%
% T_max = tmatr( end , lastMinInd(end) );
%
% if T_max > par.T_max
%     T_max = par.T_max;
% end
%== collect the recombination event times into an ordered vector

[tvect, tind] = sort(tmatr(:));

tvect = tvect(tvect<= par.T_max);
tind = tind(tvect<= par.T_max);
%% == generate the state-derivative matrix
%=  based on the order of the events:
dx_v = ones(par.Rec_max,1);
dx_v(2:2:end) =  -1; %<= [1, -1, 1, -1, 1, -1...]'

%== pick up the position of the SNP under selection:
% t0 =  par.cSNPpos;       %<= approximately


% t0ind = find(tvect > t0, 1,'first');
% t0 = tvect(t0ind);  %<= precisely

if par.Selection %<= linked    
    [~, t0ind] = min( abs(tvect - par.cSNPpos) );
    t0 = tvect(t0ind);    
    %== select the points where at time t0 the state is 1
    %== find indices of the chromosome points coming first after 't0' in each plant:
    [~, indt0] = max(tmatr>t0);   %<= abuse of the 'max' beheaviour
    % [-dx_v(indt0), 1-2*mod(indt0,2)']
    % mut_v = ( ones(1,par.Nplants) - 2*mod(indt0,2) );
    mut_v = ( -dx_v(indt0)' );
    %== include the false positives
    mut_v( 1:(par.FalsePositives) ) = - mut_v( 1:(par.FalsePositives) );
else   %<= stationary
    t0ind = 1;
    t0 = 0;
    mut_v = ones(1,par.Nplants);
    %== assign (almost each second plant) '-1'
    Nmut = ceil(par.Nplants./2);
    % Nmut = random('binom',par.Nplants, 0.5);
    mut_v(1, 1:Nmut ) = -1; %<= [1, -1, 1, -1, 1, -1...]
end

%== multiply the state-derivative and mutation vectors:
DxMatr = bsxfun(@times, dx_v, mut_v);
%== sum up the states changes in the order they occur:
cumdx = cumsum( DxMatr(tind) );

if par.Selection  %<= linked
    %== offset cumdx(t0ind)
    x = par.Nplants - par.FalsePositives - cumdx(t0ind) + cumdx;
    % varargout = {tvect(tvect<T_max), x(tvect<T_max), t0};
else    %<= stationary
    x = par.Nplants - Nmut + cumdx;
    % varargout = {tvect(tvect<T_max), x(tvect<T_max)};
end

% % %% approximate scaling
% % tvect = tvect(tvect<T_max);
% % x = x(tvect<T_max);
% % % tvect =   tvect* length(x)/ ( T_max ) ;
% % % t0 = t0 * length(x)/ ( T_max ) ;
% % T_max = max(tvect);

varargout = {tvect, x, t0, t0ind, par.T_max};