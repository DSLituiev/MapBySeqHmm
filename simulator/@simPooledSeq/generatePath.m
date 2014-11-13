function varargout = generatePath(obj)

%= tv : vector of mutation locations
%= xv : vector of the number of mutant plants
%= t0 : location of the mutation under selection
%% === Generate the path recombination times:

k_tau = 1; % propensity

tmatr = 1/k_tau * cumsum( - log( rand(obj.Rec_max, obj.pop.N, 'single') ) );

%== find the cutoff time (the time when the last recombination in any path occured)
% [~, lastMinInd] = min(tmatr,[],2);
%
% T_max = tmatr( end , lastMinInd(end) );
%
% if T_max > obj.T_max
%     T_max = obj.T_max;
% end
%== collect the recombination event times into an ordered vector

[obj.tChr, tind] = sort(tmatr(:));

obj.tChr = obj.tChr(obj.tChr<= obj.T_max);
tind = tind(obj.tChr<= obj.T_max);
%% == generate the state-derivative matrix
%=  based on the order of the events:
dx_v = ones(obj.Rec_max,1);
dx_v(2:2:end) =  -1; %<= [1, -1, 1, -1, 1, -1...]'

%== pick up the position of the SNP under selection:
% t0 =  obj.xCausativeSNP;       %<= approximately


% indCausativeSNP = find(tvect > t0, 1,'first');
% t0 = tvect(indCausativeSNP);  %<= precisely

if obj.Selection %<= linked    
    [~, obj.indCausativeSNP] = min( abs(obj.tChr - obj.xCausativeSNP) );
    obj.t0 = obj.tChr(obj.indCausativeSNP);    
    %== select the points where at time t0 the state is 1
    %== find indices of the chromosome points coming first after 't0' in each plant:
    [~, indt0] = max(tmatr>obj.t0);   %<= abuse of the 'max' beheaviour
    % [-dx_v(indt0), 1-2*mod(indt0,2)']
    % mut_v = ( ones(1,obj.Nplants) - 2*mod(indt0,2) );
    mut_v = ( -dx_v(indt0)' );
    %== include the false positives
    mut_v( 1:(obj.FalsePositives) ) = - mut_v( 1:(obj.FalsePositives) );
    cumdx = calcCumStateChange(dx_v, mut_v, tind);
    %== offset cumdx(indCausativeSNP)
    obj.kChr = obj.pop.N - obj.FalsePositives - cumdx(obj.indCausativeSNP) + cumdx;
else   %<= stationary
    obj.indCausativeSNP = 1;
    obj.t0 = 0;   
    %== assign (almost each second plant) '-1'
    Nmut = ceil(obj.pop.N./2);
    % Nmut = random('binom',obj.Nplants, 0.5);
    mut_v = ones(1, obj.pop.N);
    mut_v(1, 1:Nmut ) = -1; %<= [1, -1, 1, -1, 1, -1...]
    
    cumdx = calcCumStateChange(dx_v, mut_v, tind);
    obj.kChr = obj.pop.N - Nmut + cumdx;
end

% % %% approximate scaling
% % tvect = tvect(tvect<T_max);
% % x = x(tvect<T_max);
% % % tvect =   tvect* length(x)/ ( T_max ) ;
% % % t0 = t0 * length(x)/ ( T_max ) ;
% % T_max = max(tvect);

varargout = {obj.tChr, obj.kChr, obj.t0, obj.indCausativeSNP, obj.T_max};