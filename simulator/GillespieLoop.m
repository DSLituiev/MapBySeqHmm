%%% Gillespie's Loop Function, 1D
% Version 3 : records the paths

function [xxSt tvector] =  GillespieLoop(Stmax,X0, N, MaxT)
%xout = X0;
%% generate random variables
mylogrand = log(1./rand(1, Stmax,'single'));
myrandRxn = rand(1, Stmax,'single');

%% initialize copy number and time values
xxSt = zeros(Stmax,2, 'int32');
tvector = zeros(Stmax,1,'single');
xxSt(1) =  X0; % initial value of X
T = single(0);
%% run the Gilespie's loop
for ii = 1:(Stmax-1)
     %% compute propensities and the waiting time
     % a = c'.*hmatr(:, xxSt(ii)+1); % c = 1; may be omitted;
     a = single(flipud(xxSt(ii,:))); 
     a0 = sum(a); % a0 = 1; to be omitted
     T = T + mylogrand(ii)/a0; % time counter
     tvector(ii+1) = T;        % put into time 
     %% determine which reaction occurs this time:
     a_roulette = cumsum(a)./a0;
     myl = logical(myrandRxn(ii)<a_roulette);
     j_ii = find(myl,1,'first');
     %% 'run' the reaction
    % xout = xout + N(2,j_ii); 
     xxSt(ii+1,:)= xxSt(ii,:)+ N(:, j_ii)';  
     %% break if the time exceeds 200 s
     if T > MaxT 
         % xout = xxSt(ii+1);
         % xxSt(ii+1:end,:) = repmat(xxSt(ii+1,:),[(size(xxSt,1) - ii ), 1]);   
         % tvector(ii+1:end) = T;
         xxSt = xxSt(1:ii,:);
         tvector = tvector(1:ii);
         break
     end
end

if T < MaxT 
    xxSt(ii+1:end,:) = NaN;
end

end