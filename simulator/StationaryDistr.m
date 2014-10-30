function p = StationaryDistr(N)
% returns stationary distribution of the Ehrenfest model
p = zeros(N+1,1);
ii = 0;
p(ii+1) = 2^(-N).*nchoosek(N,ii);
% Disable the most recent warning
[~, msgid] = lastwarn;
if ~isempty(msgid)
    warning('off', msgid);
end

for ii = 1:N
    p(ii+1) = 2^(-N).*nchoosek(N,ii);
    [~, msgid] = lastwarn;
    if ~isempty(msgid)
        warning('off', msgid);
    end
end
% Enable all warnings
warning('on', 'all')
