function inds = halve(inds, minBin, varargin)

M = fix(log2(numel(inds)/fix(minBin)));

if nargin>2
    N = varargin{1};
else
    N = 2^M;
end
    
thr = fix(numel(inds)/2);

if M <= 1 % thr < minBin
    inds(1 : thr) = N-1;
    inds(thr+1 : end) = N;
    fprintf('N:%u, M:%u, thr:%u\n', N, M, thr);
    return
end

 inds = [halve(inds(1 : thr), minBin, N - 2^(M-1)); halve(inds(thr+1 : end), minBin, N)];