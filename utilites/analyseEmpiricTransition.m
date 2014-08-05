close all
chrNum = 5;

E = [];
tau = [];
jj = 1;

for cc = chrNum:-1:1    
    for ii = numel(taus{cc})-1:-1:1
        Eb(:,:,jj) = bsxfun(@times, Es{cc}(:,ii), Es{cc}(:,ii)');
        tau(jj,1) = taus{cc}(ii);
        jj = jj+1;
    end
    
end

figure
hist(tau, 40)