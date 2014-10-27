
ell = 1;
em = 1;
Qy = [-ell, ell; em, -em];

t = 0:0.1:2.4;

for ii= numel(t):-1:1
    T(:,:, ii) = expm( Qy.* t(ii) );
end

pt = squeeze(mtimesx([0,1], T));

figure
plot(t, pt(1,:))
hold all

log2(pt(1,t == .5))