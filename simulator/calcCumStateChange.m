function cumdx = calcCumStateChange(dx_v, mut_v, tind)

%== multiply the state-derivative and mutation vectors:
DxMatr = bsxfun(@times, dx_v, mut_v);
%== sum up the states changes in the order they occur:
cumdx = cumsum( DxMatr(tind) );
