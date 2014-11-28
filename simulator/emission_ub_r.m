function emission_r = emission_ub_r(~, r, pop)
%EMISSION_UB_R Summary of this function goes here
%   Detailed explanation goes here

assert(~isempty(r), 'provide read numbers per locus!')

q = (1:max(r))';
r = permute(r, [3,2,1] );

r_q = bsxfun(@minus, r, q);
r_choose_q = exp(  bsxfun(@minus,gammaln(r),gammaln(q)) - gammaln(r_q) );

Pq_z_r = bsxfun(@times, r_choose_q, ...
    bsxfun(@times,...
    bsxfun(@power, pop.f, q),...
    bsxfun(@power, ( 1 - pop.f), r_q )) );

Pq_stat = mtimesx( pop.Pstat,  Pq_z_r, 'T');

emission_r =  squeeze(mtimesx( Pq_stat, Pq_z_r ))';

end

