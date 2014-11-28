clear all; clc
% close all; 

N = 5;
R = 4;
Q = (1:R)';

mbs = MapBySeqGeneric(N);
mbs.t = [0,1]*0.01;

for ii = flipud(Q)'
    for jj = flipud(Q)'
        mbs.r = R*[1,1]';
        mbs.q = [ii, jj]';
        mbs.initHMM();
        mbs.E;
        [xPout, xkPout] = mbs.HMM.getLikelihoodOfAModel(mbs.HMM.getHiddenStateModel('sel'));
        P_km(ii, jj) = xPout(1); 
    end
end

figure; pcolor(P_km)

p_q = mbs.HMM.getHiddenStateModel('sta') * feval(mbs.emissionHandle, Q, repmat(R,R,1))';

Expect_logP = p_q*P_km*p_q'
