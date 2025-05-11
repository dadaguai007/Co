% DSP


Fc      = 193.1e12; 
% 色散恢复
% CD compensation
paramEDC    = struct();
paramEDC.L  = param.Ltotal;
paramEDC.D  = param.D;
paramEDC.Fc = Fc;
paramEDC.Fs = fs;


sigCDC = cdc(sig_RxMatch, paramEDC);
scatterplot(sigCDC(:,1))


%均衡即可
% Equ
paramEq=struct();
paramEq.nTaps=5;
paramEq.numIter = 5;
% paramEq.mu =5e-3;
paramEq.mu = [5e-3, 2e-4];
%RLS forgetting factor
paramEq.lambdaRLS=0.99;
paramEq.SpS=paramDec.SpS_out;
% coefficient taps
paramEq.H=[];
% Eq length
paramEq.L = [floor(0.2*length(d)), floor(0.8*length(d))];
% paramEq.L = [750,800];
paramEq.Hiter=[];
paramEq.storeCoeff='true';
% alg is the cell
% paramEq.alg='da-rde';
paramEq.constSymb=pnorm(GrayMapping(M, 'qam'));
paramEq.alg = {'cma','cma'};
[yEq, H, errSq, Hiter] = mimoAdaptEqualizer(x, d, paramEq);
scatterplot(yEq(:,1))
scatterplot(yEq(:,2))



% Constellation
z = (0:M-1)';
y = qammod(z,M);
% CPR frequence offest
paramCPR = struct();
paramCPR.alg = 'bps';
paramCPR.N= 85;
paramCPR.B= 64;
paramCPR.pilotInd = (1:20:length(SigRx)); 
paramCPR.Ts=Ts;
paramCPR.constSymb=y;
% paramCPR.constSymb= GrayMapping(M, 'qam') ;
[y_CPR_BPS,theta] = cpr(yEq,d,paramCPR);
y_CPR_BPS = pnorm(y_CPR_BPS);
eyediagram(y_CPR_BPS(1:10000,2),12)
title('After BPS CPR ')
scatterplot(y_CPR_BPS(:,1))
% scatterplot(y_CPR_BPS(:,2))
figure;hold on;
plot(real(y_CPR_BPS(:,1)))
plot(real(y_CPR_BPS(:,2)))





% % Equ
% paramEq=struct();
% paramEq.nTaps=5;
% paramEq.numIter = 1;
% % paramEq.mu =5e-3;
% paramEq.mu = [5e-3, 2e-4];
% %RLS forgetting factor
% paramEq.lambdaRLS=0.99;
% paramEq.SpS=paramDec.SpS_out;
% % coefficient taps
% paramEq.H=[];
% % Eq length
% paramEq.L = [floor(0.2*length(d)), floor(0.8*length(d))];
% % paramEq.L = [750,800];
% paramEq.Hiter=[];
% paramEq.storeCoeff='true';
% % alg is the cell
% % paramEq.alg='da-rde';
% paramEq.constSymb=pnorm(GrayMapping(M, 'qam'));
% paramEq.alg = {'da-rde','rde'};
% [yEq, H, errSq, Hiter] = mimoAdaptEqualizer(x, d, paramEq);
% scatterplot(yEq)

% 4th power frequency offset estimation/compensation
[Ei, ~] = fourthPowerFOE(SigRx, 1/Ts);
%norm
Ei = Ei./sqrt(mean(abs(Ei).^2));
scatterplot(Ei)
title('4th power frequency offset')


% Constellation
z = (0:M-1)';
y = qammod(z,M);
% CPR frequence offest
paramCPR = struct();
paramCPR.alg = 'bps';
paramCPR.N= 85;
paramCPR.B= 64;
paramCPR.pilotInd = (1:20:length(SigRx)); 
paramCPR.Ts=Ts;
paramCPR.constSymb=y;
% paramCPR.constSymb= GrayMapping(M, 'qam') ;
[y_CPR_BPS,theta] = cpr(x,d,paramCPR);
y_CPR_BPS = pnorm(y_CPR_BPS);
eyediagram(y_CPR_BPS(1:10000),24)
title('After BPS CPR ')
if 1
scatterplot(y_CPR_BPS)
title('After BPS CPR and downsample')
else
scatterplot(downsample(y_CPR_BPS,2))
title('After BPS CPR and downsample')
end
close all;

% CPR frequence offest
paramCPR.alg = 'ddpll';
paramCPR.tau1 = 1/(2*pi*10e3);
paramCPR.tau2 = 1/(2*pi*10e3);
paramCPR.Kv  = 0.1;
paramCPR.pilotInd = (1:20:length(SigRx)); 
paramCPR.Ts=Ts;
paramCPR.constSymb= GrayMapping(M, 'qam') ;
[y_CPR_PLL,theta1] = cpr(x,d,paramCPR);
y_CPR_PLL = pnorm(y_CPR_PLL);
eyediagram(y_CPR_PLL(1:10000),24)
title('After PLL CPR ')
if 1
scatterplot(y_CPR_BPS)
title('After PLL CPR and downsample')
else
scatterplot(downsample(y_CPR_BPS,2))
title('After PLL CPR and downsample')
end




% constellation = constellation / sqrt(bandpower(constellation));
% % [mma_out_mimo, error_mimo, w_mimo] = mimo2_2_mma(y, 21, [5e-3, 5e-3], 2^15, constellation, 2, 8);
% % mma_out_mimo = mma_out_mimo(:, 1) + 1j*mma_out_mimo(:, 2);
mf_sig = mf_sig / sqrt(bandpower(mf_sig));
[mma_out, error, w] = my_mma_lms(mf_sig, 21, 5e-3, 2^15, constellation, 2, 8);



% [vv_out, estimatedPhase] = V_V(mma_out, 10);
[fse_out, esfreq] = FSE(mma_out, Fb);
% [pll_out, esphase] = DDPLL(vv_out, 1/(2*pi*100e6), 1/(2*pi*100e6), 0.25, Fb, 0, constellation);
[pll_out, esphase_mimo] = DDPLL(fse_out, 1/(2*pi*100e6), 1/(2*pi*100e6), 0.2, Fb, 0, constellation);


[mma_out_mimo, error_mimo, w_mimo] = mimo2_2_mma([real(pll_out), imag(pll_out)], 21, [5e-3, 5e-3], 2^15, constellation, 1, 8);
[bps_out, estimatedPhase] = BPS(fse_out(1:2^15), constellation, 4, 2);
% [mma_out_mimo, error_mimo, w_mimo] = mimo2_2_mma([real(bps_out), imag(bps_out)], 21, [5e-3, 5e-3], 2^15, constellation, 1, 8);

