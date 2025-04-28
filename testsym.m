


clc;clear;close all;

addpath(genpath('Fncs'));
    
M=4;
Fb=32e9;


% addpath("device\");
addpath('lecroy\LeCroyDSO');

ipaddr = '192.168.1.17';
dso = LeCroyScope(ipaddr);
dso.ChannelConfig = [1, 2];  
dso.SamplingRate = 80e9; 
dso.Memory = 2e6;
data = dso.Capture; 

% datapath = 'Coherent_data_0421';
% if ~exist(datapath,'dir')
%     mkdir(datapath);
% end
% 
% prefix = '20250421_250mV_QPSK_8GBaud_20GSa_O_B2B';
% sT = 50;
% for idx = 1:50
%     filename = sprintf('%s\\%s-%d.mat',datapath, prefix, idx);
%     y = dso.readwaveform([1,3]);
%     save(filename,"y");
%     clear y;
% end

% y = dso.readwaveform([1,3]);
% % 
rx = resample(data.', 2*Fb, 80e9);
rx = rx(:, 1) + rx(:, 2)*1j;
% I = clk_recovery(y(:, 1).',2,4,Fb,Fs).';
% Q = clk_recovery(y(:, 2).',2,4,Fb,Fs).';
% scatterplot(y);
% y = normalize(y);
% rx=resample(y(:, 1)+1j*y(:, 2),Fb*2,Fs);
cdc_sig = cdc(rx, 40060, 16e-6, 0.08e3, 1550e-9, Fb*2);
ch = rcosdesign(0.5, 8, 2);
mf_sig = filter(ch, 1, cdc_sig);
% mf_sig = filter(ch, 1, rx);
scatterplot(cdc_sig);
constellation = qammod(0:M-1, M);
% constellation = constellation / sqrt(bandpower(constellation));
% % [mma_out_mimo, error_mimo, w_mimo] = mimo2_2_mma(y, 21, [5e-3, 5e-3], 2^15, constellation, 2, 8);
% % mma_out_mimo = mma_out_mimo(:, 1) + 1j*mma_out_mimo(:, 2);
mf_sig = mf_sig / sqrt(bandpower(mf_sig));
[mma_out, error, w] = my_mma_lms(mf_sig, 21, 5e-3, 2^15, constellation, 2, 8);
figure, plot(abs(error));
scatterplot(mma_out);
% [vv_out, estimatedPhase] = V_V(mma_out, 10);
[fse_out, esfreq] = FSE(mma_out, Fb);
% [pll_out, esphase] = DDPLL(vv_out, 1/(2*pi*100e6), 1/(2*pi*100e6), 0.25, Fb, 0, constellation);

[pll_out, esphase_mimo] = DDPLL(fse_out, 1/(2*pi*100e6), 1/(2*pi*100e6), 0.2, Fb, 0, constellation);
[mma_out_mimo, error_mimo, w_mimo] = mimo2_2_mma([real(pll_out), imag(pll_out)], 21, [5e-3, 5e-3], 2^15, constellation, 1, 8);
[bps_out, estimatedPhase] = BPS(fse_out(1:2^15), constellation, 4, 2);
% [mma_out_mimo, error_mimo, w_mimo] = mimo2_2_mma([real(bps_out), imag(bps_out)], 21, [5e-3, 5e-3], 2^15, constellation, 1, 8);
% scatterplot(bps_out);
% % scatterplot(mma_out_mimo);
% % getEVM(pll_out, constellation)
getEVM(mma_out_mimo(:, 1)+1j*mma_out_mimo(:, 2), constellation)
% % % scatterplot(downsample(wout, round(Fs/Fb)));