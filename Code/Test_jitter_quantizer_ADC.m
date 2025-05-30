clc;clear;close all;
% Test_ADC
Fs = 3200;
fc = 100;

t = (0:299999) * (1/Fs);
pi_value = pi;

% generate sinusoidal signal
sig = sin(2 * pi_value * fc * t);

figure;hold on;
plot(t, sig, '-o', 'markersize', 4);
xlim([min(t), max(t)]);

% quantizer
% 量化比特
nBits = 2;
maxV=1;
minV=-1;
sig_q = quantizer(sig, nBits,maxV, minV);
plot(t, sig_q, '-k', 'markersize', 4);
xlim([0, 10 * 1 / fc]);



% 采样频率偏移，抖动
Fs_in = Fs;
Fs_out = 1.001*Fs;
jitter_rms = 1e-9;

% time jitter and mistake
t_dec = clockSamplingInterp(t, Fs_in, Fs_out, jitter_rms);
sig_dec = clockSamplingInterp(sig, Fs_in, Fs_out, jitter_rms);
figure;hold on;
plot(t_dec, sig_dec, '-o', 'markersize', 4,'Color','r');
plot(t, sig, '-o', 'markersize', 4,'Color','k');
xlim([0, 10 * 1 / fc]);
legend('采样频率偏移','正确采样')

%% Test ADC

paramADC = struct();
paramADC.Fs_in = Fs;
paramADC.Fs_out = Fs_out;
paramADC.jitter_rms = jitter_rms;
paramADC.nBits =  2;
paramADC.Vmax = max(real(sig));
paramADC.Vmin = min(real(sig));
paramADC.AAF = 'on';
paramADC.N = 1001;

sigRx = adc(sig, paramADC);
figure;hold on;
plot(t_dec, sigRx, '-o', 'markersize', 4,'Color','k');
plot(t, sig, '-o', 'markersize', 4);
xlim([0, 10 * 1 / fc]);