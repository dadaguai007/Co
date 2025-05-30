% Test DSP
clc;close all;clear;

%% 信号生成
coherentGenerationClock;
% 添加路径
addpath('Plot\')
addpath('Dsp\')
addpath('Phase_Sync\')
addpath('Fncs\')

%% 时间轴和频率轴创建
[freq,time]   =  freq_time_set(Tx.TxPHY.NSym *Tx.TxPHY.sps,fs);
% 发射信号 
[txSig,qamSig]  =  Tx.dataOutput();


% 实部的分布
[~,percentReal] = amplitude_distribution(real(txSig));

% 虚部的分布
[value,percentImag] = amplitude_distribution(imag(txSig));


% 滤波前
figure;
mon_ESA(txSig,fs);


% 带通滤波器
LowF=4e9;
HighF=5e9;
dataout=BPF(txSig,fs,LowF,HighF);
figure;
mon_ESA(dataout,fs);
% 带阻滤波器
dataout_BSF=BSF(txSig,fs,LowF,HighF,0);
figure;
mon_ESA(dataout_BSF,fs);
% 低通滤波器
bandwidth=5e9;
dataout_LPF=LPF(txSig,fs,bandwidth);
figure;
mon_ESA(dataout_LPF,fs);

% 高通滤波器
bandwidth=5e9;
dataout_HPF=HPF(txSig,fs,bandwidth);
figure;
mon_ESA(dataout_HPF,fs);


% 低通滤波器
N=101;
lowbandwidth=3e9;
h = lowpassFIR(lowbandwidth, fs, N, 'gauss');
y = firFilter(h, txSig); % 第一种滤波器工作方式

figure;
mon_ESA(y,fs);

%% 滤波器形状
bwl = Tx.TxPHY.fb*0.44;
type='bessel';
order=4;
sps=Tx.TxPHY.sps;
outputVoltage=0;
[wout,filtercoeffreq] = efilter(txSig, type, order, bwl, fs, sps, outputVoltage);
% 滤波器工作
out=filter_work(txSig,filtercoeffreq);
figure;
mon_ESA(wout,fs);

figure;
mon_ESA(out,fs);
%% 贝塞尔滤波器
f3dB=5e9;
order=5;
% 归一化的截止频率
fcnorm = f3dB/(fs/2);
verbose=0;
filt=Besslf_filter(order,fcnorm,verbose);
y1 = firFilter(filt.h, txSig); % 第一种滤波器工作方式,时域
% 
H = filt.H(freq/fs);            % 获取频响 (注意归一化)
H_test=filt.H_test(freq/fs);    % 带延时的频响

out1=filter_work(txSig,H); % 第二种滤波器工作方式，频域

figure;
mon_ESA(out1,fs);
figure;
mon_ESA(y1,fs);

%% 高斯
nlength=length(freq);
verbose=1;
filt_Gaussian=Gaussian_filter(order,fcnorm,nlength,verbose);

y2 = firFilter(filt_Gaussian.h, txSig); % 第一种滤波器工作方式,时域

figure;
mon_ESA(y2,fs);