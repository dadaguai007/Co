% CO system Train
clc;close all;clear;
addpath('Plot\')
addpath('Dsp\')
% addpath('D:\BIT_PhD\Base_Code\Codebase_using\')
addpath('Sync\')
%% 系统初始化
% 信号生成
coherentGeneration;

% 时间轴和频率轴创建
[freq,time]   =  freq_time_set(Tx.TxPHY.NSym *Tx.TxPHY.sps,fs);

% 发射信号 双偏振
[signalDualPol,qamDualPol]  =  Tx.dataOutput();

% DSP参考信号
label=qamDualPol;
% 数据性能起始位置
skip = 0.2 * Tx.TxPHY.NSym; % 跳过初始瞬态（20%符号）
%% Tx and Rx Lo
% phase noise TX
lw   =  100e3;
RIN  =  0;
Plo_dBm  =  10;
phase_Noise  =  Tx.phaseNoiseMode(signalDualPol,lw);
% Rx LO
sigLO   = Tx.laserMode(signalDualPol,lw,RIN,Plo_dBm);
% 添加频偏
FO      = 150e6;                 % frequency offset
sigLO_Offset = Tx.addFrequencyOffset(sigLO,FO); % add frequency offset

%% 器件频响建立

% obj.Implementation.responType='Bessel';
% f3dB=20e9;
% order=5;
% verbose=0;
% [filt,H]=Rx.createFrequencyResponse(freq,order,f3dB,verbose);
% out1 = firFilter(filt.h, txSig); % 第一种滤波器工作方式,时域
% out2=filter_work(txSig,H); % 第二种滤波器工作方式，频域

%% 双偏发射机创建
Amp=0.5; % EA放大
for indMode = 1:Tx.TxPHY.Nmodes

    % X,Y偏振信号传输
    if indMode==1
        sigTxCh = iqm(phase_Noise, Amp*signalDualPol(:,indMode).',paramIQ);
    else
        sigTxCh = DP_iqm(phase_Noise, Amp*signalDualPol(:,indMode).',paramIQ);
    end
    % 每个通道的输出模式
    fprintf(' mode %d\t power: %.2f dBm\n', indMode, 10 * log10(((10 .^ (Pout_dBm(indMode) / 10)) * 1e-3/Tx.TxPHY.Nmodes) / 1e-3));
    % 设置功率
    sigTxCh=Tx.setSignalPower(sigTxCh,channelPowerType,Pout_dBm(inMode));
    power=signalpower(sigTxCh);
    fprintf(' mode %d optical signal power: %.2f dBm\n',indMode, 10 * log10(power / 1e-3));
    % 装载信号
    sigDualPolTx(:,indMode)=sigTxCh;

end

%% 添加信道延迟
% timeOffset = 25;       % 信道延迟（采样点数）
% % 延迟模块（模拟信道延迟）
% DELAY = dsp.Delay(timeOffset);
% % 信道延迟
% delaySig = step(DELAY, sigDualPolTx);    % 添加固定延迟;直接添加零

% % 延迟
% delay=[38*1e-12,35*1e-12];
% for indMode=1:Tx.TxPHY.Nmodes
%     delaySig(:,indMode)=Rx.addSkew(sigDualPolTx(:,indMode),delay(indMode));
% end

%% 信号传输
sigRxo=Rx.signalTran(sigDualPolTx,param);
if Tx.TxPHY.Nmodes==2
    fprintf('total ssfm signal power model 1: %.2f dBm\n', 10 * log10(signalpower(sigRxo(:,1))/ 1e-3));
    fprintf('total ssfm signal power model 2: %.2f dBm\n', 10 * log10(signalpower(sigRxo(:,2))/ 1e-3));
else
    fprintf('total ssfm signal power model 1: %.2f dBm\n', 10 * log10(signalpower(sigRxo(:,1))/ 1e-3));
end
% 频谱参考
mon_ESA(sigRxo,fs);

%% 相干接收
theta=0;
sig_RxE=Rx.coherentReceive(sigRxo,sigLO_Offset,theta,paramPD);


%% 应用集成的时钟偏移函数
% ppm=-0.5;
% jitter_rms=0; %1e-10
% txResamp1=Rx.addSamplingClockOffset(txSig,ppm,jitter_rms);



%% 匹配滤波
matchOut=Rx.matchFiltering(sig_RxE);


%% 定时恢复算法
% 参数初始化
Bn_Ts    = 0.01;       % 环路带宽×符号周期（归一化噪声带宽）
eta      = 1;          % 环路阻尼系数（控制收敛速度）
debug_tl_static  = 1; % 静态调试标志（1=显示最终星座图）
debug_tl_runtime = 0; % 运行时调试标志（1=显示同步过程示波器）

clockRecovery = DspSyncDecoding( ...
    fs,...         % 接收信号的采样率
    fb, ...        % 接收信号的波特率
    M,...         % 接收信号的格式
    fs/fb, ...     % 上采样率
    [],...         % 时钟信号的上采样率
    skip, ...      % 误码计算起始位置
    label,...      % 参考信号
    "MLTED");
clockRecovery.Implementation.eta       = eta;              % 环路阻尼因子（稳定性控制）
clockRecovery.Implementation.Bn_Ts     = Bn_Ts;           % 环路带宽 × 符号周期（控制同步速度）
clockRecovery.Implementation.Ex =1;
clockRecovery.Implementation.intpl=2;          % 插值方法：0=多相，1=线性，2=二次，3=三次
rollOff=Tx.Nr.psfRollOff;
rcDelay=Tx.Nr.psfLength ;
sps=Tx.TxPHY.sps;
% [ rxSync1 ] = symbolTimingSync(TED, intpl, sps, sig_RxE, matchOut, K1, K2, ...
%     const, Ksym, rollOff, rcDelay, debug_tl_static, debug_tl_runtime);
% 集成的算法
% %接收信号，为sig_RxE\txResamp1
% [rxSync,Kp] =clockRecovery.TED_recoverOptimalSamplingPoints(sig_RxE,matchOut,rollOff,rcDelay,const,Ksym);


% 搜寻时间信号的极值
% [rxSync,P,OptSampPhase,MaxCorrIndex]=clockRecovery.time_phase_Recovery(matchOut);

%% 接收机参数补齐
% 参考信号输入
Rx.Implementation.qam_signal=qamDualPol;
% 参考信号(解码使用)
Rx.createReferenceSignal()
% 创建参考星座图
[const,Ksym]=Rx.creatReferenceConstellation();

%% DSP 处理
% 色散补偿
% CD compensation
paramEDC = struct();
paramEDC.L = param.Ltotal;
paramEDC.D = param.D;
paramEDC.Fc = Fc;
outSignal=Rx.cdCompensation(matchOut,paramEDC);

scatterplot(outSignal(:,1))


%均衡即可
% Equ
paramEq=struct();
paramEq.nTaps=5;
paramEq.numIter = 5;
% paramEq.mu =5e-3;
paramEq.mu = [5e-3, 2e-4];
%RLS forgetting factor
paramEq.lambdaRLS=0.99;
paramEq.SpS=2;
% coefficient taps
paramEq.H=[];
% Eq length
paramEq.L = [floor(0.2*length(label)), floor(0.8*length(label))];
% paramEq.L = [750,800];
paramEq.Hiter=[];
paramEq.storeCoeff='true';
% alg is the cell
% paramEq.alg='da-rde';
paramEq.constSymb=const*Ksym;
paramEq.alg = {'cma','cma'};
[yEq, H, errSq, Hiter] = mimoAdaptEqualizer(outSignal, label, paramEq);
scatterplot(yEq(:,1))
scatterplot(yEq(:,2))




% 4th power frequency offset estimation/compensation
[FOE_out, ~] = fourthPowerFOE(yEq, 1/Ts);
%norm
FOE_out = pnorm(FOE_out);
scatterplot(FOE_out(:,1))
title('4th power frequency offset')


% CPR frequence offest
paramCPR = struct();
paramCPR.alg = 'bps';
paramCPR.N= 85;
paramCPR.B= 64;
paramCPR.pilotInd = (1:20:length(outSignal));
paramCPR.Ts=Ts;
paramCPR.constSymb=const*Ksym;
[y_CPR_BPS,theta_array] = cpr(yEq,label,paramCPR);
y_CPR_BPS = pnorm(y_CPR_BPS);
% eyediagram(y_CPR_BPS(1:10000,2),12)
% title('After BPS CPR ')
scatterplot(y_CPR_BPS(:,1))
scatterplot(y_CPR_BPS(:,2))


figure;hold on;
plot(real(y_CPR_BPS(:,1)))
plot(real(y_CPR_BPS(:,2)))



% % CPR frequence offest
% paramCPR.alg = 'ddpll';
% paramCPR.tau1 = 1/(2*pi*10e3);
% paramCPR.tau2 = 1/(2*pi*10e3);
% paramCPR.Kv  = 0.1;
% paramCPR.pilotInd = (1:20:length(outSignal));
% paramCPR.Ts=Ts;
% paramCPR.constSymb=const*Ksym;;
% [y_CPR_PLL,theta1] = cpr(yEq,label,paramCPR);




% cpr算法
% [vv_out, estimatedPhase] = V_V(yEq, 10);
% [pll_out, esphase] = DD_PhaseLockedLoop(vv_out, 1/(2*pi*100e6), 1/(2*pi*100e6), 0.25, fb, 0, const*Ksym);

% [fse_out, esfreq] = FSE(yEq, fb);
% [pll_out, esphase_mimo] = DD_PhaseLockedLoop(fse_out, 1/(2*pi*100e6), 1/(2*pi*100e6), 0.2, fb, 0, const*Ksym);
% [bps_out, estimatedPhase] = BlindPhaseSearch(fse_out, const*Ksym, 4, 2);


% 实部的分布
[value,percentReal] = amplitude_distribution(real(y_CPR_BPS(:,1)));
% 获取EVM
EVM = Rx.getEVM(y_CPR_BPS);
fprintf('receiver signal EVM : %.4f and %.4f\n', EVM(1),EVM(2));