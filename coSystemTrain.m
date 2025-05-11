% CO system Train
clc;close all;clear;
addpath('Plot\')
addpath('Dsp\')
addpath('D:\BIT_PhD\Base_Code\Codebase_using\')

% 信号生成
coherentGeneration;

% 时间轴和频率轴创建
[freq,time]   =  freq_time_set(Tx.TxPHY.NSym *Tx.TxPHY.sps,fs);

% 发射信号 双偏振
[signalDualPol,qamDualPol]  =  Tx.dataOutput();
% phase noise TX
lw   =  100e3;
RIN  =  0;
Plo_dBm  =  10;
phase_Noise  =  Tx.phaseNoiseMode(signalDualPol,lw);
% Rx LO
sigLO   = Tx.laserMode(signalDualPol,lw,RIN,Plo_dBm);
% 添加频偏
FO      = 150e6;                 % frequency offset
sigLO_WithOffset = Tx.addFrequencyOffset(sigLO,FO); % add frequency offset

% 双偏发射机创建
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


% 信号传输
sigRxo=Rx.signalTran(sigDualPolTx,param);
if Tx.TxPHY.Nmodes==2
    fprintf('total ssfm signal power model 1: %.2f dBm\n', 10 * log10(signalpower(sigRxo(:,1))/ 1e-3));
    fprintf('total ssfm signal power model 2: %.2f dBm\n', 10 * log10(signalpower(sigRxo(:,2))/ 1e-3));
else
    fprintf('total ssfm signal power model 1: %.2f dBm\n', 10 * log10(signalpower(sigRxo(:,1))/ 1e-3));
end
% 频谱参考
mon_ESA(sigRxo,fs);

% 相干接收
theta=0;
sig_RxE=Rx.coherentReceive(sigRxo,sigLO_WithOffset,theta,paramPD);

% 匹配滤波
sig_RxMatch=Rx.matchFiltering(sig_RxE);

% 接收机参数补齐
% 参考信号输入
Rx.Implementation.qam_signal=qamDualPol;
% 参考信号(解码使用)
Rx.createReferenceSignal()
% 创建参考星座图
Rx.creatReferenceConstellation();
