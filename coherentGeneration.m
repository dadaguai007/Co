clc;clear;close all;

% 发射机参数
Tx=CoherentTx(     ...
                   fs, ...                  % 发射信号的采样率
                   fb, ...                  % 发射信号的波特率
                   2, ...                   % 随机信号的阶数
                   15, ...                  % prbs码的阶数
                   16, ...                  % 调制格式 M
                   4, ...                   % 每符号采样点
                   1e5, ...                 % 码元数目 
                   2, ...                   % 偏振状态
                   'sqrt', ...              % 脉冲形式    
                   0.2, ...                 % 滚降系数 
                   4096, ...                % 影响长度  
                   'nrz', ...               % 用户的成型滤波器  
                   'rand', ...              % 选择模式    
                   'system');               % 成型滤波器的生成方式

% 发射信号
[signalDualPol,qamDualPol]=Tx.dataOutput();
% 成型滤波器
hsqrt = Tx.systemHsqrt();

% phase noise
Pin=Tx.phaseNoise(sigTx,lw);
% LO
sigLO = Tx.laserMode(inputSignal,lw,RIN,Plo_dBm);




