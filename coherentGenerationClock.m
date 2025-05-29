clc;clear;close all;
addpath('Fncs\')
addpath('Plot\')
addpath('Dsp\')
% addpath('D:\BIT_PhD\Base_Code\Codebase_using\')
% addpath('D:\PhD\Project\Base_Code\Base\')

% 发射机参数
Tx=CoherentTx(     ...
    32e9, ...                % 发射信号的采样率
    1e9, ...                % 发射信号的波特率
    2, ...                   % 随机信号的阶数
    15, ...                  % prbs码的阶数
    16, ...                  % 调制格式 M
    32e9/1e9, ...           % 每符号采样点
    1e5, ...                 % 码元数目
    1, ...                   % 偏振状态
    'sqrt', ...              % 脉冲形式
    0.2, ...                 % 滚降系数
    10, ...                % 影响长度
    'nrz', ...               % 用户的成型滤波器
    'rand', ...              % 选择模式
    'system');               % 成型滤波器的生成方式

fs=Tx.TxPHY.fs;
fb=Tx.TxPHY.fb;
ta  = 1/fs;
% 成型滤波器
hsqrt = Tx.systemHsqrt();

% IQM 调制器参数
paramIQ.Vpi=2;
paramIQ.VbI=-paramIQ.Vpi;
paramIQ.VbQ=-paramIQ.Vpi;
paramIQ.Vphi=paramIQ.Vpi/2;
% PD 参数
paramPD=struct();
paramPD.B =fb;
paramPD.R =1;
paramPD.type = 'ideal';
paramPD.Fs=fs;
% ssfm
param=struct();
param.Ltotal = 40; %km
param.Lspan =10;
param.hz= 0.5;
param.alpha=0.2;
param.D = 16;
param.gamma = 1.3;
param.Fc = 193.1e12;
param.NF = 4.5;
param.amp='ideal';
param.Fs=fs;
param.maxIter = 10;      % maximum number of convergence iterations per step
param.tol = 1e-5;       % error tolerance per step
param.nlprMethod = 'True'; % use adaptive step-size based o maximum nonlinear phase-shift
param.maxNlinPhaseRot = 2e-2; % maximum nonlinear phase-shift per step

% 信号功率设置模式
channelPowerType='output';
% 信号功率
for inMode=1:Tx.TxPHY.Nmodes
    Pout_dBm(inMode)=0;
end

% 接收机参数
Rx=CoherentRx(     ...
    Tx, ...                % 发射机参数
    fs, ...                % 接收机采样率
    2, ...                 % Dsp的上采样因子
    hsqrt, ...             % 匹配滤波器
    [],...                 % 接收的光电流信号
    [],...                 % 参考信号
    'on');                 % 是否选取全部信号 或者 分段选取
