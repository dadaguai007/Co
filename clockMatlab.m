% CO system Train
clc;close all;clear;
addpath('Plot\')
addpath('Dsp\')
addpath('Phase_Sync\')
addpath('Fncs\')

%% 信号生成
coherentGenerationClock;
% 时间轴和频率轴创建
[freq,time]   =  freq_time_set(Tx.TxPHY.NSym *Tx.TxPHY.sps,fs);
% 发射信号 
[txSig,QPSK]  =  Tx.dataOutput();


%% 定时恢复环路参数
Bn_Ts    = 0.01;       % 环路带宽×符号周期（归一化噪声带宽）
eta      = 1;          % 环路阻尼系数（控制收敛速度）

% 信号处理参数
timeOffset = 25;       % 信道延迟（采样点数）
fsOffsetPpm =0;       % 采样时钟频偏（ppm单位）
EsN0     = 20;         % 信噪比（Es/N0 dB）
Ex       = 1;          % 符号平均能量
TED      = 'ELTED';    % 定时误差检测器类型（关键算法选择） 'ELTED', 'ZCTED', 'GTED', 'MMTED' % 四种算法
intpl    = 2;          % 插值方法：0=多相，1=线性，2=二次，3=三次
forceZc  = 0;          % 强制零交符号模式（调试自噪声）
% 延迟模块（模拟信道延迟）
DELAY = dsp.Delay(timeOffset);

% 参考信号
M=Tx.TxPHY.M;
const = qammod(0:M-1,M );          % QAM星座图
% 调制误差率(MER)测量模块
mer = comm.MER;
mer.ReferenceSignalSource = 'Estimated from reference constellation';
mer.ReferenceConstellation = const;   % 设置参考星座

%% 定时恢复环路常数计算
% 计算定时误差检测器增益
Kp = calcTedKp(TED, Tx.Nr.psfRollOff);         % 自定义函数计算TED增益

% 调整增益（考虑符号能量）
K  = 1;           % 假设信道增益为1（实际系统需AGC）
Kp = K * Ex * Kp; % 调整后的TED增益

%% 算法创建
% TED类型映射
tedMap = containers.Map({'ELTED', 'ZCTED', 'GTED', 'MMTED'}, ...
    {'Early-Late (non-data-aided)', ...
    'Zero-Crossing (decision-directed)', ...
    'Gardner (non-data-aided)', ...
    'Mueller-Muller (decision-directed)'});
if strcmp(TED, "MLTED")  % MATLAB不支持MLTED时的兼容处理
    warning("MLTED not supported by MATLAB's synchronizer - using ZCTED");
    matlabTed = "ZCTED";
else
    matlabTed = TED;
end

% 创建符号同步对象
SYMSYNC = comm.SymbolSynchronizer(...
    'TimingErrorDetector', tedMap(matlabTed), ... % 选择TED类型
    'SamplesPerSymbol', Tx.TxPHY.sps, ...                   % 过采样倍数
    'NormalizedLoopBandwidth', Bn_Ts, ...        % 归一化环路带宽
    'DampingFactor', eta, ...                    % 阻尼系数
    'DetectorGain', Kp);                         % 检测器增益

%% Train
% 模拟采样时钟偏移（通过重采样）【可以进行优化，使用插值函数】
fsRatio = 1 + (fsOffsetPpm * 1e-6); % 计算频率比例（Rx/Tx）
tol = 1e-9;                          % 重采样容差
[P, Q] = rat(fsRatio, tol);          % 转换为有理分数
txResamp = resample(txSig, P, Q);    % 执行重采样

% 信道延迟
delaySig = step(DELAY, txResamp);    % 添加固定延迟;直接添加零

% 添加高斯白噪声（AWGN）
txSigPower = 1 / sqrt(Tx.TxPHY.sps);            % 计算信号功率
rxSeq = awgn(delaySig, EsN0, txSigPower); % 添加带限噪声

% match 
mfOut=Rx.matchFiltering(rxSeq);

%% 定时恢复处理
% 无定时恢复的直接抽取
rxNoSync = downsample(mfOut, Tx.TxPHY.sps);     % 简单抽取（存在时偏）

% 理想定时恢复（已知延迟）
rxPerfectSync = downsample(mfOut, Tx.TxPHY.sps, timeOffset); % 完美同步

% MATLAB内置定时恢复
rxSync2 = step(SYMSYNC, mfOut);      % 调用系统对象


%% 性能评估
skip = 0.2 * Tx.TxPHY.NSym; % 跳过初始瞬态（20%符号）

% 计算并显示MER结果
fprintf("\n测量MER结果:\n")
fprintf("无定时恢复: %.2f dB\n", mer(rxNoSync(skip:end)))
fprintf("理想定时恢复: %.2f dB\n", mer(rxPerfectSync(skip:end)))
fprintf("MATLAB %s算法: %.2f dB\n", matlabTed, mer(rxSync2(skip:end)))


%% 调试绘图（需要手动开启debug标志）
debug_tl_static  = 1; % 静态调试标志（1=显示最终星座图）
if debug_tl_static
     scatterplot(rxNoSync(skip:end));      % 未同步星座
    title('无定时恢复');

     scatterplot(rxPerfectSync(skip:end)); % 理想同步
    title('理想定时恢复');

     scatterplot(rxSync2(skip:end));       % MATLAB算法
    title(sprintf('MATLAB %s 算法', matlabTed));
end