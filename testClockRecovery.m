% matlab中自带的时钟恢复函数

%% 初始化环境
clearvars;   % 清空工作区变量
clc;         % 清空命令窗口
close all;   % 关闭所有图形窗口
addpath('Phase_Sync\')

%% 调试配置
debug_tl_static  = 1; % 静态调试标志（1=显示最终星座图）

%% 系统参数配置
% 基本通信参数
L        = 32;         % 过采样倍数（每符号采样数）
M        = 16;          % 调制阶数（QAM/PAM）
N        = 2;          % 符号维度（1=PAM，2=QAM）
nSymbols = 1e5;        % 发送符号数量

% 定时恢复环路参数
Bn_Ts    = 0.01;       % 环路带宽×符号周期（归一化噪声带宽）
eta      = 1;          % 环路阻尼系数（控制收敛速度）

% 信号处理参数
rollOff  = 0.2;        % 升余弦滚降系数
timeOffset = 25;       % 信道延迟（采样点数）
fsOffsetPpm = 0;       % 采样时钟频偏（ppm单位）
rcDelay  = 10;         % 升余弦滤波器时延（符号数）
EsN0     = 20;         % 信噪比（Es/N0 dB）
Ex       = 1;          % 符号平均能量
TED      = 'MMTED';    % 定时误差检测器类型（关键算法选择）
intpl    = 2;          % 插值方法：0=多相，1=线性，2=二次，3=三次
forceZc  = 0;          % 强制零交符号模式（调试自噪声）

%% 创建系统对象
% 发射端升余弦滤波器
TXFILT = comm.RaisedCosineTransmitFilter( ...
    'OutputSamplesPerSymbol', L, ...  % 每符号输出L个采样
    'RolloffFactor', rollOff, ...     % 指定滚降系数
    'FilterSpanInSymbols', rcDelay);  % 滤波器跨度

% 接收端匹配滤波器（不进行抽取）
RXFILT = comm.RaisedCosineReceiveFilter( ...
    'InputSamplesPerSymbol', L, ...   % 输入采样率与发射端一致
    'DecimationFactor', 1, ...        % 保留原始采样率（后续处理需要）
    'RolloffFactor', rollOff, ...
    'FilterSpanInSymbols', rcDelay);
%%
% 延迟模块（模拟信道延迟）
DELAY = dsp.Delay(timeOffset);

% 生成参考星座（用于MER计算）
if (N==2)
    const = qammod(0:M-1,M);          % QAM星座图
else
    const = pammod(0:M-1,M);          % PAM星座图
end
Ksym = modnorm(const, 'avpow', Ex);   % 能量归一化因子
const = Ksym * const;                 % 调整后的参考星座

% 调制误差率(MER)测量模块
mer = comm.MER;
mer.ReferenceSignalSource = 'Estimated from reference constellation';
mer.ReferenceConstellation = const;   % 设置参考星座


%% 定时恢复环路常数计算
% 计算定时误差检测器增益
Kp = calcTedKp(TED, rollOff);         % 自定义函数计算TED增益

% 调整增益（考虑符号能量）
K  = 1;           % 假设信道增益为1（实际系统需AGC）
Kp = K * Ex * Kp; % 调整后的TED增益

%% MATLAB内置符号同步器（用于性能对比）
% TED类型映射（兼容MATLAB的命名）
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
    'SamplesPerSymbol', L, ...                   % 过采样倍数
    'NormalizedLoopBandwidth', Bn_Ts, ...        % 归一化环路带宽
    'DampingFactor', eta, ...                    % 阻尼系数
    'DetectorGain', Kp);                         % 检测器增益

%% 生成发射符号
if (forceZc)  % 强制零交符号模式（调试自噪声）
    data = zeros(nSymbols, 1);     % 创建全零向量
    data(1:2:end) = M-1;           % 交替设置最大值（产生零交）
else          % 正常随机数据模式
    data = randi([0 M-1], nSymbols, 1); % 生成随机整数
end

% 调制处理
if (N==2)
    modSig = Ksym * qammod(data, M);  % QAM调制
else
    modSig = real(Ksym * pammod(data, M)); % PAM调制（取实部）
end

%% 信号处理链路仿真
% 发射滤波
txSig = step(TXFILT, modSig); % 使用升余弦滤波器成形; 相当于conv，full 函数，去除拖尾的延迟


% 模拟采样时钟偏移（通过重采样）
fsRatio = 1 + (fsOffsetPpm * 1e-6); % 计算频率比例（Rx/Tx）
tol = 1e-9;                          % 重采样容差
[P, Q] = rat(fsRatio, tol);          % 转换为有理分数
txResamp = resample(txSig, P, Q);    % 执行重采样

% 信道延迟
delaySig = step(DELAY, txResamp);    % 添加固定延迟;直接添加零

% 添加高斯白噪声（AWGN）
txSigPower = 1 / sqrt(L);            % 计算信号功率
rxSeq = awgn(delaySig, EsN0, txSigPower); % 添加带限噪声

% 接收端匹配滤波
mfOut = step(RXFILT, rxSeq);         % 升余弦匹配滤波

%% 定时恢复处理
% 无定时恢复的直接抽取
rxNoSync = downsample(mfOut, L);     % 简单抽取（存在时偏）

% 理想定时恢复（已知延迟）
rxPerfectSync = downsample(mfOut, L, timeOffset); % 完美同步

% MATLAB内置定时恢复
rxSync2 = step(SYMSYNC, mfOut);      % 调用系统对象

%% 性能评估
skip = 0.2 * nSymbols; % 跳过初始瞬态（20%符号）

% 计算并显示MER结果
fprintf("\n测量MER结果:\n")
fprintf("无定时恢复: %.2f dB\n", mer(rxNoSync(skip:end)))
fprintf("理想定时恢复: %.2f dB\n", mer(rxPerfectSync(skip:end)))
fprintf("MATLAB %s算法: %.2f dB\n", matlabTed, mer(rxSync2(skip:end)))

%% 调试绘图（需要手动开启debug标志）
if (debug_tl_static)
     scatterplot(rxNoSync(skip:end));      % 未同步星座
    title('无定时恢复');

     scatterplot(rxPerfectSync(skip:end)); % 理想同步
    title('理想定时恢复');

     scatterplot(rxSync2(skip:end));       % MATLAB算法
    title(sprintf('MATLAB %s 算法', matlabTed));
end