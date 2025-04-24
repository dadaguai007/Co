function [x, y] = genTestVector(nSymbols, sps, rolloff, rcDelay, ...
    tedChoice, interpChoice, loopBw, dampingFactor)
% GEN_TEST_VECTOR 符号定时恢复测试向量生成器
% 功能：生成带定时误差的接收信号及恢复后的符号序列，用于验证定时同步模块
% ---------------------------------------------------------------------
% 输入参数:
%   nSymbols      - 生成符号数 (如1000个QPSK符号)
%   sps           - 过采样率 (如4表示每符号4个采样点)
%   rolloff       - 升余弦滚降因子 (0~1，典型值0.2~0.5)
%   rcDelay       - 升余弦滤波器时延 (符号数，决定滤波器长度)
%   tedChoice     - 定时误差检测器类型 ('GTED','ZCTED','MMTED'等)
%   interpChoice  - 插值器类型 (0=多相,1=线性,2=二次,3=三次)
%   loopBw        - 环路带宽 (归一化值，控制收敛速度)
%   dampingFactor - 阻尼系数 (典型值1为临界阻尼)
% 输出:
%   x - 接收端匹配滤波器输入/输出信号 (测试输入向量)
%   y - 符号同步后输出序列 (期望输出向量)

%% 星座图初始化 (QPSK)
M = 4;          % 调制阶数 (QPSK)
Ex = 1;         % 符号平均能量
% 生成标准QPSK星座并做能量归一化
const = qammod(0:(M-1), M);           % 生成初始星座
Ksym = modnorm(const, 'avpow', Ex);   % 计算归一化系数
const = Ksym * const;                 % 能量归一化星座

%% 定时恢复环路参数计算
% 计算定时误差检测器增益Kp (S曲线斜率)
Kp = calcTedKp(tedChoice, rolloff);
K0 = -1;  % 环路增益符号补偿 (通常固定为-1)
% 计算PI控制器参数 (比例项K1/积分项K2)
[K1, K2] = piLoopConstants(Kp, K0, dampingFactor, loopBw, sps);

%% 根升余弦滤波器设计 (Tx/Rx匹配滤波)
% Tx成型滤波器 (脉冲整形)
htx = rcosdesign(rolloff, rcDelay, sps); 
% Rx匹配滤波器 (与Tx共轭对称)
hrx = htx;  % 根升余弦特性确保联合响应为升余弦

%% 测试信号生成流程
% 生成随机QPSK符号序列
data = randi([0 M-1], nSymbols, 1);        % 随机整数生成
test_syms = Ksym * qammod(data, M);        % QPSK调制 + 能量归一化

% 发射端处理链
test_syms_up = upsample(test_syms, sps);   % 过采样 (插入sps-1个零)
x_mf = conv(test_syms_up, htx, 'same');    % 成型滤波 (保留中间有效部分)

% 接收端处理链 (理想信道，无噪声)
y_mf = conv(x_mf, hrx, 'same');            % 匹配滤波输出

%% 符号定时恢复处理
% 调用符号同步核心算法 (闭环恢复)
y_sync = symbolTimingSync(tedChoice, interpChoice, sps, x_mf, y_mf, ...
    K1, K2, const, Ksym, rolloff, rcDelay);

%% 测试向量格式化输出
% 根据插值器类型选择输入信号源:
% - 多相插值直接处理MF输入，其他插值处理MF输出
if (interpChoice == 0)
    printComplexVec(x_mf, "in");   % 输出MF输入序列
else
    printComplexVec(y_mf, "in");   % 输出MF输出序列
end
printComplexVec(y_sync, "out");    % 输出同步后符号

% 返回变量 (兼容函数调用)
x = x_mf;  
y = y_sync;
end

%% 辅助函数：复数向量格式化输出
function [] = printComplexVec(x, label)
% 功能：将复数向量格式化为 (a+bj), (c+dj)... 形式
% 用于生成其他语言(如Python/C)可直接解析的测试数据
    fprintf("%s = [", label);
    % 遍历前N-1个元素，添加逗号分隔
    for i = 1:(length(x) - 1)
        fprintf('(%f%+fj), ', real(x(i)), imag(x(i)));
    end
    % 最后一个元素单独处理 (避免末尾逗号)
    i = length(x);
    fprintf('(%f%+fj)]\n', real(x(i)), imag(x(i)));
end
