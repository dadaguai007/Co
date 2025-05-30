function [out, estimatedPhase] = V_V(in, ML)
% Viterbi-Viterbi 相位估计算法 (适用于QPSK等正交调制)
% 输入参数:
%   in: 输入信号矩阵 [L×N], L=帧长, N=通道数
%   ML: 滑动窗口半宽 (窗长=2*ML+1)
%
% 输出参数:
%   out: 相位补偿后的信号 [L×N]
%   estimatedPhase: 估计的相位噪声轨迹 [L×N]

[L, N] = size(in);
% === 步骤1: 相位差提取与调制消除 ===
% 计算相邻样本的共轭乘积 (获得相对相位差)
d = in .* conj(circshift(in, 1));
% 四次方变换: 消除QPSK的90°相位跳变 (0°,90°,180°,270° → 0°)
d4 = d.^4;
% === 步骤2: 滑动窗口平均 ===
block = zeros(2*ML+1, N);
block(1+ML:end, :) = d4(1:1+ML, :);% 填充初始数据 (前ML+1个样本)
% 处理起始段 (前ML个样本)
for j = 1:ML
    S(j, :) = mean(block, 1); % 计算窗口均值
    % 更新窗口: 移除首行，添加新数据
    block = [block(2:end, :); d4(j+ML+1, :)];
end
% 处理中间段 (完整窗口)
for j = ML+1:L-ML
    % 取以j为中心的窗 (j-ML 到 j+ML)
    block = d4(j-ML:j+ML, :);
    S(j, :) = mean(block, 1);% 计算窗内均值
end
% 处理结束段 (后ML个样本)
for j = L-ML+1:L
    % 更新窗口: 移除首行，末尾补零
    block = [block(2:end, :); zeros(1, N)];
    S(j, :) = mean(block, 1); % 均值
end
% === 步骤3: 相位重建与补偿 ===
estimatedPhase = zeros(L+1, N);
% 通过相位差累积重建绝对相位:
%   angle(S)/4: 恢复原始相位差 (消除四次方变换影响)
%   cumsum: 从相位差累积得到绝对相位
estimatedPhase(2:end, :) = cumsum(angle(S)/4);
estimatedPhase = estimatedPhase(1:end-1, :);
out = in .* exp(-1j*estimatedPhase);
end