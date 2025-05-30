function [wout, estimatedPhase] = BlindPhaseSearch(win, constellation, ML, steptimes)
% 盲相位搜索算法 (Blind Phase Search)
% 输入参数:
%   win: 输入信号矩阵 [L×N], L=帧长, N=通道数
%   constellation: 调制星座点 (如QPSK: [1+1i, -1+1i, -1-1i, 1-1i]/sqrt(2))
%   ML: 滑动窗口半宽 (窗长=2*ML+1)
%   steptimes: 相位搜索步进倍数 (默认=2)
%
% 输出参数:
%   wout: 相位补偿后的信号 [L×N]
%   estimatedPhase: 最优相位估计值 [L×N]

if nargin < 4
    steptimes = 2;
end
win = win ./ sqrt(bandpower(win));
constellation = constellation / sqrt(bandpower(constellation));
[L, N] = size(win);
B = length(constellation) * steptimes;% 候选相位总数(相位分辨率)
% 生成候选相位偏移向量 (范围: -π/4 ~ π/4)
phaseshift = (-B/2:B/2-1)/B * pi / 2;
% 初始化输出矩阵
estimatedPhase = zeros(L, N);
wout = zeros(L, N);

for i = 1:N
    % === 相位旋转测试 ===
    % 生成候选相位集合: 每个样本旋转B个不同相位
    temp = win(:, i) * exp(1j*phaseshift);
    % 星座判决: 找到每个旋转后信号最近的星座点
    judge = decision(temp, constellation);
    % 计算平方误差:
    d = temp - judge;
    d2 = d .* conj(d);

    % === 滑动窗口误差累积 ===
    block = zeros(2*ML+1, B);
    block(1+ML:end, :) = d2(1:1+ML, :);
    % 1. 处理起始段 (前ML个样本)
    for j = 1:ML
        S(j, :) = sum(block, 1); %  当前窗口误差和
        % 更新窗口: 移除首行，添加新数据
        block = [block(2:end, :); d2(j+ML+1, :)];
    end
    % 2. 处理中间段 (完整窗口)
    for j = ML+1:L-ML
        block = d2(j-ML:j+ML, :);% 取以j为中心的窗
        S(j, :) = sum(block, 1);% 计算窗内误差和
    end
    % 3. 处理结束段 (后ML个样本)
    for j = L-ML+1:L
        % 更新窗口: 移除首行，末尾补零
        block = [block(2:end, :); zeros(1, B)];
        S(j, :) = sum(block, 1);% 误差和
    end
    % === 选择最优相位 ===
    % 找到每个样本点最小误差对应的相位索引
    [~, index] = min(transpose(S));
    % 获取最优相位值
    estimatedPhase(:, i) = phaseshift(index).';
    % 相位补偿输出
    wout(:, i) = win(:, i) .* exp(1j*estimatedPhase(:, i));
end

end