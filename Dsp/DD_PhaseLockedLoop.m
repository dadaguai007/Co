function [wout, esphase] = DD_PhaseLockedLoop(win, tau1, tau2, k, Fs, phase, constellation)
% 判决导向数字锁相环（Decision-Directed Phase Locked Loop），用于载波相位同步。
% 核心功能是校正输入信号的相位偏差，使其与标准星座图对齐。
% 输入参数:
%   win: 输入信号矩阵 [L×N], L=帧长, N=通道数
%   tau1, tau2: 环路滤波器时间常数
%   k: 环路增益 (控制收敛速度)
%   Fs: 采样率 (Hz)
%   phase: 初始相位估计 (弧度)
%   constellation: 调制星座点 
%
% 输出参数:
%   wout: 相位校正后的输出信号 [L×N]
%   esphase: 相位估计轨迹 [L×N]

% tau1 主导积分分量（消除稳态误差）
% tau2 主导比例分量（提高响应速度）

[L, N] = size(win);
win = win ./ sqrt(bandpower(win));
% 基于时间常数 tau1（积分时间）、tau2（比例时间）和采样率 fs 设计二阶比例积分滤波器
filtercoeff = [1, 1/(2*tau1*Fs)*(1+[-1, 1]/tan(1/(2*tau2*Fs)))];
% 初始化相位估计矩阵 (首行为初始相位)
esphase = phase * ones(L+1, N);
wout = zeros(L, N);
for i = 1:N
    phaseout = 0; % 当前相位误差
    filterout = 0; % 滤波器输出缓存
    % 逐样本处理信号
    for idx = 1:L
        prephaseout = phaseout;% 保存上一时刻相位误差
        % 1. 相位旋转补偿
        wout(idx, i) = win(idx, i) * exp(-1j*esphase(idx, i));
        % 2. 星座判决: 找到最近星座点
        [~, index] = min(abs(constellation-wout(idx, i)));
        % 3. 计算相位误差 (当前信号与理想符号的夹角)
        phaseout = angle(wout(idx, i)*conj(constellation(index)));
        % 4. 更新环路滤波器 (二阶比例积分)
        filterout = filtercoeff * [filterout; prephaseout; phaseout];
         % 5. 更新下一时刻相位估计
        esphase(idx+1, i) = esphase(idx, i) + k*filterout;
    end
end
esphase = esphase(1:end-1, :);
end