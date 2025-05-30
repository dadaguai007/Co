function [out, esfreq] = FSE(in, Fs)
% 频率偏移估计与补偿算法 (Frequency Shift Estimator)
% 输入参数:
%   in: 输入信号矩阵 [L×N], L=帧长, N=通道数
%   Fs: 采样率 (Hz)
%
% 输出参数:
%   out: 频率补偿后的信号 [L×N]
%   esfreq: 估计的频率偏移量 [1×N] (单位: Hz)


L = length(in);
% maxL = round((maxfo/Fs*L)*4)+1;

% 计算信号相位的FFT幅度谱 (用于频率偏移估计)
judge = abs(fft(angle(in)));
for i = 1:size(in, 2)
     % === 频率偏移估计 ===
    % 寻找相位FFT幅度谱的最大值位置 (主峰)
    idx = find(judge(:, i) == max(judge(:, i)));
    % 排除直流分量干扰
    if idx(1) == 1
        judge(1) = 0;
        idx = find(judge(:, i) == max(judge(:, i)));
    end
    frac = (idx(1)-1)/4;
    esfreq(i) = frac * Fs / L;
    % === 频率偏移补偿 ===
    % 生成线性相位旋转向量 (模拟频率偏移)
    fshift = 2*pi*esfreq(i)/Fs*(0:L-1).';
    % 负向补偿频率偏移
    out(:, i) = in(:, i) .* exp(-1j*fshift);
    % === 补偿方向验证 ===
    % 计算补偿后信号的相位FFT
    dirc_judge = abs(fft(angle(out)));
    % 寻找新频谱的主峰位置
    dirc_idx = find(dirc_judge(:, i) == max(dirc_judge(:, i)));
    % 计算剩余频率偏移
    dirc_frac = (dirc_idx(1)-1)/4;
    dirc_esfreq = dirc_frac * Fs / L;
    % 若补偿后偏移量更大说明补偿方向错误
    if dirc_esfreq > esfreq(i)
         % 反转补偿方向 (正向补偿)
        out(:, i) = in(:, i) .* exp(1j*fshift);
    end
end

end