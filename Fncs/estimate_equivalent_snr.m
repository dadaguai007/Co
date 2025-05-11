function [SNR_mmse_le, SNR_mmse_dfe] = estimate_equivalent_snr(x, Nfft, ...
    OS_filt, OS, ro, varargin)
% 基于电频谱估计理想MMSE均衡后的等效信噪比
% 本函数通过分析接收信号的功率谱密度(PSD)，计算采用理想最小均方误差(MMSE)
% 线性均衡(LE)和判决反馈均衡(DFE)后的理论等效信噪比(SNR)
% 理论依据详见: https://doi.org/10.1002/0471439002.ch2
%
% 输入参数：
%   x        - 接收电信号矩阵，每列为一个通道，至少2倍过采样
%   Nfft     - 计算功率谱密度(PSD)使用的FFT点数
%   OS_filt  - 额外过采样因子(通常设为4以提高频谱分辨率)
%   OS       - 输入信号的过采样因子(例如2表示2倍过采样)
%   ro       - 升余弦滚降系数(范围0~1)
%   start_filt - (可选)光滤波器起始频点索引
%   end_filt   - (可选)光滤波器截止频点索引
%
% 输出参数：
%   SNR_mmse_le  - MMSE线性均衡后的等效SNR(dB)
%   SNR_mmse_dfe - MMSE-DFE均衡后的等效SNR(dB)，理论上可达信道容量

% 光纤通信系统性能评估示例
% Nfft = 4096;         % 高分辨率FFT
% received_signal = ... % 从光接收机获取的实际信号
% 
% % 估算理论性能
% [SNR_le, SNR_dfe] = estimate_equivalent_snr(...
%     received_signal, ...  % 接收信号
%     Nfft, ...            % 4096点FFT
%     4, ...               % OS_filt=4
%     2, ...               % OS=2(2倍过采样)
%     0.1 ...             % 滚降系数ro=0.1
% );


%% 估计信号和噪声功率谱密度(PSD) ----------------------------------------------
N_ch = size(x, 2);  % 获取信号通道数
pxx = pwelch(x, Nfft);  % 使用Welch方法估计功率谱密度，默认50%重叠汉宁窗

% 确定噪声估计频段范围(避开信号主瓣)
start_nois = ceil(length(pxx)/2/OS*(1+ro)) + round(Nfft/100);  % 噪声起始位置
end_nois = floor(-start_nois + length(pxx)) - round(Nfft/100); % 噪声结束位置

%% 去除光滤波器影响(可选参数处理) --------------------------------------------
if nargin > 5  % 如果提供了光滤波器参数
    start_filt = varargin{1};  % 获取光滤波器起始索引
    validateattributes(start_filt, {'numeric'},...
        {'scalar','integer','nonnegative'}, '', 'start_filt', 6);
    end_filt = varargin{2};    % 获取光滤波器结束索引
    validateattributes(end_filt, {'numeric'},...
        {'scalar','integer','>=',start_filt}, '', 'end_filt', 7);
else  % 未提供光滤波器参数时的默认处理
    start_filt = end_nois;    % 使用噪声结束位置作为默认值
    end_filt = end_nois;
end

%% 计算信号和噪声PSD ------------------------------------------------------
% 计算噪声平均功率谱密度(避开信号和滤波器频段)
Pnois = mean([pxx(start_nois:start_filt,:);...  % 前段噪声
    pxx(end_filt:end_nois,:)]);                % 后段噪声

% 计算信号功率谱密度(总PSD减去噪声基底)
Ps = max(eps, pxx - Pnois);  % 确保非负，保持数值稳定性

%% 上采样并叠加频域响应 ----------------------------------------------------
% 扩展频谱用于后续处理
Ps = [Ps(1:Nfft/2,:);...        % 保留正频率部分
    zeros(Nfft*(OS_filt-1),N_ch);...  % 中间插入零值实现上采样
    Ps(Nfft/2+1:end,:)];        % 保留负频率部分

% 初始化累加器
num = zeros(OS_filt*Nfft, N_ch);

% 定义频移量(覆盖±2个频点)
mu = (-2:2)';  

% 循环移位叠加频响(模拟信道与均衡器联合效应)
for i = 1:length(mu)
    num = num + circshift(Ps, [mu(i)*Nfft/OS, 0]);  % 循环移位并累加
end

%% 计算等效SNR ----------------------------------------------------------
% 计算每个频率点的信噪比
snr_per_freq = num ./ Pnois;  

% 截取有效频段(去除边缘效应)
snr_per_freq = [snr_per_freq(1:Nfft/4,:);...  % 正频率有效部分
    snr_per_freq(end-Nfft/4+1:end,:)];       % 负频率有效部分

% 计算MMSE线性均衡等效SNR(分贝)
SNR_mmse_le = -10*log10(2/Nfft*sum(1./(snr_per_freq+1))).'; 

% 计算MMSE-DFE等效SNR(基于信道容量公式转换)
SNR_mmse_dfe = 10*log10(exp(1))*2/Nfft*sum(log(snr_per_freq+1)).'; 
end
