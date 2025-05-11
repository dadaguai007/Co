function [y,phn] = phase_recovery_data_aided(x,a,Parameters)
%PHASE_RECOVERY_DATA_AIDED 数据辅助型相位恢复算法
% 本函数利用已知发送序列进行精准相位恢复，适用于离线分析或仿真验证场景
% 主要功能：
%   1. 基于功率互相关的信号对齐
%   2. 相位误差估计与滑动平均滤波
%   3. 多通道共轭补偿处理
% 注意事项：不能自动检测信号共轭问题，需手动配置参数

%% 参数提取 ===============================================================
Nmem = Parameters.cpe_memory;      % 相位记忆深度：滑动平均窗口长度（样本数）
sps = Parameters.cpe_sps;          % 输入信号过采样率（1或2）
is_conj = Parameters.cpe_conj;     % 共轭补偿标志（布尔向量，按通道设置）
Ntx = Parameters.Ntx;              % 发射序列数量（多输入场景）
ph = 1;                            % 下采样相位选择（适用于2SPS场景）

%% 信号预处理 =============================================================
x_ds = x(ph:sps:end,:); % 下采样至1SPS（选择指定相位点）
% 示例：当sps=2时，取奇数索引样本(1,3,5...)

%% 信号对齐处理 ===========================================================
% 通过功率互相关实现收发序列时域对齐
align_samples = size(a,1)+1:2*size(a,1); % 取中间段避免边界效应
a = power_tx_sequence_align(x_ds(align_samples,:),a); 

%% 序列扩展与共轭补偿 =====================================================
a = repmat(a,ceil(size(x_ds,1)/size(a,1)),1); % 循环扩展发送序列
a = a(1:size(x_ds,1),:); % 截断至接收信号长度

% 多通道共轭补偿（应对硬件相位反转）
for n = 1:Ntx 
    if is_conj(n) % 若当前通道需要共轭补偿
        % 间隔Ntx选取对应通道数据（支持多偏振场景）
        a(:,n:Ntx:end) = conj(a(:,n:Ntx:end)); 
    end
end

%% 相位误差估计 ===========================================================
% 计算瞬时相位差（接收信号与理想符号的相位偏差）
phase_diff = angle(x_ds) - angle(a); % 直接相位相减（弧度）
phn = exp(1j*phase_diff); % 转换为复数形式（保留幅度信息）

% 滑动平均滤波（抑制相位噪声）
phn = movmean(phn, Nmem); % 使用移动窗口平均

%% 上采样处理（适用于2SPS输入）===========================================
if sps == 2 
    % 相位误差矩阵重构：将1SPS结果扩展为2SPS
    phn = reshape(repmat(phn, 1, 2).', [], size(phn,2));
end

%% 相位补偿应用 ===========================================================
% 通过复数共轭乘法实现相位旋转补偿
% sign(phn)提取相位方向，conj反转相位偏差
y = x .* conj(sign(phn)); 
end

%% 辅助函数1：发射序列对齐 ================================================
function a = power_tx_sequence_align(x,a)
% 基于信号功率互相关的时域对齐方法
% 输入：
%   x - 接收信号片段（用于对齐参考）
%   a - 原始发送序列
% 输出：
%   a - 时域对齐后的发送序列

Lf = size(a,1);    % 训练序列长度
Nrx = size(x,2);   % 接收通道数
Ntx = size(a,2);   % 发射通道数
xc = NaN(Lf,Nrx*Ntx); % 预分配互相关矩阵

% 多通道互相关计算
for n = 1:Ntx 
    an = abs(a(:,n)) - mean(a(:,n)); % 发送序列功率（去均值）
    for m = 1:size(x,2)
        rm = abs(x(1:Lf,m)) - mean(x(1:Lf,m)); % 接收信号功率（去均值）
        xc(:,(n-1)*Ntx+m) = corrx(an,rm); % 计算循环互相关
    end
end

% 寻找最佳对齐位置
xc = abs(sum(xc,2)); % 多通道互相关结果叠加
[~,t0] = max(xc);    % 寻找互相关峰值
a = circshift(a,Lf-t0); % 循环移位对齐
end

%% 辅助函数2：快速互相关计算 ==============================================
function y = corrx(x,h)
% 基于FFT的快速互相关计算（替代MATLAB xcorr）
% 原理：y = IFFT(FFT(x) .* conj(FFT(h)))
N = max(size(x,1),size(h,1)); % 计算FFT长度
y = ifft(fft(x,N) .* fft(flip(conj(h),1),N)); % 频域相乘+逆变换
end
