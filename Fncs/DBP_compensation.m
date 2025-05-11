function x = DBP_compensation(x, Parameters, varargin)
% 数字反向传播非线性补偿算法
% 通过逆向仿真光信号在光纤中的传输过程，补偿色散(CD)和非线性效应(SPM)
% 主要功能：
%   1. 频域色散补偿(反向传播模型)
%   2. 分步非线性相位补偿
%   3. 多跨段光纤传输逆向仿真
% 参考理论：非线性薛定谔方程逆向求解


% Lspan = 100e3;
% adB = 0.2e-3;
% gam = 1.3e-3;
% b2 = -21.27e-27;
% DBPsteps = 2;
% csi = 1;
% Pch = 0;
% Nspan = 5;
% DBP_compensation(x,Parameters,2,false);

%% 处理可选输入参数
if nargin > 2
    SPS = varargin{1};  % 获取用户指定的符号过采样率(Samples Per Symbol)
else
    SPS = 2;            % 默认2倍过采样(符合标准处理要求)
end
if nargin > 3
    cut = varargin{2};  % 获取是否截断输出的标志
else
    cut = true;         % 默认启用输出截断(去除边缘失真)
end

%% 从参数结构体提取光纤传输参数
fs = SPS * Parameters.Rs;       % 计算实际采样率 = 符号率 × 每符号样本数
b2 = Parameters.b2;             % 光纤色散系数(单位: s/(Hz·m))
Lspan = Parameters.Lspan;       % 单跨段光纤长度(单位: 米)
Pch = Parameters.Pch;           % 单信道发射功率(单位: dBm)
Nspan = Parameters.Nspan;       % 光纤跨段总数
adB = Parameters.adB;           % 光纤衰减系数(单位: dB/m)
csi = Parameters.csi;           % 非线性补偿强度调节因子(经验系数)  校准因子，用于调节非线性补偿强度
gam = Parameters.gam;           % 光纤非线性系数(单位: 1/(W·m))
Nstep = Parameters.DBPsteps;    % 每跨段分割的补偿步数

% 异常处理：如果跨段数<=0直接返回
if Nspan <= 0
    return;
end


%% 计算色散延迟并确定填充长度
L_orig = size(x, 1);                        % 原始信号长度(样本数)
max_d = 4*Nspan*abs(b2)*pi*(fs^2)*Lspan;    % 预估最大色散延迟(样本数)
L_pad = 2^nextpow2(L_orig + 1 + max_d);     % 扩展至最近的二次幂长度(优化FFT效率)

% 确定输出长度
if cut
    L_out = floor(L_orig - max_d);          % 截断输出长度(去除色散拖尾)
else
    L_out = L_orig;                         % 保留原始长度(可能含无效样本)
end

%% 对称填充信号(保持信号连续性)
% 生成填充内容(重复原始信号片段)
padding = repmat(x, ceil((L_pad-L_orig)/2/L_orig), 1);
% 构造填充后的信号矩阵
x_pad = [padding(1+size(padding,1)-floor((L_pad-L_orig)/2):end, :);...  % 前填充
    x;...                                                          % 原始信号
    padding(1:ceil((L_pad-L_orig)/2), :)];                         % 后填充

%% 生成频域坐标
f_step = fs / L_pad;                                % 频率分辨率(Hz)
f = fftshift((-fs/2:f_step:fs/2-f_step)).';         % 生成对称频率轴(-Fs/2 ~ Fs/2)

%% 计算分步补偿参数
att = adB / 20 / log10(exp(1));                     % 将衰减系数转换为Neper/m
delta = (1 - exp(-2*att*Lspan)) / Nstep;            % 等效跨段长度增量
% 计算最优分步长度(根据指数衰减模型)
L1 = log((1-((1:Nstep)'-1)*delta) ./ (1-(1:Nstep)'*delta)) / 2 / att;

%% 逆向跨段传播处理(核心算法)
for ii = 1:Nspan                                    % 遍历每个光纤跨段
    % 计算跨段末端功率(考虑衰减)
    P = 10^(Pch/10) * 1e-3 * exp(-2*att*Lspan);     % 转换为线性功率(W)

    for k = 1:Nstep                                 % 每跨段分步处理  每跨段分步计算，提升精度
        %% 计算当前步长参数
        L_step = L1(Nstep - k + 1);                 % 逆向选择步长(从后往前)
        Leff = (1 - exp(-2*att*L_step)) / 2 / att;  % 计算有效非线性作用长度

        %% 线性色散补偿(频域相位旋转)
        x_pad = ifft(fft(x_pad) .* ...
            exp(+1j*2*pi^2*b2*L_step*f.^2));        % +号表示逆向补偿色散

        %% 非线性相位补偿(时域相位旋转)
        P = P * exp(2*att*L_step);                  % 逆向功率递增(对应实际衰减)
        % 应用非线性相位旋转(SPM补偿)
        x_pad = x_pad .* exp(1j * csi * gam * Leff * P / SPS * sum(abs(x_pad).^2, 2));
    end
end

%% 移除填充并截断输出
x = double(x_pad(floor((L_pad-L_out)/2)+1 : ...     % 计算有效起始位置
    L_pad - ceil((L_pad-L_out)/2), :)); % 计算有效结束位置
end
