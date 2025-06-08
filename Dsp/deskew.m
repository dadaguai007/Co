function sigRx = deskew(rIn, SpS, Rs, N, ParamSkew)
% I/Q支路时延偏差补偿（基于拉格朗日插值法）
%
% 输入参数:
%   rIn : 复数数组 - 输入信号（ADC采样后的复数信号）
%   SpS : 整数 - 每符号采样数 (Samples per Symbol)
%   Rs  : 整数 - 符号速率 (Baud rate, 单位: symbols/s)
%   N   : 整数 - 拉格朗日插值多项式阶数
%   ParamSkew : 结构体 - 时延参数:
%        - TauIV: I分量时延(秒)
%        - TauQV: Q分量时延(秒)
%
% 输出参数:
%   sigRx : 复数数组 - 时延补偿后的信号
%
% 参考文献:
%   [1] Digital Coherent Optical Systems, Architecture and Algorithms
%   [2] Tanimura et al., ECOC 2009 (简单的数字skew补偿器)

% 确保输入信号为列向量
rIn = rIn(:);

% 分离实部(I)和虚部(Q)
sigRx = [real(rIn), imag(rIn)];

% 计算ADC采样周期 T = 1/(采样率) = 1/(SpS*Rs)
TADC = 1 / (SpS * Rs);

% 计算以采样周期为单位的时延偏差
Skew = [ParamSkew.TauIV / TADC; ParamSkew.TauQV / TADC];

% 以最小时延为基准进行归一化 (使最小延迟=0)
Skew = Skew - min(Skew);

% 分解时延为整数部分和小数部分
nTADC = floor(Skew);      % 整数采样偏移
muTADC = -(Skew - nTADC); % 亚采样偏移（用于插值）

taps = N + 1;  % FIR滤波器抽头数 = 插值器阶数+1
buffer = [];   % 存储插值结果

% 创建插值索引 [-N/2, ..., 0, ..., N/2]
index = (0:N) - floor((N + 1)/2);

% 分别处理I支路(1)和Q支路(2)
for i = 1:2
    % 初始化拉格朗日系数数组
    L = zeros(1, taps);
    
    % 计算每个插值点的拉格朗日系数
    for indL = 1:taps
        n = index(indL) + nTADC(i);
        
        % 获取除当前点外的所有插值点
        m_points = index + nTADC(i);
        m_points(indL) = [];  % 移除当前点
        
        % 计算拉格朗日基函数 (核心公式)
        numerator = muTADC(i) - m_points;
        denominator = n - m_points;
        L(indL) = prod(numerator ./ denominator);
    end
    
    % 信号边界填充（防止卷积越界）
    pad_len = floor(taps/2);
    padded_sig = [zeros(pad_len, 1); sigRx(:, i); zeros(pad_len, 1)];
    
    % 创建卷积矩阵（等效于FIR滤波器）
    C = convmtx(padded_sig, taps);
    C = C(:, taps:end-1);  % 裁剪有效部分
    
    % 执行插值：系数行向量 × 卷积矩阵
    interp = L * C;
    
    buffer = [buffer; interp];
end

% 重组复数信号：I支路 + j*Q支路
sigRx = buffer(1, :).' + 1j * buffer(2, :).';
end
