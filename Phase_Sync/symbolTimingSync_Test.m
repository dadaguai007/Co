function [ xI ] = symbolTimingSync(TED, intpl, L, mfIn, mfOut, K1, K2, ...
    const, Ksym, rollOff, rcDelay, debug_s, debug_r)
% 符号定时恢复闭环控制系统
% ---------------------
% 输入参数:
% TED     - 定时误差检测器类型: 'MLTED'/'ELTED'/'ZCTED'/'GTED'/'MMTED'
% intpl   - 插值方法: 0)多相 1)线性 2)二次 3)三次
% L       - 过采样因子（每符号采样数）
% mfIn    - 匹配滤波器输入序列（多相插值时使用）
% mfOut   - 匹配滤波器输出序列（非多相插值时使用）
% K1/K2   - PI控制器比例/积分增益
% const   - 星座图（用于符号判决）
% Ksym    - 符号归一化因子（接收信号需除以该因子后判决）
% rollOff - 升余弦滚降因子（0~1）
% rcDelay - 升余弦滤波器时延（符号数，通常为Tx/Rx滤波器总时延）
% debug_s - 静态调试绘图开关（0关闭，1开启）
% debug_r - 实时调试示波器开关（0关闭，1开启）
% 输出:
% xI      - 同步后的符号序列（列向量）



%% 参数校验与默认值设置
if (nargin < 12)
    debug_s = 0; % 默认关闭静态调试绘图
end
if (nargin < 13)
    debug_r = 0; % 默认关闭实时调试示波器
end

%% 中点偏移计算（兼容奇偶过采样因子L）
midpointOffset = ceil(L / 2);       % 符号间中点偏移量（如L=4时为2）
muOffset = midpointOffset - L/2;   % μ偏移补偿（奇数L时为0.5）

%% 调制阶数提取与输入信号格式化
M = numel(const); % 根据星座图获取调制阶数（如16-QAM对应M=16）

% 强制输入信号转为列向量（确保后续处理一致性）
if (size(mfIn, 1) == 1)
    mfIn = mfIn(:);
end
if (size(mfOut, 1) == 1)
    mfOut = mfOut(:);
end

%% 选择输入信号源（多相插值用MF输入，其他用MF输出）
if (intpl == 0)
    inVec = mfIn;   % 多相插值：直接处理未滤波信号（联合完成MF）
else
    inVec = mfOut;  % 其他插值：需外部预先完成MF
end


%% 实时调试工具配置（星座图+定时误差监视）
if (debug_r)
    % 星座图显示器（显示解调符号分布）
    hScope = comm.ConstellationDiagram(...
        'ReferenceConstellation', Ksym * const, ... % 参考星座（归一化后）
        'XLimits', [-1.5 1.5]*max(real(const)), ... % X轴显示范围
        'YLimits', [-1.5 1.5]*max(imag(const)));     % Y轴显示范围

    % 分数间隔μ实时监视器（观察定时误差收敛过程）
    hTScopeCounter = dsp.TimeScope(...
        'Title', 'Fractional Interval μ', ...
        'YLimits', [-1 1], ...       % μ范围[-1,1]
        'TimeSpan', 1e4);            % 显示窗口长度
end


%% MLTED专用导数滤波器生成（非多相插值时）
if (intpl ~= 0 && strcmp(TED, 'MLTED'))
    % 设计根升余弦匹配滤波器
    mf = rcosdesign(rollOff, rcDelay, L);
    % 生成导数滤波器（中心差分法）
    dmf = derivativeMf(mf, L);
    % 对MF输入信号进行导数滤波（得到dMF输出）
    dMfOut = filter(dmf, 1, mfIn);
end

%% 多相插值滤波器组设计（intpl=0时启用）
if (intpl == 0)
    polyInterpFactor = 128; % 多相插值因子（细分128个相位）

    % 设计高倍过采样升余弦滤波器（联合实现插值与匹配滤波）
    interpMf = sqrt(polyInterpFactor) * ...
        rcosdesign(rollOff, rcDelay, L * polyInterpFactor);

    % 多相分解：将长滤波器拆分为polyInterpFactor个子滤波器
    polyMf = polyDecomp(interpMf, polyInterpFactor);
    polyMf = fliplr(polyMf); % 翻转滤波器系数（便于卷积计算）

    % 生成多相导数滤波器组（用于MLTED）
    polyDMf = zeros(size(polyMf));
    for i = 1:polyInterpFactor
        polyDMf(i, :) = derivativeMf(polyMf(i, :), L); % 逐行求导
    end
    polyDMf = fliplr(polyDMf);
else
    polyMf = []; % 非多相插值时置空
end

%% Farrow结构系数矩阵（二次/三次插值）
if (intpl == 2)
    % 二次插值系数（α=0.5优化参数）
    alpha = 0.5;
    b_mtx = flipud(fliplr(...
        [+alpha,      -alpha, 0; ...
         -alpha, (1 + alpha), 0; ...
         -alpha, (alpha - 1), 1; ...
         +alpha,      -alpha, 0]));
elseif (intpl == 3)
    % Table 8.4.2
    % 三次插值系数（标准Farrow系数）
    b_mtx = flipud(fliplr(...
        [+1/6,    0, -1/6, 0; ...
         -1/2, +1/2,   +1, 0; ...
         +1/2,   -1, -1/2, 1; ...
         -1/6, +1/2, -1/3, 0]));
else
    b_mtx = []; % 线性/多相插值无需系数
end

%% 定时恢复主循环初始化
nSamples = length(inVec); 
nSymbols = ceil(nSamples / L); 
xI = zeros(nSymbols, 1);   % 输出符号缓存
mu = zeros(nSymbols, 1);   % 分数间隔μ缓存
v  = zeros(nSamples, 1);   % PI控制器输出缓存
e  = zeros(nSamples, 1);   % 定时误差缓存

% 状态变量初始化
k = 0;          % 当前符号索引
strobe = 0;     % 插值触发标志（1表示需要插值）
cnt = 1;        % 模1计数器（初始值设为1以对齐首符号）
vi = 0;         % PI积分器状态
last_xI = inVec(1); % 上一个插值符号（初始化为第一个采样）

%% 确定循环边界（防止数组越界）
if strcmp(TED, 'ELTED')
    n_end = nSamples - L;    % ELTED需预留超前采样
elseif intpl > 1
    n_end = nSamples - 1;    % 二次/三次插值需预留两个采样
else
    n_end = nSamples;        % 其他情况处理全部采样
end

%% 确定循环起始点（多相插值时需考虑滤波器长度）
if (intpl == 0)
    poly_branch_len = size(polyMf, 2); % 多相子滤波器长度
    n_start = max(1, poly_branch_len - L); % 确保足够历史采样
    if (strcmp(TED, 'ELTED') || strcmp(TED, 'ZCTED') || strcmp(TED, 'GTED'))
        n_start = n_start + ceil(L/2); % 补偿零交叉点偏移
    end
else
    n_start = 1; % 其他插值方法从第一个采样开始
end

%% 主循环：逐个采样处理
for n = n_start:n_end 
    if strobe == 1
        % ================= 插值操作 =================
        xI(k) = interpolate(intpl, inVec, m_k, mu(k), b_mtx, polyMf);
        
        % ================= 定时误差检测 =================
        a_hat_k = Ksym * slice(xI(k)/Ksym, M); % 符号判决（去归一化后切片）
        switch (TED)
            case 'MLTED' % 最大似然TED（需导数插值）
                if intpl == 0
                    xdotI = interpolate(0, mfIn, m_k, mu(k), [], polyDMf);
                else
                    xdotI = interpolate(intpl, dMfOut, m_k, mu(k), b_mtx);
                end
                e(n) = real(a_hat_k)*real(xdotI) + imag(a_hat_k)*imag(xdotI);
                
            case 'ELTED' % 早-迟TED
                % 计算超前与滞后插值点差值
                early_idx = m_k + midpointOffset;
                late_idx = m_k - midpointOffset;
                x_early = interpolate(intpl, inVec, early_idx, mu(k)-muOffset, b_mtx, polyMf);
                x_late  = interpolate(intpl, inVec, late_idx, mu(k)+muOffset, b_mtx, polyMf);
                e(n) = real(a_hat_k)*(real(x_early)-real(x_late)) + ... 
                       imag(a_hat_k)*(imag(x_early)-imag(x_late));
                       
            case 'ZCTED' % 过零TED
                a_hat_prev = Ksym * slice(last_xI/Ksym, M); % 前一符号判决
                zc_idx = m_k - midpointOffset;
                x_zc = interpolate(intpl, inVec, zc_idx, mu(k)+muOffset, b_mtx, polyMf);
                e(n) = real(x_zc)*(real(a_hat_prev)-real(a_hat_k)) + ...
                       imag(x_zc)*(imag(a_hat_prev)-imag(a_hat_k));
                       
            case 'GTED' % Gardner TED（非数据辅助）
                zc_idx = m_k - midpointOffset;
                x_zc = interpolate(intpl, inVec, zc_idx, mu(k)+muOffset, b_mtx, polyMf);
                e(n) = real(x_zc)*(real(last_xI)-real(xI(k))) + ...
                       imag(x_zc)*(imag(last_xI)-imag(xI(k)));
                       
            case 'MMTED' % Mueller-Muller TED
                a_hat_prev = Ksym * slice(last_xI/Ksym, M); 
                e(n) = real(a_hat_prev)*real(xI(k)) - real(a_hat_k)*real(last_xI) + ...
                       imag(a_hat_prev)*imag(xI(k)) - imag(a_hat_k)*imag(last_xI);
        end
        last_xI = xI(k); % 更新历史插值符号
        
        % ================ 实时调试 ================
        if debug_r
            step(hScope, xI(k));       % 更新星座图
            step(hTScopeCounter, mu(k)); % 更新μ值曲线
        end
    else
        e(n) = 0; % 非插值时刻误差置零
    end
    
    % ================ PI控制器更新 ================
    vp = K1 * e(n);     % 比例项
    vi = vi + K2 * e(n);% 积分项（累积）
    v(n) = vp + vi;     % 控制器输出
    
    % ================ 模1计数器控制 ================
    W = 1/L + v(n);     % 动态步长（基础步长+控制量）
    strobe = cnt < W;   % 检查是否下溢（触发插值）
    if strobe
        k = k + 1;      % 更新符号索引
        m_k = n;        % 记录当前基带点索引
        mu(k) = cnt / W;% 计算分数间隔μ
    end
    cnt = mod(cnt - W, 1); % 更新计数器
end

% 裁剪输出序列（去除未完成的插值）
if strobe
    xI = xI(1:k-1); 
else
    xI = xI(1:k); 
end


%% 静态调试：绘制误差曲线、控制量、μ值
if (debug_s)
    figure; 
    plot(e); 
    title('定时误差 e(n)'); 
    xlabel('采样点 n'); 
    ylabel('误差值');
    
    figure; 
    plot(v); 
    title('PI控制器输出 v(n)'); 
    xlabel('采样点 n'); 
    ylabel('控制量');
    
    figure; 
    plot(mu, '.'); 
    title('分数间隔 μ(k)'); 
    xlabel('符号索引 k'); 
    ylabel('μ值');
end


end
