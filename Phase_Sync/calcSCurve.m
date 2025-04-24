function [ normTauE, g ] = calcSCurve(TED, rollOff, rcDelay)
%% 功能：通过解析公式计算定时误差检测器（TED）的S曲线（数据辅助模式）
% 输入:
%   TED      - TED类型: 'MLTED'/'ELTED'/'ZCTED'/'GTED'/'MMTED'
%   rollOff  - 升余弦滚降因子（0到1）
%   rcDelay  - 升余弦滤波器时延（符号数，默认10）
% 输出:
%   normTauE - 归一化定时误差向量（范围[-0.5, 0.5]，单位符号周期）
%   g        - S曲线值，表示TED输出随定时误差的变化

%% 参数校验与初始化
if nargin < 3
    rcDelay = 10; % 默认滤波器时延10符号
end

%% 全局参数设置
L = 1e3;    % 过采样因子（仅用于高精度生成S曲线，与实际系统无关）
K = 1;      % 假设信道增益归一化（需自动增益控制）
Ex = 1;     % 符号平均能量归一化

%% 生成升余弦自相关脉冲（等效Tx+Rx级联响应）
% r_p为升余弦滤波器的自相关函数（连续时间脉冲的离散采样）
% 参数说明：
%   1      - 符号速率（归一化）
%   L      - 每符号采样点数（仅用于生成高分辨率脉冲）
%   'normal' - 模式，生成非平方根升余弦脉冲
%   rollOff - 滚降因子
%   rcDelay - 滤波器时延（符号数）
r_p = rcosine(1, L, 'normal', rollOff, rcDelay);

%% 计算升余弦脉冲的导数（解析方法）
% 通过中心差分法近似导数：r_p'(t) ≈ [r_p(t+Δt) - r_p(t-Δt)] / (2Δt)
% 离散实现：使用滤波器系数h = L*[0.5, 0, -0.5]，Δt=1/L
h = L * [0.5 0 -0.5];            % 中心差分滤波器系数
r_p_diff = conv(h, r_p);          % 卷积计算导数脉冲
r_p_diff = r_p_diff(2:end-1);     % 截断首尾无效点（保持长度一致）

%% 定义定时误差范围与索引
tau_e = -L/2 : L/2;               % 定时误差（采样点单位，范围[-L/2, L/2]）
normTauE = tau_e / L;             % 归一化为符号周期单位（范围[-0.5, 0.5]）
% 计算中心索引（对应脉冲对称中心）
% rcDelay*L为脉冲中心位置，fliplr反转tau_e顺序以匹配导数方向
idx = L * rcDelay + 1 + fliplr(tau_e); 

%% 根据TED类型计算S曲线
switch (TED)
    case 'MLTED' % 最大似然TED（式8.30）
        % S曲线正比于导数脉冲值：g(τ_e) = K*Ex*r_p'(τ_e)
        g = K * Ex * r_p_diff(idx); 
        
    case {'ELTED', 'ZCTED'} % 早-迟/过零TED（式8.35和8.42）
        % S曲线为超前与滞后采样点差值：g(τ_e) = K*Ex*[r_p(τ_e+0.5T_s) - r_p(τ_e-0.5T_s)]
        % 注意：ELTED与ZCTED的S曲线解析式相同
        g = K * Ex * (r_p(idx + L/2) - r_p(idx - L/2));
        
    case 'GTED' % Gardner TED（式8.45）
        % 近似正弦特性：g(τ_e) = 4K²Ex*C*sin(2πτ_e/T_s)
        % 其中C = sin(πβ/2)/(4π(1-β²/4)), β=rollOff
        C = sin(pi * rollOff / 2) / (4 * pi * (1 - (rollOff^2 / 4)));
        g = (4 * K^2 * Ex) * C * sin(2 * pi * tau_e / L);
        
    case 'MMTED' % Mueller-Muller TED（式8.50）
        % 利用相邻符号间隔脉冲差值：g(τ_e) = K*Ex*[r_p(τ_e+T_s) - r_p(τ_e-T_s)]
        g = K * Ex * (r_p(idx + L) - r_p(idx - L));
end
end
