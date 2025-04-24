function [ Kp ] = calcTedKp(TED, rollOff, method, plotSCurve, rcDelay)
% 计算定时误差检测器（TED）的增益Kp（S曲线原点处斜率）
% 输入参数：
%   TED        - TED类型：'MLTED'/'ELTED'/'ZCTED'/'GTED'/'MMTED'
%   rollOff    - 升余弦滚降因子（0~1）
%   method     - 计算方法：'analytic'（解析公式）或'simulated'（仿真模拟），默认解析
%   plotSCurve - 是否绘制S曲线，默认false
%   rcDelay    - 升余弦滤波器时延（符号数），默认10
% 输出参数：
%   Kp         - TED增益，即S曲线在τ_e=0处的斜率（单位：电压/符号周期）

%% 参数默认值处理
if nargin < 3
    method = 'analytic'; % 默认使用解析方法计算S曲线
end

if nargin < 4
    plotSCurve = false;   % 默认不绘制S曲线
end

if nargin < 5
    rcDelay = 10;         % 默认滤波器时延10符号
end

%% 获取S曲线数据（两种方法）
if (strcmp(method, 'simulated'))
    % 仿真方法：通过发送测试信号生成S曲线
    [ normTauE, g ] = simSCurve(TED, rollOff, rcDelay);
else
    % 解析方法：直接通过数学公式计算S曲线
    [ normTauE, g ] = calcSCurve(TED, rollOff, rcDelay);
end

%% 计算TED增益Kp（S曲线原点处导数）
L = length(g) - 1;             % 过采样因子（根据S曲线点数推算）
midPoint = L/2 + 1;           % 中心点索引（对应τ_e=0）
assert(normTauE(midPoint) == 0); % 验证中心点是否为0

% 计算S曲线在τ_e=0附近的差分斜率
delta_y = g(midPoint + 1) - g(midPoint - 1);  % Y轴变化量（两侧点差值）
delta_x = normTauE(midPoint + 1) - normTauE(midPoint - 1); % X轴变化量（固定为2/L）
Kp = delta_y / delta_x;        % 斜率即为增益Kp（单位：电压/符号周期）

%% 可选：绘制S曲线（调试用）
if (plotSCurve)
    method(1) = upper(method(1)); % 首字母大写（美化标题）
    figure
    plot(normTauE, g, 'LineWidth', 1.5)
    xlabel("归一化定时误差 $\tau_e / T_s$", 'Interpreter', 'latex')
    ylabel("TED输出 $g(\tau_e)$", 'Interpreter', 'latex')
    title(sprintf("%s方法 - %s的S曲线（滚降系数=%.2f）", ...
        method, TED, rollOff), 'Interpreter', 'latex')
    grid on
    ax = gca;
    ax.FontSize = 12;          % 设置字体大小
    ax.GridLineStyle = '--';   % 网格线样式
end
end
