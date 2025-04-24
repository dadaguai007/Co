function [ K1, K2 ] = piLoopConstants(Kp, K0, eta, Bn_Ts, L)
%% 功能：计算符号定时环路PI控制器的比例和积分增益
% 输入：
%   Kp    - 定时误差检测器增益（TED输出斜率）
%   K0    - 插值控制器增益（控制相位调整速度）
%   eta   - 阻尼因子（默认1，确保临界阻尼）
%   Bn_Ts - 归一化环路噪声带宽（Bn*Ts，无量纲）
%   L     - 过采样因子（每符号L个采样）
% 输出：
%   K1, K2 - PI控制器的比例和积分增益

%% 计算流程

% 将符号周期归一化的噪声带宽转换为采样周期归一化：Bn*T = Bn*Ts / L
% 公式依据：T = Ts / L（T为采样间隔）
Bn_times_T = Bn_Ts / L; 

% 计算等效参数Theta_n（基于阻尼因子和噪声带宽）
% 公式来源：经典控制理论中的环路带宽与阻尼关系（参见文献 Eq. C.57）
Theta_n = (Bn_times_T) / (eta + (1 / (4 * eta))); 

% 通过连续系统传递函数类比计算中间变量（参见文献 Eq. C.56）
% Kp*K0*K1 和 Kp*K0*K2 的表达式对应标准二阶系统的系数
Kp_K0_K1 = (4 * eta * Theta_n) / (1 + 2*eta*Theta_n + Theta_n^2);
Kp_K0_K2 = (4 * Theta_n^2) / (1 + 2*eta*Theta_n + Theta_n^2);

% 计算最终的K1和K2（解耦TED和计数器增益的影响）
K1 = Kp_K0_K1 / (Kp * K0); % 比例增益
K2 = Kp_K0_K2 / (Kp * K0); % 积分增益

end
