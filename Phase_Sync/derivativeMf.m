function [dmf] = derivativeMf(mf, L)
%% 功能：基于中心差分法计算匹配滤波器的导数滤波器（dMF）
% 输入：
%   mf - 匹配滤波器系数向量
%   L  - 过采样因子（采样间隔 T = 1/L）
% 输出：
%   dmf - 导数匹配滤波器系数向量（与mf等长）

% 定义中心差分滤波器系数（基于数值导数公式）
% 公式来源：中心差分法 (f(x+Δx) - f(x-Δx)) / (2Δx)
% 其中 Δx = T = 1/L，因此分母为 2*(1/L) = 2/L
% 最终系数为 L/2 * [1, 0, -1]（等效于 L*[0.5, 0, -0.5]）
h = L * [0.5, 0, -0.5]; 

% 将中心差分滤波器与匹配滤波器进行卷积
% 结果长度 = length(mf) + length(h) - 1 = length(mf) + 2
central_diff_mf = conv(h, mf); 

% 截断首尾各一个元素，使导数滤波器长度与原匹配滤波器一致
% 最终长度 = length(mf) + 2 - 2 = length(mf)
dmf = central_diff_mf(2:end-1); 

end
