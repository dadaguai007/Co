function [ polyFiltBank ] = polyDecomp(filt, L)
%% 多相分解：将FIR滤波器系数按相位拆分为L个子滤波器
% 输入:
%   filt - 原始滤波器系数向量（需为L带滤波器）
%   L    - 分解因子（子滤波器数量 = 插值因子）
% 输出:
%   polyFiltBank - L行矩阵，每行为一个子滤波器系数

%% 步骤1：补零对齐（确保滤波器长度为L的整数倍）
if (mod(length(filt), L) == 0)
    paddedFilt = filt; % 长度已对齐，无需补零
else
    % 计算需补零数量（使总长度 = L的整数倍）
    nZerosToPad = L - mod(length(filt), L); 
    % 在滤波器末尾补零
    paddedFilt = [filt zeros(1, nZerosToPad)]; 
end

%% 步骤2：多相分解（按列重组为L个子滤波器）
lenSubfilt = length(paddedFilt) / L; % 每个子滤波器长度
% 重组：将补零后的滤波器向量转换为L行矩阵，每行为一个子滤波器
polyFiltBank = reshape(paddedFilt, L, lenSubfilt); 

end
