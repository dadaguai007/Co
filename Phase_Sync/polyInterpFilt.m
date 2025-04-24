function [ polyFiltBank ] = polyInterpFilt(I, P, Blim)
% POLYINTERPFILT 多相插值滤波器组生成器
%   该函数设计并分解Lth带抗镜像滤波器，生成高效多相插值滤波器组
%
% 输入参数:
%   I    - 插值因子 (上采样倍数，如2表示2倍插值)
%   P    - 邻域样本数 (单侧样本数，总样本数2P，控制滤波器长度)
%   Blim - 带宽限制因子 (0~1，1=全带宽Nyquist，0.5=半带宽)
%
% 输出:
%   polyFiltBank - 多相滤波器组矩阵，I行×N列，每行为一个子滤波器
%
% 示例: 
%   polyFiltBank = polyInterpFilt(2, 2, 0.5)
%   设计2倍插值滤波器，使用4个邻域样本(P=2)，信号带宽为Nyquist的50%

%% 抗镜像滤波器设计 (Lth带滤波器)
% 设计参数:
%   I    - 零交叉间隔=插值因子，保证插值时原样点保留
%   P    - 滤波器长度=2*I*P-1 (例：I=2,P=2 → 长度7)
%   Blim - 控制过渡带宽，值越小抗混叠能力越强但通带纹波可能增加
interpFilt = intfilt(I, P, Blim); 

%% 多相分解
% 将长滤波器分解为I个并行的子滤波器(多相分支)
% 分解原理:
%   原滤波器h(n)分解为I个子滤波器h_k(m) = h(k + mI), 
%   其中k=0,1,...,I-1，m为整数索引
% 优势:
%   允许在原始采样率(Fs)处理数据，避免上采样后的高计算量
polyFiltBank = polyDecomp(interpFilt, I);

%% 滤波器组结构示例 (I=2时)
% 假设原型滤波器系数: h = [h0 h1 h2 h3 h4 h5 h6]
% 分解结果:
%   polyFiltBank = [h0 h2 h4 h6;  ← 偶数相位分支(处理原始样点)
%                  h1 h3 h5  0]  ← 奇数相位分支(处理插值样点)
% 注: 末尾补零使各分支长度一致

end
