function sigRx = gsop(rLn)
% Gram-Schmidt正交化处理 (用于补偿I/Q不平衡)
%
% 输入参数:
%   rLn : 复数数组，输入信号（实部为I分量，虚部为Q分量）
%
% 输出参数:
%   sigRx : 正交化后的复数信号
%
% 参考文献:
%   [1] Digital Coherent Optical Systems, Architecture and Algorithms. 
%   [2] I. Fatadin, S.J. Savory, D. Ives, Compensation of quadrature imbalance 
%       in an optical QPSK coherent receiver. IEEE Photon. Technol. Lett. 20(20), 1733–1735 (2008)

% 输入信号为列向量
% 将复数信号拆分为实部(I)和虚部(Q)，组成两列矩阵
% Rin的第一列是I分量，第二列是Q分量
Rin = [real(rLn), imag(rLn)];

% 步骤1: 对Q分量进行归一化（除以Q分量的均方根值）
rQOrt = Rin(:,2) / sqrt(mean(Rin(:,2).^2));

% 步骤2: 消除I分量中的Q相关性（正交化核心操作）
% 计算投影系数：E(I·Q)/E(Q²)
coeff = mean(Rin(:,1) .* Rin(:,2)) / mean(Rin(:,2).^2);
% 从I分量中减去其在Q分量上的投影
rIInt = Rin(:,1) - coeff * Rin(:,2);

% 步骤3: 对正交化后的I分量进行归一化
rIOrt = rIInt / sqrt(mean(rIInt.^2));

% 组合为复数信号：正交化I分量 + j*归一化Q分量
sigRx = rIOrt + 1j * rQOrt;
end
