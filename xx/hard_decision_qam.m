function rxSignal=hard_decision_qam(M,R)
% 硬判决
const = qammod(0:M-1,M);
% 归一化
const=pnorm(const);
for i = 1:length(R)
    distances = abs(R(i) - const).^2; % 计算接收信号点到所有星座点的距离
    [~, index] = min(distances); % 找到距离最近的星座点的索引
    rxSymbols(i) = index; % 记录判决结果
end
% 恢复出的信号
rxSignal= const(rxSymbols);
end