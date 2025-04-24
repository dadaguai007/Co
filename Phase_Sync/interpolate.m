function [xI] = interpolate(method, x, m_k, mu, b_mtx, poly_f)
% 多方法插值计算
% 输入:
%   method - 插值方法: 0)多相 1)线性 2)二次 3)三次
%   x      - 输入信号向量
%   m_k    - 基带点索引（插值区间左端）
%   mu     - 分数间隔（0≤mu<1）
%   b_mtx  - 多项式插值系数矩阵（二次/三次时需提供）
%   poly_f - 多相滤波器组（多相插值时需提供）
% 输出:
%   xI     - 插值结果

    % μ越界处理（支持奇数L的偏移补偿）
    if mu < 0
        m_k = m_k - 1; 
        mu = mu + 1;
    elseif mu >= 1
        m_k = m_k + 1; 
        mu = mu - 1;
    end
    assert(mu >= 0 && mu < 1); % 确保μ在合法范围

    switch method
        case 0 % 多相插值
            polyInterpFactor = size(poly_f, 1);
            phase = floor(mu * polyInterpFactor) + 1; % 选择相位
            taps = poly_f(phase, :); % 获取对应子滤波器
            xI = taps * x(m_k - length(taps) + 1 : m_k); % 卷积计算
            
        case 1 % 线性插值
            xI = mu * x(m_k + 1) + (1 - mu) * x(m_k);
            
        case 2 % 二次插值（Farrow结构）
            v_l = x(m_k - 1 : m_k + 2) * b_mtx; 
            xI = (v_l(3)*mu + v_l(2)) * mu + v_l(1);
            
        case 3 % 三次插值（Farrow结构）
            v_l = x(m_k - 1 : m_k + 2) * b_mtx;
            xI = ((v_l(4)*mu + v_l(3)) * mu + v_l(2)) * mu + v_l(1);
    end
end
