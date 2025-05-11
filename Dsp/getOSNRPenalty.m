function reqOSNR = getOSNRPenalty(OSNRs, BERs, tBER)
%OSNR_PENALTY 光信噪比代价计算函数
% 本函数通过线性插值法在BER-OSNR曲线上计算达到目标误码率所需的OSNR值
% 适用场景：光通信系统性能评估、接收机灵敏度测试等
% 参考理论：基于误码率与光信噪比的单调关系进行线性插值
% 可以处理多条BER曲线
% tBER为目标 误码率

%% 初始化输出矩阵为NaN（处理无效曲线情况）
reqOSNR = NaN(size(BERs,2), 1); % M×1列向量，M为BER曲线数量

%% 数据预处理（提高计算鲁棒性）
if tBER > 0
    % 过滤异常高BER值：当目标BER>0时，认为BER>0.3的数据不可靠
    BERs(BERs > 0.3) = NaN; 
end
% 处理无穷大值（可能由误码率计算错误导致）
BERs(abs(BERs) == Inf) = NaN; 

%% 遍历每条BER曲线进行处理
for n = 1:size(BERs,2)
    % 提取当前处理曲线的BER数据
    test_ber = BERs(:,n); 
    
    % ========== 确定曲线单调性方向 ==========
    % 计算BER差分均值判断整体趋势：
    % 正常系统BER应随OSNR增加而降低（负差分）
    if mean(diff(test_ber)) <= 0 
        % 情况1：单调递减（正常情况）
        % 寻找最后一个高于目标BER的点（左侧插值点）
        i_before = find(test_ber > tBER, 1, 'last'); 
        % 寻找第一个低于目标BER的点（右侧插值点）
        i_after  = find(test_ber < tBER, 1);         
    else
        % 情况2：单调递增（异常情况，需检查测试数据方向）
        % 寻找最后一个低于目标BER的点
        i_before = find(test_ber < tBER, 1, 'last'); 
        % 寻找第一个高于目标BER的点
        i_after  = find(test_ber > tBER, 1);         
    end
    
    % ========== 边界条件检查 ==========
    % 检查插值点是否存在且有效
    if isempty(i_after) || isempty(i_before) ||...       % 未找到插值点
       i_before >= size(test_ber,1) ||...                % 左侧点超出索引
       i_after > size(test_ber,1)                        % 右侧点超出索引
        continue; % 跳过当前曲线（保持NaN）
    end
    
    % ========== 线性插值计算 ==========
    % 公式：OSNR_req = OSNR_before + (OSNR_after - OSNR_before) * 
    %                (tBER - BER_before)/(BER_after - BER_before)
    reqOSNR(n) = OSNRs(i_before) + ...                   % 基准OSNR值
        (OSNRs(i_after) - OSNRs(i_before)) * ...         % OSNR差值
        (tBER - test_ber(i_before)) / ...                % BER差值分子
        (test_ber(i_after) - test_ber(i_before));        % BER差值分母
end
