function y = channelDelayCompensation(x,skew,fADC)

% Delay_COMpensation     Compensate for the skew between channels
% 多通道（Mode） 或者 信号（实部与虚部之间的skew）
% 输入x必须为多通道参量

dt = skew;
fs = fADC;


y = NaN(size(x));  % 初始化输出数组为NaN
for ii = 1:size(x,2)  % 遍历每个通道
    if dt(ii) == 0
        y(:,ii) = x(:,ii);  % 无时滞，直接复制信号
    else
        % 原始时间轴（单位：秒）
        original_time = (0:size(x,1)-1).'/fs;
        % 平移后的时间轴（补偿时滞）
        shifted_time = original_time + dt(ii);
        % 三次插值重采样
        y(:,ii) = interp1(original_time, x(:,ii), shifted_time, 'pchip');
    end
end
