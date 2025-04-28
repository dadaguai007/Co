function [Eo, timing_values] = Test_gardnerClockRecovery(Ei, param)
% Perform clock recovery using Gardner's algorithm with a loop PI filter.

% Check and set default values for input parameters
% kp = getfield(param, 'kp', 1e-3);
% ki = getfield(param, 'ki', 1e-6);
% isNyquist = getfield(param, 'isNyquist', true);
% returnTiming = getfield(param, 'returnTiming', false);

% 初始化参数
kp = 1e-3;
ki = 1e-6;
isNyquist = 'true';
returnTiming = 'false';
lpad = 1;
maxPPM = 500;


if isfield(param, 'kp')
    kp = param.kp;
end

if isfield(param, 'ki')
    ki = param.ki;
end

if isfield(param, 'isNyquist')
    isNyquist = param.isNyquist;
end

if isfield(param, 'returnTiming')
    returnTiming = param.returnTiming;
end

if isfield(param, 'lpad')
    lpad = param.lpad;
end
if isfield(param, 'maxPPM')
    %     预期的最大时钟频率偏差，以PPM表示，默认值为500。
    maxPPM = param.maxPPM;
end


if size(Ei,1)<size(Ei,2)
    % 保证每一列为一个模式的信号
    Ei=Ei.';
end
% Initializing variables:
%    padarray
%    在信号的末尾填零
Ei = padarray(Ei, [lpad,0], 0, 'post');
[nSamples, nModes] = size(Ei);


% Gardner时钟恢复函数

% 初始化输出向量
% 使用1e6来将ppm转换为比例值（因为 1 ppm 等于 1/1,000,000 或 1e-6）。
% 1 - maxPPM / 1e6 计算的是最大时钟偏差对应的样本数占总样本数的比例。
% 然后，这个比例乘以 nSamples 得到可能的最大样本偏差数。
% 最后，使用 round() 函数将结果转换为整数，确保输出数组的大小是一个整数值。
Ln = round((1 - maxPPM / 1e6) * nSamples);
Eo = zeros(Ln, nModes);

% 初始化t_nco值
t_nco_values = zeros(size(Eo));
last_n = 0;

% 进行时钟恢复
for indMode = 1:nModes
    intPart = 0;%积分部分的初始值
    t_nco = 0;% NCO 时钟偏差的初始值
    n = 3;%样本索引
    m = 3;%插值样本索引
    while n < Ln && m < nSamples - 1
        %插值
        Eo(n, indMode) = interpolator(Ei(m-2:m+1, indMode), t_nco);
        if mod(n, 2) == 1 % 为奇数
            if strcmp(isNyquist,'True')
                ted = gardnerTEDnyquist(Ei(n-2:n, indMode));
            else
                ted = gardnerTED(Ei(n-2:n, indMode));
            end
            % 循环PI滤波器
            intPart = ki * ted + intPart;
            propPart = kp * ted;
            loopFilterOut = propPart + intPart;
            t_nco = t_nco - loopFilterOut;
        end
        % NCO时钟间隔
        if t_nco > 1
            %NCO的时钟比接收信号的时钟快，需要回退一个样本以对齐信号
            t_nco = t_nco - 1;
            n = n - 1;
        elseif t_nco < -1
            %时钟慢，前进一个样本
            t_nco = t_nco + 1;
            n = n + 2;
            m = m + 1;
        else
            % 时钟可接受
            n = n + 1;
            m = m + 1;
        end
        % 记录
        t_nco_values(n, indMode) = t_nco;
    end
    if n > last_n
        last_n = n;
    end
end
Eo = Eo(1:last_n, :);

% 计算算法估计的采样频谱偏移
ppm = calcClockDrift(t_nco_values);
for i=1:nModes
    fprintf('Estimated clock drift: %.2f ppm\n',ppm);
end

if strcmp(returnTiming,'True')
    timing_values=t_nco_values;
end

% Gardner TED函数
    function ted = gardnerTED(x)
        ted = real(conj(x(2)) .* (x(3) - x(1)));
    end

% Gardner TED Nyquist函数
    function ted = gardnerTEDnyquist(x)
        ted = abs(x(2)).^2 .* (abs(x(1)).^2 - abs(x(3)).^2);
    end

% Farrow结构的三次插值器
    function y = interpolator(x, t)
        y = x(1) * (-1/6 * t.^3 + 1/6 * t) ...
            + x(2) * (1/2 * t.^3 + 1/2 * t.^2 - t) ...
            + x(3) * (-1/2 * t.^3 - t.^2 + 1/2 * t + 1) ...
            + x(4) * (1/6 * t.^3 + 1/2 * t.^2 + 1/3 * t);
    end


    function ppm = calcClockDrift(t_nco_values)

        % 计算时钟漂移
        timingError = t_nco_values - mean(t_nco_values, 1);% 对列进行求平均
        t = (0:length(t_nco_values)-1)';

        ppm = zeros(nModes, 1);
        for indMode = 1:nModes
            %计算时间误差序列的连续差分，得到时间误差的变化率。峰值大于0.5
            [~,locs] = findpeaks(abs(diff(timingError(:,indMode))), 'MinPeakHeight', 0.2);
            mean_period = mean(diff(t(locs))); % 平均周期
            %频率
            fo = 1 / mean_period;
            % 时间偏差平均趋势
            ppm(indMode) = fo * 1e6 * sign(mean(t_nco_values(:, indMode)));
        end
    end




end