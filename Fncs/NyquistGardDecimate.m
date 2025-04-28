function [out, idx]= NyquistGardDecimate(Ei, param)


if isfield(param, 'SpS_in')
    SpS_in = param.SpS_in;
end
if isfield(param, 'SpS_out')
    SpS_out = param.SpS_out;
end

% 确认输入信号的形式
[M,N1]=size(Ei);
if M < N1
    error('the Ei should be the N × models')
end


r = SpS_in/SpS_out;


if mod(SpS_in, 2), error('Input Nss must be even'); end


N = SpS_in*fix(numel(Ei)/SpS_in);
x = Ei(1:N);
%将数据转化为矩阵，将信号划分为列（每列表示一个符号）
%将输入信号 x 划分成一系列非重叠的窗口，每个窗口的大小是 2*Nss_in + 1，用于后续计算
symbols = buffer(x,2*SpS_in+1, 1, 'nodelay'); % reshape signal into columns (column=symbol)
cstart=SpS_in/2;
%时序误差度量 (err_gard) 和方差误差度量 ， 找最佳采样点
for jj=1:SpS_in
    %表示第一个窗口的信号样本，第二个窗口的信号样本 之间进行 相减 ，并乘上当前采样点的误差， 进行均值误差的计算
    err_gard(jj) = mean(real((symbols(cstart+jj-SpS_in/2, :).*conj(symbols(cstart+jj-SpS_in/2, :))...
        -symbols(cstart+jj+SpS_in/2, :).*conj(symbols(cstart+jj+SpS_in/2, :)))...
        .*conj(symbols(cstart+jj, :)).*symbols(cstart+jj, :)));
    %conj(symbols(cstart+jj, :)).*symbols(cstart+jj, :) 表示当前采样点的幅度平方。
%     var_err_gard(jj) = var(real((symbols(cstart+jj-SpS_in/2, :).*conj(symbols(cstart+jj-SpS_in/2, :))...
%         -symbols(cstart+jj+SpS_in/2, :).*conj(symbols(cstart+jj+SpS_in/2, :)))...
%         .*conj(symbols(cstart+jj, :)).*symbols(cstart+jj, :)));
end

%figure, plot(err_gard, '-o');
%figure, plot(var_err_gard, 'r-o');
%var_gard = sum(var_err_gard)
%             比较时序误差度量的变化，找到时序误差发生零交叉的点（crossPoint）。
% 通过检测误差度量 err_gard 的过零点（即符号切换点）来调整采样点的位置，以找到最佳的采样点位置，并进行下采样
crossPoint = (1:1:N);
%用 符号进行 定义了 ，不再使用误差值进行计算
signum = sign(err_gard);	               % get sign of data
%x 值为0的元素设置为1，这是为了处理 err_gard 中的精确零值（正零）。
signum(x==0) = 1;	                       % set sign of exact data zeros to positive
%用于找到 signum 中出现变化的地方，即 diff(signum)~=0
% 表示 signum 中相邻元素不相等的位置，然后取这些位置中最大的一个。
% 这个值表示了 err_gard 中的零交叉点，即 err_gard 从正数变成负数或从负数变成正数的位置。
ptr = max(crossPoint(diff(signum)~=0));	   % get zero crossings by diff ~= 0

% 在零交叉点附近微调最佳采样点，以获得更精确的采样位置。这部分包括了一些差值和细微调整的操作。
% Fine sample adjust in the crossing point
if abs(err_gard(ptr)) > abs(err_gard(ptr + 1))
    ptr = ptr + 1;
elseif ptr > 1
    if abs(err_gard(ptr)) > abs(err_gard(ptr - 1))
        ptr = ptr - 1;
    end
end
%如果后一个值的绝对值比前一个更大，就将零交叉点位置向后移动一个位置，以找到更接近零的位置。
% 如果前一个值的绝对值更大，就将零交叉点位置向前移动一个位置。
% 最终确定输出信号的采样点索引 idx，以及降采样后的输出信号 out。

% 如果没有找到 零点 的 交叉点 ， 或者 ptr 移位到0，则需要重新判断
if isempty(ptr)
    %找到 err_gard 中的最小值，并将其对应的索引作为采样点位置。
    [~,ptr] = min(abs(err_gard)); 
end
%根据ptr移位
x = circshift(x, -ptr);

idx = 1:r:N;
out = x(idx);
end