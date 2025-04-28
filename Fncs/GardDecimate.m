function [out, idx]= GardDecimate(Ei,param)

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

%输入采样率和输出采样率的比值，不是整数进行报错
r = SpS_in/SpS_out;

if mod(SpS_in, 2), error('Input Nss must be even'); end


%对输入信号进行截取，将输入信号长度约束为 SpS_in 的整数倍
N = SpS_in*fix(numel(Ei)/SpS_in);
x = Ei(1:N);
symbols = buffer(x,2*SpS_in+1, 1, 'nodelay'); % reshape signal into columns (column=symbol)

%符号长度的一半，起始位置
cstart=SpS_in/2;
%通过比较两个相邻符号样本的特定差值来计算的。这个值表示符号之间的时序误差
%计算了两个符号样本之间的差值的实部的平均值，这是 Gardner 标准的一部分
for jj=1:SpS_in
    err_gard(jj) = mean(real((symbols(cstart+jj-SpS_in/2, :)-symbols(cstart+jj+SpS_in/2, :)).*conj(symbols(cstart+jj, :))));
    % var_err_gard(jj) = var(real((symbols(cstart+jj-Nss_in/2, :)-symbols(cstart+jj+Nss_in/2, :)).*conj(symbols(cstart+jj, :))));
end
%找到最佳采样点
[~,ptr] = min(abs(err_gard)); % find maximum variance point

ptr = mod(ptr-1,r)+1;
%根据找到的最佳采样点和采样率比例，进行降采样
idx = ptr:r:N;
out = x(idx);
symbols = symbols.';
end