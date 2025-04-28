function [out, idx]= SLNDecimate(Ei, param)
%基于SLN（Squared-Linear Norm）标准的信号降采样功能

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


% 降采样
r = SpS_in/SpS_out;

if mod(SpS_in, 4), error('Input Nss must be a multiple of 4'); end

%对输入信号进行截取，将输入信号长度约束为 SpS_in 的整数倍
N = SpS_in*fix(numel(Ei)/SpS_in);
x= Ei(1:N);
symbols = buffer(x,2*SpS_in+1, 1, 'nodelay').'; % reshape signal into columns (column=symbol)
%SLN 标准中的误差度量
%权重向量，包含了四个复数权重。这些权重是 SLN 标准中用于非线性处理的部分
%angle非线性平方处理后的信号的相位角度
%所选样本的相位角度的平均值
%计算了 SLN 标准的误差度量 err_SLN。它使用权重向量和非线性处理来计算每个符号的误差。
% 具体来说，它计算每个符号中不同位置的加权相位角度的平均值。
for jj=1:SpS_in
    %symbols(:, jj:Nss_in/4:jj+Nss_in-1)：这一部分从 symbols 中选择一列符号，并从中提取一定间隔的子序列。这里使用 Nss_in/4 作为间隔，从 jj 开始，以 jj+Nss_in-1 结束。
    % 这样做是为了选择每个符号中的不同位置，以进行非线性处理。
    %                 abs(...).^2：对上述选择的子序列执行绝对值平方操作，以获得每个位置处的信号强度的平方
    %exp 是 计算一个复数权重向量。它创建了一个包含四个复数权重的向量
    %在 SLN 标准的误差度量中，相位角度是重要的，因为反映了信号在不同位置上的相位信息
    err_SLN(jj) = mean(angle(abs(symbols(:, jj:SpS_in/4:jj+SpS_in-1)).^2*exp(-1i*.5*pi*(0:3)).'));
end
[~,ptr] = min(abs(err_SLN)); % find maximum variance point


ptr = mod(ptr-1,r)+1;
idx = ptr:r:N;
out = x(idx);
end