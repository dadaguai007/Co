function [wout,filtercoeffreq] = efilter(win, type, order, bwl, Fs, sps, outputVoltage)
% 示例，bwl为 Rs*0.44;
% param.efilter = struct('type', 'gauss', 'bw', bwl, 'order', 5);
[L, N] = size(win);
switch lower(type)
    case 'window'
        filtercoeftime = (2*bwl/(Fs))*sinc(2*sps*(bwl/(Fs))*linspace(-floor(L/2-1)/(sps), floor(L/2-1)/(sps), floor(L/2-1)*2+1)).';
        coefshift = -floor(length(filtercoeftime)/2);
        filtercoeftime(L) = 0;
        filtercoeffreq = fftshift(fft(circshift(filtercoeftime(:), [coefshift 0])));
    case 'gauss'
        filtercoeffreq = exp(-log(sqrt(2))*((linspace(-0.5,0.5,L)/(bwl/(Fs))).').^(2*order));
    case 'bessel'
        [b,a] = besself(order, bwl);
        besselcoef = polyval(b, Fs*2j*pi*linspace(-0.5,0.5,L))./polyval(a, Fs*2j*pi*linspace(-0.5,0.5,L));
        filtercoeffreq = exp(1j*angle(besselcoef(:)));
end
% 应用滤波器
for idx = 1:N
    wout(:, idx) = fftshift(fft(win(:, idx)));
    wout(:, idx) = ifft(ifftshift(win(:, idx).*filtercoeffreq(:)));
end
% 输出信号电压缩放
for idx = 1:N
    if outputVoltage
        maxVoltage = max(max(abs([real(wout) imag(wout)])));
        wout(:, idx) = (real(wout(:, idx)) + 1j*imag(wout(:, idx)))*outputVoltage/maxVoltage;
    else
        wout(:, idx) = (real(wout(:, idx)) + 1j*imag(wout(:, idx)));
    end
end

end


%这个MATLAB函数`efilter`的功能是对输入信号`win`应用不同类型的滤波器进行滤波，并可选择对输出信号进行电压缩放。

% 该函数接受7个输入参数：
% - `win`：输入信号矩阵，大小为`L×N`，其中`L`是信号长度，`N`是信号通道数。
% - `type`：滤波器类型，字符串，可选值为`'window'`、`'gauss'`、`'bessel'`。
% - `order`：滤波器阶数，不同类型滤波器的阶数含义不同。
% - `bwl`：滤波器带宽。
% - `Fs`：采样频率。
% - `sps`：每符号采样数。
% - `outputVoltage`：布尔值，用于决定是否对输出信号进行电压缩放。
% 
% ### 计算滤波器系数
% ```
% - 获取输入信号`win`的长度`L`和通道数`N`。
% - 根据`type`选择不同的滤波器类型：
%     - `'window'`：使用`sin(x)/x`（`sinc`）函数生成时域滤波器系数`filtercoeftime`，然后对其进行移位、加零处理，并通过傅里叶变换转换到频域得到`filtercoeffreq`。
%     - `'gauss'`：生成一个基于高斯函数的频域滤波器系数`filtercoeffreq`，通过调整`order`和`bwl`可以控制滤波器的形状和截止频率。
%     - `'bessel'`：使用`besself`函数设计一个贝塞尔滤波器，得到滤波器的分子`b`和分母`a`系数。然后通过`polyval`函数计算在给定频率点上的滤波器响应`besselcoef`，最后取其相位得到`filtercoeffreq`。
% 
% ### 对输入信号进行滤波
% ```
% 通过循环遍历输入信号`win`的每一列：
% - 对每一列信号进行傅里叶变换并移位到零频率在中心的形式。
% - 将频域信号与滤波器系数`filtercoeffreq`相乘。
% - 对乘积结果进行逆傅里叶变换并移回原始频率顺序，得到滤波后的信号。
% 
% ### 输出信号电压缩放
% ```
% 再次循环遍历每一列滤波后的信号：
% - 如果`outputVoltage`为真，计算滤波后信号实部和虚部绝对值的最大值`maxVoltage`，然后将信号乘以`outputVoltage/maxVoltage`进行缩放，使得信号的最大幅值为`outputVoltage`。
% - 如果`outputVoltage`为假，则直接返回滤波后的信号。
