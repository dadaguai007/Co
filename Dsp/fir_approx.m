function [hfir, Ntaps] = fir_approx(Hch, fs, Efraction, N, Mct)
%% Calculate FIR approximation of channel given by anonymous function Hch
% - Inputs:
% - Hch: anonymous function of the channel frequency response and the Hch
% can be model as Discrete data
% - fs: sampling frequency of FIR filter
% - Efraction: fraction of energy that must be contained in samples of the
% channel impulse response
% - N (optional, default=1024): number of samples to use in calculations
% - Mct (optional, default=5): oversampling ratio to emulate continuous
% time

if nargin < 4
    N = 1024;
    Mct = 5;
end

% 频率向量和逆FFT变化，到时域响应
f = freq_time_set(N, Mct*fs);
ht = fftshift(ifft(ifftshift(Hch(f))));
t0 = N/2+1;
% 正半样本和负半样本
hsamp_pos = ht(t0:Mct:end);
hsamp_neg = ht(t0:-Mct:1);
% 累积能量计算
Sp = cumsum(abs(hsamp_pos).^2)/sum(abs(hsamp_pos).^2);
Sn = cumsum(abs(hsamp_neg).^2)/sum(abs(hsamp_neg).^2);
% 找到能量分数的索引
idxp = find(Sp >= Efraction, 1);
idxn = find(Sn >= Efraction, 1);
%确定Tap数和FIR滤波器系数
Ntaps = idxp + idxn + 1;
hfir = ht(t0 + Mct*(-idxn:idxp));


%
% 这段代码是一个MATLAB函数，名为`fir_approx`，用于计算给定匿名函数`Hch`表示的通道频率响应的FIR（有限冲激响应）滤波器近似。该函数的目的是设计一个FIR滤波器，其冲激响应逼近于给定的通道响应。
% 函数`fir_approx`接受五个输入参数：
% - `Hch`：一个匿名函数，表示通道的频率响应。这个函数接受频率向量作为输入，并返回对应的通道频率响应值。
% - `fs`：FIR滤波器的采样频率。
% - `Efraction`：冲激响应样本中必须包含的能量分数。这是一个介于0和1之间的数值，表示在设计FIR滤波器时，希望保留的通道能量的比例。
% - `N`（可选，默认值为1024）：用于计算的样本数。
% - `Mct`（可选，默认值为5）：模拟连续时间的过采样比率。
% 函数返回两个输出：
% - `hfir`：设计的FIR滤波器的系数。
% - `Ntaps`：FIR滤波器的抽头数。
% 函数的工作流程如下：
% 1. 检查输入参数的数量，如果少于4个，则将`N`和`Mct`设置为默认值。
% 2. 使用`freq_time`函数（未在代码中定义，可能是一个辅助函数）生成频率向量`f`，并根据`Hch`计算通道的时域冲激响应`ht`。
% 3. 计算`ht`的正半样本`hsamp_pos`和负半样本`hsamp_neg`，这通过对冲激响应进行采样来实现，采样间隔由`Mct`决定。
% 4. 计算正半样本和负半样本的累积能量`Sp`和`Sn`，这是通过计算样本的平方的累加和，然后除以总能量来实现的。
% 5. 找到能量分数`Efraction`的索引`idxp`和`idxn`，这些索引表示在累积能量超过`Efraction`时的样本位置。
% 6. 确定FIR滤波器的抽头数`Ntaps`，这是通过将正半样本和负半样本的索引相加再加1来实现的。
% 7. 根据确定的抽头数`Ntaps`提取FIR滤波器的系数`hfir`，这些系数是从冲激响应`ht`中提取的，从`t0 + Mct*(-idxn)`开始，到`t0 + Mct*idxp`结束。
% 总的来说，这个函数通过计算通道的时域冲激响应，并根据给定的能量分数确定FIR滤波器的抽头数和系数，从而设计出一个逼近通道响应的FIR滤波器。
