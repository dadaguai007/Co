function [H, delay] = H_grpdelay(H, f)
%% Calculate and remove group delay from transfer function H
phiH = unwrap(angle(H));
df = abs(f(1)-f(2));
delay = -diff([phiH 0])/df; % in s
delay = interp1(f, delay, 0);
% Remove the delay
H = H.*exp(1j*2*pi*f*delay);



%
% 这段代码是一个MATLAB函数，名为`H_grpdelay`，用于计算传递函数`H`的群延迟，并将其从传递函数中移除。群延迟是相位响应的导数，它表示信号通过系统时不同频率成分的时间延迟差异。
% 函数`H_grpdelay`接受两个输入参数：
% - `H`：传递函数，通常是一个复数数组，表示系统的频率响应。
% - `f`：频率向量，与传递函数`H`相对应。
% 函数返回两个输出：
% - `H`：移除了群延迟的传递函数。
% - `delay`：群延迟的值，单位是秒。
% 函数的工作流程如下：
% 1. 计算传递函数`H`的相位响应`phiH`，使用`unwrap`函数来去除相位包裹，确保相位变化是连续的。
% 2. 计算 frequency bin 的宽度`df`，这是频率向量`f`中相邻元素之间的差值。
% 3. 计算群延迟`delay`。群延迟是相位响应`phiH`的差分（即变化率），除以频率 bin 的宽度`df`，并乘以-1（因为相位增加表示时间延迟，所以需要取反）。差分操作通过`diff`函数实现，并在相位响应的末尾添加一个0来保持维度一致。
% 4. 使用`interp1`函数在频率为0的位置对群延迟进行插值，因为群延迟通常是在直流（DC）频率处的值。
% 5. 移除传递函数`H`中的群延迟。这通过将`H`乘以一个复指数`exp(1j*2*pi*f*delay)`来实现，其中`1j`是MATLAB中表示虚数单位`i`的符号，`2*pi*f*delay`是相位旋转因子，用于补偿群延迟引起的时间延迟。
% 6. 返回移除了群延迟的传递函数`H`和群延迟的值`delay`。
% 总的来说，这个函数用于分析和校正传递函数的群延迟，以便在信号处理中更准确地模拟系统的行为。
