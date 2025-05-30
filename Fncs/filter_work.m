% 具体使用，根据之前设计的频域滤波器
% 逻辑没有错误，需要对应使用相应的，从频域设计的滤波器
function     out_signal=filter_work(in,HFilter)
%filter 工作函数：用于从频域设计的滤波器对信号产生作用(不需要对滤波器进行相应的变换)
in = fftshift(fft(in(:)));
out_signal = ifft(ifftshift(in.*HFilter(:)));
end