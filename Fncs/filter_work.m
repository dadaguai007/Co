% 具体使用，根据之前设计的频域滤波器
function     out_signal=filter_work(in,HFilter)
%filter 工作函数：用于从频域设计的滤波器对信号产生作用
in = fftshift(fft(in(:)));
out_signal = ifft(ifftshift(in.*HFilter(:)));
end