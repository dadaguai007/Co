  function[data] = clk_recovery(sig,sps,osr,fb,fs)
%数字时钟恢复函数
%sig 输入信号
%sps sample per symbol in dsp
%osr oversampling rate (compared with baud rate) 过采样率
%fb baud rate
%fs sampling rate

data = resample(sig,osr*fb,fs);  %四倍上采样，为了满足后续时钟恢复频率
f_clk_rec = cr(data,fb,osr,0); % 注意：f_clk_rec是对fb的估计，不是对fb*osr的估计！
t_cur = 1/osr/fb*(0:length(data)-1); % 注意：插值前data1的采样率是osr*fb
t_new = 0:1/f_clk_rec/osr:t_cur(end); % 注意：插值后data2的采样率是osr*f_clk_rec
data = interp1(t_cur,data,t_new);
data = resample(data,sps,osr);
end