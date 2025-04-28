function dataout=LPF(datain,samplerate,bandwidth,flag_plot)
if nargin<4
    flag_plot=1;
end
% LPF 函数用于实现低通滤波器。
% 输入参数：
% datain - 输入的信号数据
% samplerate - 采样率，即每秒采样的次数
% bandwidth - 低通滤波器的截止频率
% 输出参数：
% dataout - 滤波后的信号数据
N=numel(datain);
data=reshape(datain,1,numel(datain));
t=1:1:N;
% 计算频率向量
freq=samplerate*t/N;
% 初始化滤波器响应为零
filterResponse=zeros(size(freq));
position=find(freq<=bandwidth);
% 设置滤波器响应在截止频率以下的值为1
filterResponse(min(position):max(position))=1;
filterResponse(N-max(position)-1:N)=1;
dataout=ifft(fft(data).*filterResponse);
dataout=reshape(dataout,size(datain));

if flag_plot
    figure;
    plot(freq./1e9,filterResponse,'b');
    xlabel('Frequency (GHz)');
    ylabel('Magnitude (dB)');
    box on;
    set(gca, 'FontName', 'Arial', 'FontSize', 14);
    set(gca, 'LineWidth', 1.25);
end
end