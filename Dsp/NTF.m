function dataout=NTF(datain,samplerate,LowF,HighF,flag_plot)
if nargin<5
    flag_plot=1;
end
% NOTCHFilter 函数用于实现带通滤波器。
% 输入参数：
% datain - 输入的信号数据
% samplerate - 采样率，即每秒采样的次数
% LowF - Notch滤波器的低频截止频率
% HighF - Notch滤波器的高频截止频率
% 输出参数：
% dataout - 滤波后的信号数据
N=numel(datain);
data=reshape(datain,1,numel(datain));
t=1:1:N;
% 计算频率向量
freq=samplerate*t/N;
% 初始化滤波器响应为零
filterResponse=ones(size(freq));
position=find(freq<=HighF&freq>=LowF);
% 设置滤波器响应在LowF和HighF之间的值为0
filterResponse(min(position):max(position))=0;
filterResponse(N-max(position)-1:N-min(position)-1)=0;
% 应用滤波器
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