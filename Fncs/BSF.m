function dataout=BSF(datain,samplerate,LowF,HighF,value)
N = numel(datain);
data = reshape(datain, 1, numel(datain));
t = 1:1:N;
freq = samplerate * t / N;
filterResponse = ones(size(freq)); % 初始化滤波器响应为1
% 定义带阻滤波器的频率范围 (LowF 和 HighF 之间的频率将被阻止)
% 将带阻范围内的频率成分置为0
position = find(freq>=LowF&freq<= HighF);
filterResponse(min(position):max(position)) = value;
filterResponse(N-max(position)-1:N-min(position)-1)=value;
% 应用滤波器
dataout = ifft(fft(data) .* filterResponse);
dataout = reshape(dataout, size(datain));
end
