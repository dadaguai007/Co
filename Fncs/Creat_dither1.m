function out = Creat_dither1(SampleRate,SignalFrequency,N)

% 设置采样率和信号频率
% SampleRate = 60e9; % 采样率为60 Hz
% SignalFrequency = 10e6; % 信号频率为10 MHz

% 信号时长
% N=0.88*25;%Samples
Duration = N/SignalFrequency; 
% 主要计算Duration的长度
% 创建时间网格
dt = 1 / SampleRate;

% 创建频率网格
df = 1 / Duration;

% 设置时间起点
t0 = 0;

% 设置信号时长（网格点单位）
% T = 1 / (dt *df);
T = SampleRate/SignalFrequency*N;
% 设置结构的带宽（网格点单位）
fs = round(SampleRate / df);
% 生成时间序列
t = (t0:T-1) * dt;
% 生成正弦信号
amplitude = 1; % 振幅为1
out = amplitude * sin(2 * pi * SignalFrequency * t);


end