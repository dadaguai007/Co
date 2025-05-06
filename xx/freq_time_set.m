function [f, t] = freq_time_set(N, fs)
% Generate frequency and time measures
% N: data length
dt = 1/fs;
t = 0:dt:(N-1)*dt;
df = 1/(dt*N);
f = -fs/2:df:fs/2-df;

end