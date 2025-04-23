function [ Efout, freq ] = mon_ESA(signal,fs,flag_plot)
% 调用mon_ESA(signal,fs,1)，即可
if nargin < 3
    flag_plot = 1;
end
%UNTITLED Summary of this function goes hfere
%   Detailed explanation goes here
N =  length(signal);
if mod(N,2)
    N = N-1;
    signal(end) = [];
end
freq = fs/N.*[0:N/2-1,-N/2:-1].';
Efout = fftshift(fft(signal));

If = real(Efout).^2+imag(Efout).^2;
IfdBm = 10*log10(If./(max(abs(If))));
if flag_plot
    plot(fftshift(freq)./1e9,IfdBm);
    xlabel('Frequency (GHz)');
    ylabel('Magnitude (dB)');
end

