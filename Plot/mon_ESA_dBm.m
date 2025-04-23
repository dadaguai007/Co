function [ IfdBm, Fre ] = mon_ESA_dBm(signal,fs,flag_plot)
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
% If=Efout .* conj(Efout);
If=pnorm(If);
IfdBm = 10*log10(If);
Fre=fftshift(freq)./1e9;
if flag_plot
    figure;
    plot(fftshift(freq)./1e9,IfdBm,'b');
    xlabel('Frequency (GHz)');
    ylabel('Magnitude (dB)');
    box on;
    set(gca, 'FontName', 'Arial', 'FontSize', 14);
    set(gcf,'Position', [0, 0, 480, 400]);
    set(gca, 'LineWidth', 1.25);
    set(gca,'XLim',[-fs/2/1e9 fs/2/1e9],'YLim',[-50 30]);
end

