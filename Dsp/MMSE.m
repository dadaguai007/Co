% MMSE 均衡
function S1_mmse=MMSE(H1,y_fft,Eb_N0_dB,s1)
sigs12=var(s1);% 调制信号的方差
sig2b=10^(-Eb_N0_dB/10);
% Egaliseur MMSE
W1_mmse = conj(H1)./(abs(H1).^2 + (sig2b/sigs12)); % Egaliseur signal 1
% 信道系数的共轭， 噪声功率，信号方差
% MMSE 均衡器的目标是最小化均衡器输出信号与原始信号之间的均方误差
% Egalisation
S1_mmse = diag(W1_mmse)*y_fft;

end