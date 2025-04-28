function n1=EbN0_dB(z1,Eb_N0_dB)

%------------使用EBN0加入噪声----------------------%
% 加噪声
    sig2b=10^(-Eb_N0_dB/10); 
    n1 = sqrt(sig2b/2)*randn(1,length(z1))+1j*sqrt(sig2b/2)*randn(1,length(z1)); % white gaussian noise, 
end