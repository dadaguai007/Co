function [compSig] = KK_Bo(rxSig)
%bo kk
% c=1;
% rxSig=DeWaveform+c;

CC       = mean(rxSig);
data1    = sqrt(rxSig)/CC;
data2    = rxSig/(2*(CC^2));
data3    = 2*data1-data2;
data5    = (abs(data3))./2;
SSB_real = (data1-data5)*CC;
compSig = hilbert(SSB_real);

end