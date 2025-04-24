%信号同步
function[bit_sync,ref_sync,ff] = sync(bitSeq,refSeq)
bitSeq = bitSeq(:);
refSeq = refSeq(:);
[d,lags] = xcorr(bitSeq-mean(bitSeq),refSeq-mean(refSeq));  %减去均值
% [psk,locs] = findpeaks(abs(d),'MinPeakHeight',max(abs(d))*0.90);  %寻找峰值
th = 0.9;
locs = find(abs(d)>th*max(abs(d)));
% figure;plot(lags,abs(d))
f = locs-max(length(bitSeq),length(refSeq))+1;   %峰值index
f_idx = find(f>0);
ff = f(f_idx);
bit_sync = bitSeq(ff(1):end);
ref_sync = refSeq;
end