%信号同步
function[bit_sync,ref_sync] = sync(bitSeq,refSeq)
bitSeq = bitSeq(:);
refSeq = refSeq(:);
%使用互相关函数 xcorr 计算输入信号序列与参考信号序列的互相关。
% 通过减去均值，可以消除直流分量的影响。
[d,lags] = xcorr(bitSeq,refSeq);  %减去均值
% [psk,locs] = findpeaks(abs(d),'MinPeakHeight',max(abs(d))*0.90);  %寻找峰值
%设置阈值为 0.9，用于确定互相关结果中的峰值。
th = 0.9;
%找出互相关结果中绝对值大于阈值乘以最大绝对值的部分，即大于阈值的峰值位置
locs = find(abs(d)>th*max(abs(d)));
% figure;plot(lags,abs(d))
%计算峰值位置与输入信号序列长度的差，并加1，以获得信号同步的偏移量。
f = locs-max(length(bitSeq),length(refSeq))+1;   %峰值index
%找出偏移量中大于零的部分，即有效的同步偏移量。


if f>0
f_idx = find(f>0);
ff = f(f_idx);
%从输入的信号序列中，根据第一个有效同步偏移量开始截取同步后的信号序列。
bit_sync = bitSeq(ff(1):end);
ref_sync = refSeq;
else
bit_sync = bitSeq;
d =length(refSeq);
ref_sync = refSeq(d-locs+1:end);

end

%xcorr [C, lags] = xcorr(x, y)
% C为互相关结果，是一个向量，长度为 (2*maxlag+1)，maxlag 是互相关时使用的最大滞后值。
% 互相关结果中的每个元素表示 x 与 y 之间在不同滞后值下的相关性。
%lags 是对应的滞后值向量，它表示每个互相关结果的滞后值。
