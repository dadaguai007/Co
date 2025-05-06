%计算误码率
%bitSeq待统计ber的序列；refSeq参考序列；ber误码率；nErr错误的比特数;g段数
function[ber,nErr,error_location] = CalcBER(bitSeq,refSeq)
bitSeq = bitSeq(:);
refSeq = refSeq(:);
a = min(length(bitSeq),length(refSeq));
bitSeq = bitSeq(end-a+1:end);
refSeq = refSeq(end-a+1:end);
error_location = find(bitSeq ~= refSeq);
nErr = nnz(refSeq-bitSeq);
ber = nErr/length(refSeq);
end