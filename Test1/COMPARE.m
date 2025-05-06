clear
clc
txBits = randi([0,1],1,100000);
M = 16;
C = qammod([0:M-1],M);
SNR = 7:0.5:8;
SNR=[SNR,8.5:0.1:8.6];
SNR=[SNR,8.6:0.01:8.65];
SNR = [SNR,9:15];
FEC.rate = 2/3;
error = ones(1,numel(SNR));
for nEN = 1 : numel(SNR)
    %% LDPC编码后的误码率
   [encBits,Stx,txSyms,txBits,FEC] = LDPC_encoder_QAM(txBits,FEC,C);
   Srx = awgn(Stx,SNR(nEN),'measured');
   N0 = getN0_MMSE(Stx,Srx);
   Prx = var(Srx.');
   N0 = N0 .* Prx;
   [LLRs] = LLR_eval(Srx,N0,C);
   [rxBits,nBlocks] = LDPC_decoder(LLRs,FEC.LDPC_enc.ParityCheckMatrix,FEC.idx,50);
   error(nEN) = sum(rxBits~=txBits)/numel(LLRs);
   %% 未变码误码率计算
   Stx_uncode = qammod(txBits',M,'InputType','bit');
   Srx_uncode = awgn(Stx_uncode,SNR(nEN),'measured');
   rxBits = qamdemod(Srx_uncode,M,'OutputType','bit');
   error_uncode(nEN) = sum(rxBits'~=txBits)/numel(txBits);
   
end
semilogy(SNR,error,'+-')  
hold on
semilogy(SNR,error_uncode,'*-')
