%%
nn=OFDMQAMN();
nn.DataType='rand';%两个选项：prbs，rand

nn.fft_size = 1024;
nn.nPkts = 100;
nn.nCP = 32;
nn.nModCarriers = 316;
nn.NSym = nn.nModCarriers*nn.nPkts;
nn.nOffsetSub = 4;
nn.order = 2;
nn.M = 16;
nn.prbsOrder = 15;
nn.Rs = 64e9;
nn.Fs = 64e9;
nn.Nsam =nn.Fs/nn.Rs ;
% nn.Nsam=1;
nn.psfRollOff=0.01;
nn.psfLength=256;
nn.psfShape='sqrt';
nn.psfshape='Raised Cosine';
nn.len= (nn.fft_size+nn.nCP)*nn.nPkts; % OFDM 符号长度
nn.dataCarrierIndex=nn.nOffsetSub+(1:nn.nModCarriers);
[y1,y2,upSig,qam_signal,postiveCarrierIndex]=nn.Output();
% 纠正延时
% delay_ps=38.5;
% result = iqdelay(upSig, nn.Fs, delay_ps*1e-12).';
% 创建0,1信号
Signal=zeros(length(upSig),2);
Y=zeros(length(upSig),1);
Y(1:length(Signal)/2)=1;
label = nn.ofdm(qam_signal);