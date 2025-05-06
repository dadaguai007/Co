%%
%Initialization ��ʼ��
clear;
clear global;
close all;
clc;
%%
% Load Libraries:    ���ؿ⣺
%addpath(genpath('..\..\..\'));
addpath(genpath('C:\Users\Administrator\Desktop\work\emulation\3_OptDSP_lite-master - 2'));
%%
global PROG;               %������PROG����������Ϊȫ��
PROG.showMessagesLevel = 2;
initProg();
RGB = fancyColors();
%%
%Input Parameters  �������
SIG.M = 64;                 % QAM constellation size             QAM ������С
% SIG.M = 16; %�Լ��Ķ�1
SIG.symRate = 60e9;         % total symbol-rate of the signal    �źŵ��ܷ�����
SIG.bitRate_net = 200e9;    % net bit-rate [Gbps]                �������� [Gbps]
SIG.modulation = 'QAM';     % modulation type [QAM/PAM]          �������� [QAM / PAM]
SIG.rollOff = 0.1;          % roll-off factor                    ��������
SIG.nPol = 1;               % number of polarizations            ������
SIG.nSyms = 2^17;           % total number of simulated symbols  ģ���������
nSpS = 2;                   % number of samples per symbol       ÿ�����ŵ�������
laserLW = 0e6;              % laser linewidth [Hz]               �����߿� [Hz]

FEC_rate = 5/6;             % FEC rate                           FEC��

pilotRate = 1;              % Pilot rate                         �Ե���
useCPE2 = false;            % flag to decide whether to use 2nd CPE or not  �����Ƿ�ʹ�õڶ��� CPE �ı�־

SNR_dB = 12;                % SNR in dB                          ����� (dB)
%%
%Set Transmitter Parameters ���÷��Ͷ˲���
%Signal Parameters:         �źŲ�����
nBpS_net = SIG.bitRate_net / (SIG.nPol*SIG.symRate*FEC_rate*pilotRate);

TX.SIG = setSignalParams('symRate',SIG.symRate,'M',SIG.M,...
    'nPol',SIG.nPol,'nBpS',nBpS_net,'nSyms',SIG.nSyms,...
    'roll-off',SIG.rollOff,'modulation',SIG.modulation);

%Modulation Parameters:    ���Ʋ�����
TX.QAM = QAM_config(TX.SIG);

%Bit Parameters:  λ������
TX.BIT.source = 'randi';
TX.BIT.seed = 1;

%Pulse Shaping Filter Parameters:
TX.PS.type = 'RRC';
TX.PS.rollOff = TX.SIG.rollOff;
TX.PS.nTaps = 128;

%DAC Parameters:
TX.DAC.RESAMP.sampRate = nSpS*TX.SIG.symRate;

%LASER Parameters:
TX.LASER.linewidth = laserLW;

%DSP Pilots Parameters:8  DSP ��Ƶ������8
TX.PILOTS.active = true;
TX.PILOTS.rate = pilotRate;
TX.PILOTS.option = 'outerQPSK';
% TX.PILOTS.scaleFactor = 1;

%FEC Parameters:  ǰ����������
TX.FEC.active = false;
TX.FEC.rate = FEC_rate;
TX.FEC.nIter = 50;

%PCS Parameters:
TX.PCS.method = 'CCDM';

%Generate Tx Bits  ���� Tx ����s
TX.BIT.txBits = Tx_generateBits(SIG.nSyms,TX.QAM.M,TX.QAM.nPol,TX.BIT);%����������TX.BIT.txBits

%Generate Transmitted Symbols               
[S.tx,txSyms,TX.QAM] = Tx_ProbShaping(TX.BIT.txBits,TX.QAM,TX.SIG,TX.FEC.rate);%����������S.tx��txSyms��Tx.QAM.symProb��TX.QAM�ɲ鿴���

%Plot Transmitted Constellation Symbols:
plotProbShaping_PDF_const('const',TX.QAM.IQmap,'symbols',txSyms(:),'color',RGB.itred);  

%Add DSP Pilots  ��� DSP �Ե�
if isfield(TX,'PILOTS') && TX.PILOTS.active && TX.PILOTS.rate < 1
    [S.tx,TX.PILOTS] = Tx_addPilots(S.tx,TX.PILOTS,TX.QAM.IQmap);
    TX.QAM.meanConstPower = mean(abs(S.tx).^2,2);
    TX.QAM.maxConstPower = max(abs(S.tx).^2,[],2);
end

%Pulse Shaping
[S.txSC,TX.PS] = pulseShaper(S.tx,nSpS,TX.PS);

%Plot Pulse Shaping Taps �����������γ�ͷ
figure();
plot(TX.PS.W);
title('Pulse Shaping Taps')

%Plot Signal Spectrum
figure();
pwelch(S.txSC(1,:),1e4,[],[],TX.DAC.RESAMP.sampRate,'centered')

%Digital-to-Analog Converter
Fs_SC = nSpS * SIG.symRate;
[S.txSC,TX.DAC] = Tx_DAC(S.txSC,TX.DAC,Fs_SC);

%Apply Laser
nSamples = size(S.tx,2) * nSpS;
[Scw,TX.LASER] = laserCW(TX.LASER,Fs_SC,nSamples);
nPol = size(S.txSC,1);
if nPol == 2
    Scw = Scw/sqrt(2);      % divide CW power by 2 if dual-pol
end
for n = 1:nPol
    S.txSC(n,:) = Scw.*S.txSC(n,:);
end

%Plot Laser Phase Noise
figure();
plot(TX.LASER.phaseNoise)
title('Transmitted Laser Phase Noise','interp','latex')

%Set SNR
S.rx = setSNR(S.txSC,SNR_dB,TX.DAC.RESAMP.sampRate,SIG.symRate);

%Plot Received Signal Spectrum:
Sf = pwelch(S.rx,1e3,[],[],Fs_SC,'centered','psd').';
[SNR_rx_dB,Ps,Pn] = SNR_estimateFromSpectrum(Sf,Fs_SC,TX.SIG.symRate,...
    TX.SIG.rollOff,[-0.9 -0.7; 0.7 0.9],0,true);

%Set DSP Parameters
%% Matched-Filter:
DSP.MF.type = 'RRC';
DSP.MF.rollOff = TX.SIG.rollOff;

%% Carrier-Phase Estimation - 1st Stage:
DSP.CPE1.method = 'pilot-based:optimized';
DSP.CPE1.decision = 'data-aided';
% DSP.CPE1.nTaps = 15;
DSP.CPE1.nTaps_min = 1;
DSP.CPE1.nTaps_max = 201;
DSP.CPE1.PILOTS = TX.PILOTS;

%% Carrier-Phase Estimation - 1st Stage:
% DSP.CPE2.method = 'BPS:optimized';
DSP.CPE2.method = 'BPS';
DSP.CPE2.nTaps = 22;
DSP.CPE2.nTaps_min = 1;
DSP.CPE2.nTaps_max = 501;
DSP.CPE2.nTestPhases = 10;
DSP.CPE2.angleInterval = pi/8;

%% Demapper:
DSP.DEMAPPER.normMethod = 'MMSE';

%Apply Matched Filter
S.rx = LPF_apply(S.rx,DSP.MF,TX.DAC.RESAMP.sampRate,SIG.symRate);

%Downsample to 1 Sample/Symbol
S.rx_1sps = S.rx(1:nSpS:end);

%Apply 1st-Stage CPE: Pilot-Based
if pilotRate < 1
    [S.rx_1sps,DSP.CPE1] = carrierPhaseEstimation(S.rx_1sps,S.tx,DSP.CPE1);
end

%Apply 2nd-Stage CPE: Blind-Phase Search
C = TX.QAM.IQmap;
if useCPE2
    [S.rx_1sps,DSP.CPE2] = carrierPhaseEstimation(S.rx_1sps,S.tx,DSP.CPE2,C);
end

%Remove DSP Pilots
if pilotRate < 1
    [S.rx_1sps,S.tx] = pilotSymbols_rmv(S.rx_1sps,S.tx,DSP.CPE1.PILOTS);
end

%Apply Symbol Demapper
[DSP.DEMAPPER,S.tx] = symDemapper(S.rx_1sps,S.tx,C,DSP.DEMAPPER);

%Calculate Error Vector Magnitude (EVM):
EVM = EVM_eval(S.rx_1sps,S.tx);

% ����5��ԭ����ע�� ����
% %Calculate Mutual Information (MI) and Generalized Mutual Information (GMI):
% MI = MI_eval(S.rx_1sps,S.tx,DSP.DEMAPPER.C,DSP.DEMAPPER.N0,TX.QAM.symProb);
% [GMI,NGMI] = GMI_eval(S.rx_1sps,DSP.DEMAPPER.txBits,DSP.DEMAPPER.C,DSP.DEMAPPER.N0,TX.QAM.symProb);
% %Plot Constellation after Demapping:
% scatterPlot(S.rx_1sps,'EVM',EVM,'const',DSP.DEMAPPER.C);

%AWGN Theory
EVMt = 10^(-SNR_dB/20)*100;
[NGMIt,GMIt,MIt] = snr2gmi(SNR_dB,SIG.M,TX.QAM.entropy);

%Print Results
entranceMsg('RESULTS');
fprintf('\nEVM (theory): %1.3f%%',EVMt);
fprintf('\nEVM: %1.3f%%',EVM);
fprintf('\n\nMI (theory): %1.3f',MIt);
% fprintf('\nMI: %1.3f\n',MI);%�����޷���� ����
fprintf('\nGMI (theory): %1.3f\n',GMIt);
% fprintf('\nGMI: %1.3f\n',GMI);%�����޷���� ����
fprintf('\nNGMI (theory): %1.3f\n',NGMIt);
% fprintf('\nNGMI: %1.3f\n',NGMI);

%Exit PROG
exitProg();