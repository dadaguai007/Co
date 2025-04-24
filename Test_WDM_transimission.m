%Transmitter parameters:
clc;clear;close all;
M   = 4;          % order of the modulation format
constType = 'qam';  % constellation type
Rs  = 20e9;         % symbol rate [baud]
SpS = 16 ;           % samples per symbol
Nbits = (log2(M)*1e5);   % total number of bits per polarization
pulse = 'rrc';      % pulse shaping filter
Ntaps = 1024;       % number of pulse shaping filter coefficients
alphaRRC = 0.01;    % RRC rolloff
% the power is the input or output
Pch_dBm = 0;        % power per WDM channel [dBm]
Channel_power_type='output'; % input or output
Nch     = 5;        % number of WDM channels
Fc      = 193.1e12; % central optical frequency of the WDM spectrum
lw      = 100e3;    % laser linewidth
Nmodes = 2 ;        % number of signal modes [2 for polarization multiplexed signals]
freqSpac=50e9;
Tx_laser_type='ideal';
%generate WDM signal
% sigTx, symbTx_, paramTx = simpleWDMTx(paramTx);

% signal parameter
Fs = Rs*SpS;
Ts  = 1/Rs ;
Ta  = 1/Fs;

% pulse
psfRollOff=0.01;
psfLength=100;
psfShape='sqrt';

%IQM
paramIQ.Vpi=2;
paramIQ.VbI=-paramIQ.Vpi;
paramIQ.VbQ=-paramIQ.Vpi;
paramIQ.Vphi=paramIQ.Vpi/2;
%LO
Plo_dBm  = 10;
Plo =10.^(Plo_dBm/10)*1e-3 ;
paramLO=struct();
paramLO.P = Plo;
paramLO.lw = 100e3;          % laser linewidth
paramLO.RIN_var = 0;
paramLO.Fs = Fs;
%PD
paramPD=struct();
paramPD.B =Rs;
paramPD.R =0.95;
paramPD.type = 'ideal';
paramPD.Fs=Fs;
%ssfm
param=struct();
param.Ltotal = 20; %km
param.Lspan =20;
param.hz= 0.5;
param.alpha=0.2;
param.D = 16;
param.gamma = 1.3;
param.Fc = 193.1e12;
param.NF = 4.5;
param.amp='edfa';
param.Fs=Fs;
param.maxIter = 10;      % maximum number of convergence iterations per step
param.tol = 1e-5;       % error tolerance per step
param.nlprMethod = 'True'; % use adaptive step-size based o maximum nonlinear phase-shift
param.maxNlinPhaseRot = 2e-2; % maximum nonlinear phase-shift per step






% WDM each channel frequence
freqGrid = (-floor((Nch) / 2):floor((Nch-1 )/ 2)) * freqSpac;

% channel number
if mod(Nch, 2) == 0
    freqGrid = freqGrid + freqSpac / 2;
end

%power
if iscell(Pch_dBm)
    if length(Pch_dBm) == Nch
        Pch = (10.^ (cell2mat(Pch_dBm) / 10)) * 1e-3;
        % Optical signal power per WDM channel
    end
else
    % dBm to W
    Pch = (10 .^ (Pch_dBm / 10)) * 1e-3;
    Pch = Pch * ones(Nch,1);
end

rng(1);
% Time array
t = (0:((Nbits / log2(M))*SpS - 1));

%设计根升余弦脉冲成型滤波器
hsqrt = rcosdesign(psfRollOff,psfLength,SpS,psfShape);

% pulse shape
if strcmp(pulse, 'nrz')
    pulse = pulseShape('nrz', SpS);
elseif strcmp(pulse, 'rrc')
    pulse = pulseShape('rrc', SpS, 4096, psfRollOff, Ts);
end
pulse = pulse / max(abs(pulse));

% Allocate arrays
sigTxWDM = complex(zeros(length(t), Nmodes));
%symbol
symbTxWDM = complex(zeros(length(t) / SpS, Nmodes, Nch));
% total power of WDM signal，初始值为0。
Psig = 0;

for indCh = 1:Nch
    %通道数
    fprintf('channel %d\t fc : %3.4f THz\n', indCh, (Fc + freqGrid(indCh)) / 1e12);

    Pmode = 0;
    for indMode = 1:Nmodes
        %模式数,每个模式的输入功率
        fprintf(' mode %d\t power: %.2f dBm\n', indMode, 10 * log10((Pch(indCh)/Nmodes) / 1e-3));
        rng(10);
        % Generate random bits
        bitsTx = randi([0, 1], log2(M), Nbits/log2(M));

        % Map bits to constellation symbols
        symbTx=qammod(bitsTx,M,'InputType','bit','UnitAveragePower',1);
        symbTx = pnorm(symbTx);
        % Normalize symbols energy to 1
        symbTxWDM(:, indMode, indCh) = symbTx;

        % Upsampling
        symbolsUp = upsample(symbTx, SpS);
        % Pulse shaping
        sigTx=conv(symbolsUp,hsqrt,'same');
        %         sigTx = firFilter(pulse, symbolsUp);
        % Optical modulation
        if indMode == 1
            if strcmp(Tx_laser_type,'phasenoise')
                % Generate LO field with phase noise
                phi_pn_lo = phaseNoise(lw, length(sigTx), Ta);
                sigLO = exp(1i * phi_pn_lo);
                Pin=sigLO;
            else
                % this is important，理想的输入模式
                Pin=1;
            end
        end
        %调制
        Amp=0.5;
        if indMode==1
            sigTxCh = iqm(Pin, Amp*sigTx,paramIQ);
        else
            sigTxCh = DP_iqm(Pin, Amp*sigTx,paramIQ);
        end
        % 每个通道的发射功率是相对于所有偏振模式总功率，不是相对于单个偏振模式的功率。
        % 将发射功率除以偏振模式的数量可以得到每个通道在单个偏振模式下的功率。
        if strcmp(Channel_power_type,'output')
            % 这句的作用是，输出的光功率 就是 设置的功率，不用进行放大了，也即前面设置的功率为输出的光功率
            sigTxCh = sqrt(Pch(indCh) / Nmodes) * pnorm(sigTxCh);
        elseif strcmp(Channel_power_type,'input')
            % 这里的作用是， 前面设置的即为输入光功率为多少，
            sigTxCh = sqrt(Pch(indCh) / Nmodes)* (sigTxCh) ;
        end
        power=signalpower(sigTxCh);
        fprintf(' mode %d optical signal power: %.2f dBm\n',indMode, 10 * log10(power / 1e-3));
        %mode and channel sum
        phi_sig=sigTxCh.*exp(1i*2*pi*(freqGrid(indCh)/Fs)* t);
        sigTxWDM(:, indMode) = sigTxWDM(:, indMode) + phi_sig.';
        %total power is sum(abs(x.^2))
        Pmode = Pmode + signalpower(sigTxCh);
        %         fprintf('optical mode %d\t power: %.2f dBm\n', indMode, 10 * log10(signalpower(sigTxCh) / 1e-3));
    end

    Psig = Psig + Pmode;
    %每个通道的功率 dBm
    fprintf('channel %d\t power: %.2f dBm\n', indCh, 10 * log10(Pmode / 1e-3));
end
%整个WDM输出的功率 dBm
fprintf('total WDM signal power: %.2f dBm\n', 10 * log10(Psig / 1e-3));

%ssfm

sigRxoWDM = manakovssfm(sigTxWDM, param);
plot_spectrum(sigRxoWDM,Fs);
fprintf('total ssfm signal power model 1: %.2f dBm\n', 10 * log10(signalpower(sigRxoWDM(:,1))/ 1e-3));
fprintf('total ssfm signal power model 2: %.2f dBm\n', 10 * log10(signalpower(sigRxoWDM(:,2)) / 1e-3));

[~,t_up]=freq_time_set(Nbits/log2(M)*SpS,Fs);
% index
chIndex=1;

% LO
paramLO.N = length(sigTxWDM);
sigLO = basicLaserModel(paramLO);
FO      = 150e6;                 % frequency offset
f_lo   = freqGrid(chIndex)+FO;  % shift of the channel to be demodulated
sigLO = sigLO.*exp(1j*2*pi*f_lo*t_up); % add frequency offset

%PD
theta=0;
Sig_Rx=pdmCoherentReceiver(sigTxWDM, sigLO, theta, paramPD);
scatterplot(Sig_Rx(:,1))


%match filtering
for indMode = 1:Nmodes
Sig_Rx(:,indMode)=conv(Sig_Rx(:,indMode),hsqrt,'same');
Sig_Rx(:,indMode)=pnorm(Sig_Rx(:,indMode));
end
plot_spectrum(Sig_Rx,Fs);
% CD compensation
paramEDC = struct();
paramEDC.L = param.Ltotal;
paramEDC.D = param.D;
paramEDC.Fc = Fc-f_lo;
paramEDC.Fs = Fs;
sigRx = cdc(Sig_Rx, paramEDC);
scatterplot(Sig_Rx(:,1))
% 最佳采样
paramDec = struct();
paramDec.SpS_in  = SpS;
paramDec.SpS_out = 2;
%最佳采样相位,最佳采样点为delay+1
[SigRx,delay] = decimate(Sig_Rx, paramDec);

% 同步：
symbtx=symbTxWDM(:,:,chIndex);
% syns;
% for indMode = 1:Nmodes
% [sigrx(:,indMode),~,ff] = sync(SigRx(:,indMode),symbtx(:,indMode));
% end
symbRx = symbolSync(SigRx, symbtx, 2);

%norm
x = pnorm(SigRx);
d = pnorm(symbRx);
%均衡即可
% Equ
paramEq=struct();
paramEq.nTaps=5;
paramEq.numIter = 5;
% paramEq.mu =5e-3;
paramEq.mu = [5e-3, 2e-4];
%RLS forgetting factor
paramEq.lambdaRLS=0.99;
paramEq.SpS=paramDec.SpS_out;
% coefficient taps
paramEq.H=[];
% Eq length
paramEq.L = [floor(0.2*length(d)), floor(0.8*length(d))];
% paramEq.L = [750,800];
paramEq.Hiter=[];
paramEq.storeCoeff='true';
% alg is the cell
% paramEq.alg='da-rde';
paramEq.constSymb=pnorm(GrayMapping(M, 'qam'));
paramEq.alg = {'cma','cma'};
[yEq, H, errSq, Hiter] = mimoAdaptEqualizer(x, d, paramEq);
scatterplot(yEq(:,1))
scatterplot(yEq(:,2))



% Constellation
z = (0:M-1)';
y = qammod(z,M);
% CPR frequence offest
paramCPR = struct();
paramCPR.alg = 'bps';
paramCPR.N= 85;
paramCPR.B= 64;
paramCPR.pilotInd = (1:20:length(SigRx)); 
paramCPR.Ts=Ts;
paramCPR.constSymb=y;
% paramCPR.constSymb= GrayMapping(M, 'qam') ;
[y_CPR_BPS,theta] = cpr(yEq,d,paramCPR);
y_CPR_BPS = pnorm(y_CPR_BPS);
eyediagram(y_CPR_BPS(1:10000,2),12)
title('After BPS CPR ')
scatterplot(y_CPR_BPS(:,1))
% scatterplot(y_CPR_BPS(:,2))
figure;hold on;
plot(real(y_CPR_BPS(:,1)))
plot(real(y_CPR_BPS(:,2)))