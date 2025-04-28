function  Eo=ssfm(Ei,param)
% 想要加上EDFA的噪声干扰，需要更改该函数的EDFA选项
% Default parameters
Ltotal = 400; %km
Lspan = 80;
hz= 0.5;
alpha=0.2;
D = 16;
gamma = 1.3;
Fc = 193.1e12;
NF = 4.5;
amp='edfa';

%Ei should be row
if ~isrow(Ei)
    Ei=Ei.';
end

if isfield(param, 'Fs')
    Fs = param.Fs;
end
if isfield(param, 'Ltotal')
    Ltotal = param.Ltotal;
end
if isfield(param, 'Lspan')
    Lspan = param.Lspan;
end
if isfield(param, 'hz')
    hz = param.hz;
end
if isfield(param, 'alpha')
    alpha = param.alpha;
end
if isfield(param, 'D')
    D = param.D;
end
if isfield(param, 'gamma')
    gamma = param.gamma;
end
if isfield(param, 'Fc')
    Fc = param.Fc;
end
if isfield(param, 'NF')
    NF = param.NF;
end
if isfield(param, 'amp')
    amp = param.amp;
end

c = 299792458;
c_kms = c/ 1e3 ;   % speed of light (vacuum) in km/s
lamba = c_kms / Fc ;
% lamba=1.550e-9;
Alpha = alpha / (10 * log10(exp(1)));
% Alpha =alpha;
beta2 = -(D * lamba.^2) / (2 * pi * c_kms);

%Amp
paramAmp = struct();
paramAmp.G = alpha * Lspan;
paramAmp.NF = NF;
paramAmp.Fc = Fc;
paramAmp.Fs = Fs;
paramAmp.type='noise'; % 是否加噪声决定于此
%step
Nspans = floor(Ltotal / Lspan);
Nsteps = floor(Lspan / hz);
%length
Nfft = length(Ei);
%omega
w = 2 * pi * Fs * (-0.5:1/Nfft:0.5-1/Nfft);
w=ifftshift(w);
% linear
linOperator = exp(-(Alpha/2 ) * (hz / 2) - 1j * (beta2/2) * (w.^2) * (hz / 2));
Ech=Ei;
% cycle
for i = 1:Nspans
    Ech = fft(Ech);
    for j= 1:Nsteps
        %  linear
        Ech = Ech.*linOperator;
        %  non linear
        Ech = ifft(Ech);
        Ech = Ech .* exp(-1i * gamma * (Ech.*conj(Ech)) * hz);
        %  linear
        Ech = fft(Ech);
        Ech = Ech.* linOperator;
    end
    Ech=ifft(Ech);
    % amplification step
    if strcmp(amp, 'edfa')
        Ech = edfa(Ech, paramAmp);
    elseif strcmp(amp, 'ideal')
        Ech = Ech .* exp(Alpha / 2 * Nsteps * hz);
    elseif strcmp(amp, 'none')
        Ech = Ech .* exp(0);
    end

end

%output
Eo=Ech;
end
