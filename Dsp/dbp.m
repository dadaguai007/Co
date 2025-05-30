function  Eo=dbp(Ei,param)


Ltotal = 400; %km
Lspan = 80;
hz= 0.5;
alpha=0.2;
D = 16;
gamma = 1.3;
Fc = 193.1e12;


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

c = 299792458;
c_kms = c/ 1e3 ;   % speed of light (vacuum) in km/s
lamba = c_kms / Fc ;
% 取反
Alpha = -alpha / (10 * log10(exp(1)));
beta2 = (D * lamba^2) / (2 * pi * c_kms);


%step
Nspans = floor(Ltotal / Lspan);
Nsteps = floor(Lspan / hz);
%length
Nfft = length(E);
%omega
w = 2 * pi * Fs * (-0.5:1/Nfft:0.5-1/Nfft);
% linear
linOperator = exp(-(Alpha / 2) * (hz / 2) + 1j * (beta2 / 2) * (w.^2) * (hz / 2));
Ech=Ei;
% cycle
for i = 1:Nspans
    Ech = fft(Ech);
    for j= 1:Nsteps
        %  linear
        Ech = Ech.* linOperator;
        %  non linear
        Ech = ifft(Ech);
        Ech = Ech .* exp(1i * gamma * (abs(Ech).^2) * hz);
        %  linear
        Ech = fft(Ech);
        Ech = Ech.* linOperator;
    end
    Ech=ifft(Ech);
end
%output
Eo=Ech;
end
