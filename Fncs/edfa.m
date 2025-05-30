function Eo = edfa(Ei, param)
% EDFA model.
% Default parameters
h = 6.62607004e-34;
Fs = param.Fs;
G = 20;          % Amplifier gain in dB
NF = 4.5;        % EDFA noise figure in dB
Fc = 193.1e12;   % Central optical frequency

%input parameters
if isfield(param, 'G')
    G = param.G;
end
if isfield(param, 'NF')
    NF = param.NF;
end
if isfield(param, 'Fc')
    Fc = param.Fc;
end
if isfield(param, 'type')
    type = param.type;
end
if G < 0
    disp('EDFA gain should be a positive scalar');
end

if NF <= 3
    disp('The minimal EDFA noise figure is 3 dB');
end


%lin mw
NF_lin = 10^(NF / 10);
G_lin = 10^(G / 10);
nsp = (G_lin * NF_lin - 1) / (2 * (G_lin - 1));

% ASE noise power calculation
N_ase = (G_lin - 1)*nsp*h*Fc;
p_noise = N_ase * Fs;

%noise
if real(Ei)
    noise = gaussianRealNoise(size(Ei), p_noise);
else
    noise = gaussianComplexNoise(size(Ei), p_noise);
end

if strcmp(type, 'noise')
    %output
    Eo = Ei * sqrt(G_lin) + noise;
else
    Eo = Ei * sqrt(G_lin);
end
end
