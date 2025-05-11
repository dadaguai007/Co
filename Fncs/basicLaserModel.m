function optical_signal = basicLaserModel(param)
% Laser model with Maxwellian random walk phase noise and RIN.

% Input arguments:
% param: Parameters of the laser (struct).
%   - param.P: laser power [W] [default: 10 dBm]
%   - param.lw: laser linewidth [Hz] [default: 1 kHz]
%   - param.RIN_var: variance of the RIN noise [default: 1e-20]
%   - param.Fs: sampling rate [samples/s]
%   - param.Ns: number of signal samples [default: 1e3]

% Output:
% optical_signal: Optical signal with phase noise and RIN.

lw=1e3;
RIN_var = 1e-20;
if isfield(param, 'Fs')
    Fs = param.Fs;
end
if isfield(param, 'P')
    P = param.P; % Laser power in W
end

if isfield(param,'lw')
    lw=param.lw;
end
if isfield(param,'RIN_var')
    RIN_var = param.RIN_var;
end
if isfield(param,'N')
    N=param.N;
end

t = (0:N-1) * 1 / Fs;

% Simulate Gaussion random walk phase noise
pn = phaseNoise(lw, N, 1/Fs);

% Simulate relative intensity noise (RIN) 
% may be note the real signal noise
deltaP = gaussianComplexNoise(size(pn), RIN_var);

% Return optical signal
% this optical lasear is not define the frequence offest
% 2*pi*f_lo*t 
% this should be the  input of the modulator
optical_signal = sqrt(P) * exp(1i * pn) + deltaP;

end
