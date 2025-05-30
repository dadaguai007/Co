function ipd = pd(E, param)
% Pin photodiode (PD).
%
% Parameters
% ----------
% E : Input optical field.
% param : parameter object (struct)

if nargin < 2
    param = struct();
end

% Parameters default
q  = 1.60217663e-19; % Elementary charge in coulombs
R=1;
Id=5e-9;
B=30e9;

kB  = 1.380649e-23; % Boltzmann constant (in J/K)
Tc=25;
RL = 50;           % RL in Ohms
type = 'ideal';
% I_max = 5e3; %饱和电流
N = 8001;
fType = 'rect';


% Check input parameters
if isfield(param, 'R')
    R = param.R;
end
if isfield(param, 'Id')
    Id = param.Id;
end
if isfield(param, 'B')
    B = param.B;
end
if isfield(param, 'Tc')
    Tc = param.Tc;
end
if isfield(param, 'RL')
    RL = param.RL;
end
% if isfield(param, 'I_max')
%     I_max = param.I_max;
% end
if isfield(param, 'type')
    type = param.type;
end
if isfield(param, 'N')
    N = param.N;
end
if isfield(param, 'fType')
    fType = param.fType;
end
if isfield(param, 'Fs')
    Fs = param.Fs;
end

% Ideal photocurrent
ipd = R * E .* conj(E);

% this can be modify
% 判断 type 是否为 ideal

if strcmp(type, 'ideal')
    ipd = R * E .* conj(E);  % Ideal photocurrent
%     ipd1 = R * abs(E).^2;
else
    %  ipd(ipd > I_max) = Ipd_sat;
    ipd_mean = mean(ipd);

    % Shot noise
    sigma_s = 2 * q * (ipd_mean + Id) * B;  % Shot noise variance

    % Thermal noise
    T = Tc + 273.15;  % Temperature in Kelvin
    sigma_T = 4 * kB * T * B / RL;  % Thermal noise variance

    % Add noise sources to the PIN receiver
    Is = normrnd(0, sqrt(Fs * (sigma_s / (2 * B))), size(ipd));
    It = normrnd(0, sqrt(Fs * (sigma_T / (2 * B))), size(ipd));

    ipd = ipd + Is + It;

    % Lowpass filtering
    h = lowpassFIR(B, Fs, N, fType);
    ipd = firFilter(h, ipd);
end
ipd = real(ipd);
end
