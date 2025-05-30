function Ao = iqm(Ai, u, param)
% Optical In-Phase/Quadrature Modulator (IQM).

% Default parameter values
Vpi = 2;
VbI = -2;
VbQ = -2;
Vphi = 1;

% Check if 'param' is provided and update parameters if available
if nargin < 3
    param = struct();
end

if isfield(param, 'Vpi')
    Vpi = param.Vpi;
end

if isfield(param, 'VbI')
    VbI = param.VbI;
end

if isfield(param, 'VbQ')
    VbQ = param.VbQ;
end

if isfield(param, 'Vphi')
    Vphi = param.Vphi;
end

% Ensure 'u' is a complex-valued array
if ~isnumeric(u)
    u = complex(u);
end

% Ensure 'Ai' and 'u' have the same dimensions
if isscalar(Ai)
    Ai = Ai * ones(size(u));
end

% Define parameters for the I-MZM and Q-MZM
paramI.Vpi = Vpi;
paramI.Vb = VbI;

paramQ.Vpi = Vpi;
paramQ.Vb = VbQ;

% Perform IQM modulation
outputI = mzm(Ai / sqrt(2), real(u), paramI.Vb,paramI.Vpi );
outputQ = mzm(Ai / sqrt(2), imag(u), paramQ.Vb,paramQ.Vpi);
outputPM = pm(outputQ, Vphi * ones(size(u)), Vpi);

% Combine the I-MZM and PM modulation
Ao = outputI + outputPM;
end
