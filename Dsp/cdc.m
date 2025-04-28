function Eo = cdc(Ei, param)
% Electronic chromatic dispersion compensation (CDC).
paramCh=struct();

if isfield(param, 'Fs')
    paramCh.Fs = param.Fs;
else
    error('The sampling rate is missing ')
end
if isfield(param, 'L')
    paramCh.L = param.L;
else
    paramCh.L = 50;
end
if isfield(param, 'D')
    paramCh.D = -param.D;
else
    paramCh.D = -16;
end
if isfield(param, 'Fc')
    paramCh.Fc = param.Fc;
else
    paramCh.Fc  = 193.1e12;
end

paramCh.alpha = 0;

Eo = linearChannel(Ei, paramCh);


end
