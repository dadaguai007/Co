function Eo = linearChannel(Ei, param)

Fs = param.Fs;
L = param.L;
alpha = param.alpha;
D = param.D;
Fc = param.Fc;

c = 299792458; % speed of light [m/s] (vacuum)
c_kms = c / 1e3;
lambda = c_kms / Fc;
alpha = alpha / (10 * log(exp(1)));
beta2 = -(D * lambda^2) / (2 * pi * c_kms);

Nfft = length(Ei);

omega = 2 * pi * Fs * (-0.5:1/Nfft:0.5-1/Nfft);
omega=ifftshift(omega);

% transport to  N × 1
omega = reshape(omega, [length(omega), 1]);


% make Ei to N × （1 or 2）
[r,l]=size(Ei);
if r<l
    %row < lin , to invert the martix
  Ei= Ei.';
end

% check model N × （1 or 2）
Nmodes = size(Ei, 2);
Eo=zeros(size(Ei));
for n=1:Nmodes
    Eo(:,n) = ifft(fft(Ei(:,n)) .* exp( -alpha / 2 * L - 1i * (beta2 / 2) * (omega.^2) * L));
end

% N × 1
if Nmodes == 1
    Eo = reshape(Eo, [length(Eo), 1]);
end

% N × 1 to repmat N × 2 or N × 1
% omega = repmat(omega, [1, Nmodes]);
% Eo = ifft(fft(Ei) .* exp(-alpha / 2 * L + 1i * (beta2 / 2) * (omega.^2) * L));

end
