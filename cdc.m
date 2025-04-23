function out = cdc(in, Len, D, S, lambda, Fs)
L = length(in);
c = 299792458;
f = (-round(L/2)+1:(L-round(L/2)))*Fs/L;
w = 2*pi*f';

beta2 = -D * lambda^2 / (2*pi*c);   
beta3 = (lambda / (2*pi*c))^2 * (lambda^2*S+2*lambda*D);

op = beta2 / factorial(2)*w.^(2) + beta3/factorial(3)*w.^(3);
DTF = exp(1j*op*Len);
out = ifft(fft(in).*DTF);
end