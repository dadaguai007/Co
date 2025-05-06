function y=pnorm(x)
%normalizes the average power of each component of the input signal x. 
% It divides each component's power by the square root of the average power of the signal, 
% and then returns the normalized signal.
    y= x ./ sqrt(mean(real(x .* conj(x))));
end