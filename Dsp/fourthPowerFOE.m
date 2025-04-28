function [Eo, fo] = fourthPowerFOE(Ei, Fs)
    % Estimate the frequency offset (FO) with the 4th-power method.
    % The Ei must to be the N × (1,2)
    %row
    Nfft = size(Ei, 1);

    f = Fs *(-0.5:1/Nfft:0.5-1/Nfft);
    %line
    nModes = size(Ei, 2);
    
    Eo = zeros(size(Ei));
    % time array
    t = (0:(Nfft-1)) * (1/Fs);

    for n = 1:nModes
        f4 = 10 * log10(abs(fftshift(fft(Ei(:, n).^ 4))));
        [~, indFO] = max(f4);
        fo = f(indFO) / 4;
        phi=exp(-1i * 2 * pi * fo * t).';
        Eo(:, n) = Ei(:, n) .* phi;
    end

end


