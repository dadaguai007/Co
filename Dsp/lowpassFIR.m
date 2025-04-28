function h = lowpassFIR(fc, fs, N, typeF)
    % Calculate FIR coefficients of a lowpass filter.
    % Parameters:
    %   fc: Cutoff frequency.
    %   fa: Sampling frequency.
    %   N: Number of filter coefficients.
    %   typeF: Type of response ('rect', 'gauss').
    
    %norm cut off frequence fu （截止频率除以采样频率）
    fu = fc / fs;
    % N should be odd number , so the d represents the index in the middle position 
    d = (N - 1) / 2;
    % 系数对称
    n = 0:N-1;

    % Calculate filter coefficients
    if strcmp(typeF, 'rect')
        h = (2 * fu) * sinc(2 * fu * (n - d));
    elseif strcmp(typeF, 'gauss')
        h = sqrt(2 * pi / log(2)) * fu * exp(-(2 / log(2)) * (pi * fu * (n - d)).^2);
    end
end
