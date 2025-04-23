function [out, esfreq] = FSE(in, Fs)
L = length(in);
% maxL = round((maxfo/Fs*L)*4)+1;
judge = abs(fft(angle(in)));
for i = 1:size(in, 2)
    idx = find(judge(:, i) == max(judge(:, i)));
    if idx(1) == 1
        judge(1) = 0;
        idx = find(judge(:, i) == max(judge(:, i)));
    end
    frac = (idx(1)-1)/4;
    esfreq(i) = frac * Fs / L;
    fshift = 2*pi*esfreq(i)/Fs*(0:L-1).';
    out(:, i) = in(:, i) .* exp(-1j*fshift);
    
    dirc_judge = abs(fft(angle(out)));
    dirc_idx = find(dirc_judge(:, i) == max(dirc_judge(:, i)));
    dirc_frac = (dirc_idx(1)-1)/4;
    dirc_esfreq = dirc_frac * Fs / L;
    
    if dirc_esfreq > esfreq(i)
        out(:, i) = in(:, i) .* exp(1j*fshift);
    end
end

end