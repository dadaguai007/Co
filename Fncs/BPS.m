function [wout, estimatedPhase] = BPS(win, constellation, ML, steptimes)
if nargin < 4
    steptimes = 2;
end
win = win ./ sqrt(bandpower(win));
constellation = constellation / sqrt(bandpower(constellation));
[L, N] = size(win);
B = length(constellation) * steptimes;
phaseshift = (-B/2:B/2-1)/B * pi / 2;
estimatedPhase = zeros(L, N);
wout = zeros(L, N);

for i = 1:N
    temp = win(:, i) * exp(1j*phaseshift);
    judge = decision(temp, constellation);
    d = temp - judge;
    d2 = d .* conj(d);
    block = zeros(2*ML+1, B);
    block(1+ML:end, :) = d2(1:1+ML, :);
    for j = 1:ML
        S(j, :) = sum(block, 1);
        block = [block(2:end, :); d2(j+ML+1, :)];
    end
    for j = ML+1:L-ML
        block = d2(j-ML:j+ML, :);
        S(j, :) = sum(block, 1);
    end
    for j = L-ML+1:L
        block = [block(2:end, :); zeros(1, B)];
        S(j, :) = sum(block, 1);
    end
    [~, index] = min(transpose(S));
    estimatedPhase(:, i) = phaseshift(index).';
    wout(:, i) = win(:, i) .* exp(1j*estimatedPhase(:, i));
end

end