function [out, estimatedPhase] = V_V(in, ML)

[L, N] = size(in);
d = in .* conj(circshift(in, 1));
d4 = d.^4;

block = zeros(2*ML+1, N);
block(1+ML:end, :) = d4(1:1+ML, :);
for j = 1:ML
    S(j, :) = mean(block, 1);
    block = [block(2:end, :); d4(j+ML+1, :)];
end
for j = ML+1:L-ML
    block = d4(j-ML:j+ML, :);
    S(j, :) = mean(block, 1);
end
for j = L-ML+1:L
    block = [block(2:end, :); zeros(1, N)];
    S(j, :) = mean(block, 1);
end
estimatedPhase = zeros(L+1, N);
estimatedPhase(2:end, :) = cumsum(angle(S)/4);
estimatedPhase = estimatedPhase(1:end-1, :);
out = in .* exp(-1j*estimatedPhase);
end