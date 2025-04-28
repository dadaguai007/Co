function phaseError = viterbi(Ei, N, M)
    % Viterbi & Viterbi carrier phase recovery algorithm.

    if nargin < 2
        N = 35;
    end
    if nargin < 3
        M = 4;
    end

    phaseError = -unwrap(angle(movingAverage(Ei.^M, N)) / M, 2*pi / M) - pi / 4;
end
