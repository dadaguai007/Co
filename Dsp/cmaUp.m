function [H, err_sq] = cmaUp(x, R, outEq, mu, H, nModes)
    % Coefficient update with the CMA algorithm.

    indMode = 1:nModes;
    % out onces transpose to   1 × modes
    outEq = outEq.';
    err = R - abs(outEq).^2;  % calculate output error for the CMA algorithm

    prodErrOut = diag(err) * diag(outEq);  % define diagonal matrix

    % update equalizer taps coffcient
    for N = 1:nModes
        indUpdTaps = indMode + (N - 1) * nModes;  % simplify indexing
        inAdapt = x(:, N).';  % input samples
        inAdaptPar = repmat(inAdapt, nModes, 1);  % expand input to parallelize tap adaptation
        H(indUpdTaps, :) = H(indUpdTaps, :) + mu * prodErrOut * conj(inAdaptPar);  % gradient descent update
    end

    err_sq = abs(err).^2;

end
