function [H, err_sq] = rdeUp(x, R, outEq, mu, H, nModes)
    % Coefficient update with the RDE algorithm.

    indMode = 1:nModes;
    % the outEq is the 1×N
    outEq=outEq.';
    decidedR = zeros(size(outEq));

    % find closest constellation radius
    for k = 1:nModes
        [~, indR] = min(abs(R - abs(outEq(:, k))));
        decidedR(:, k) = R(indR);
    end

    err = decidedR.^2 - abs(outEq).^2;  % calculate output error for the RDE algorithm

    prodErrOut = diag(err) * diag(outEq);  % define diagonal matrix

    % update equalizer taps
    for N = 1:nModes
        indUpdTaps = indMode + (N - 1) * nModes;  % simplify indexing
        inAdapt = x(:, N).';  % input samples
        inAdaptPar = repmat(inAdapt, nModes, 1);  % expand input to parallelize tap adaptation
        H(indUpdTaps, :) = H(indUpdTaps, :) + mu * prodErrOut * conj(inAdaptPar);  % gradient descent update
    end

    err_sq = abs(err).^2;

end
