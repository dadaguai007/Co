function [H, err_sq] = ddlmsUp(x, constSymb, outEq, mu, H, nModes)
    % Coefficient update with the DD-LMS algorithm.

    indMode = 1:nModes;
    % the outEq is the 1×N
    outEq=outEq.';
    decided = zeros(size(outEq));

    for k = 1:nModes
        [~, indSymb] = min(abs(outEq(:, k) - constSymb));
        decided(:, k) = constSymb(indSymb);
    end

    err = decided - outEq;  % calculate output error for the DD-LMS algorithm

    errDiag = diag(err);  % define diagonal matrix from error array

    % update equalizer taps
    for N = 1:nModes
        indUpdTaps = indMode + (N - 1) * nModes;  % simplify indexing
        inAdapt = x(:, N).';  % input samples
        inAdaptPar = repmat(inAdapt, nModes, 1);  % expand input to parallelize tap adaptation
        H(indUpdTaps, :) = H(indUpdTaps, :) + mu * errDiag * conj(inAdaptPar);  % gradient descent update
    end

    err_sq = abs(err).^2;

end
