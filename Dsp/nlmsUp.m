function [H, err_sq] = nlmsUp(x, dx, outEq, mu, H, nModes)
    % Coefficient update with the NLMS algorithm.

    indMode = 1:nModes;
    % the outEq is the 1×N
    outEq = outEq.';
    err = dx - outEq; % calculate output error for NLMS algorithm

    errDiag = diag(err); % define diagonal matrix from error array

    % update equalizer taps
    for N = 1:nModes
        indUpdTaps = indMode + (N - 1) * nModes; % simplify indexing and improve speed
        inAdapt = x(:, N).'/ norm(x(:, N)).^2; % NLMS normalization
        inAdaptPar = repmat(inAdapt, nModes,1); % expand input to parallelize tap adaptation
        H(indUpdTaps, :) = H(indUpdTaps, :) + mu * errDiag * conj(inAdaptPar); % gradient descent update
    end

    err_sq = abs(err).^2;

end
