function [H, Sd, err_sq] = rlsUp(x, dx, outEq, lambda, H, Sd, nModes)
    % Coefficient update with the RLS algorithm.

    nTaps = size(H, 2);
    indMode = 1:nModes;
    indTaps = 1:nTaps;
    % the outEq is the 1×N
    outEq=outEq.';
    err = dx - outEq; % calculate output error for the RLS algorithm
    errDiag = diag(err); % define diagonal matrix from error array
    % update equalizer taps and inverse correlation matrix
    for N = 1:nModes
        indUpdModes = indMode + (N - 1) * nModes;
        indUpdTaps = indTaps + (N - 1) * nTaps;

        Sd_ = Sd(indUpdTaps, :);

        inAdapt = conj(x(:, N).'); % input samples
        inAdaptPar = repmat(inAdapt, nModes, 1); % expand input to parallelize tap adaptation
        
        % should exam the size ×
        Sd_ = (1 / lambda) * (Sd_ - (Sd_ * (inAdapt'* conj(inAdapt))* Sd_) / (lambda + conj(inAdapt) * Sd_ * inAdapt.'));

        H(indUpdModes, :) = H(indUpdModes, :) + errDiag * (Sd_ * inAdaptPar.').';

        Sd(indUpdTaps, :) = Sd_;
    end

    err_sq = abs(err).^2;

end
