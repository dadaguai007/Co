function phi = phaseNoise(lw, Nsamples, Ts)
    % Generate realization of a random-walk phase-noise process.
    % 
    % Parameters
    % ----------
    % lw : 
    %     Laser linewidth.
    % Nsamples : double
    %     Number of samples to be drawn.
    % Ts : 
    %     Sampling period.

    sigma2 = 2 * pi * lw * Ts;
    phi = zeros(1, Nsamples);

    for index = 1:(Nsamples - 1)
        phi(index + 1) = phi(index) + normrnd(0, sqrt(sigma2));
    end
end
