function noise = gaussianComplexNoise(shapeOut, sigma2)
    % Generate complex circular Gaussian noise.

    if nargin < 2
        sigma2 = 1.0;
    end

    noise = sqrt(sigma2/2) * (randn(shapeOut) + 1i * randn(shapeOut));
end
