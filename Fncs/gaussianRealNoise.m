function noise=gaussianRealNoise(shapeOut, sigma2)
% use the size to ensure the shapeOut
% Generate real circular Gaussian noise.
    noise = sqrt(sigma2) * (randn(shapeOut));
end
