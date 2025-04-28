function [yEq, H, errSq, Hiter] = coreAdaptEq(x, dx, SpS, H, L, mu, lambdaRLS, nTaps, storeCoeff, alg, constSymb)
    % Adaptive equalizer core processing function

    % the data should be the N × (1,2)
    nModes = size(x, 2);
    indTaps = 1:nTaps;
    indMode = 1:nModes;

    errSq = zeros(nModes, L);
    yEq = zeros(L, nModes);
    outEq = zeros(nModes, 1);

    %record the  coefficient of the taps value of all signal symbols
    % 历史的抽头系数
    if strcmp(storeCoeff,'true')
        Hiter = zeros(nModes.^2, nTaps, L);
    else
        Hiter = zeros(nModes.^2, nTaps, 1);
    end

    if strcmp(alg, 'rls')
        Sd = eye(nTaps);
        for i = 1:nTaps - 1
            Sd = [Sd; eye(nTaps)];
        end
    end

    % Radii cma, rde
    %complex
    Rcma = complex(mean(abs(constSymb).^4) / mean(abs(constSymb).^2) * ones(1, nModes) );

    Rrde = unique(abs(constSymb));

    for ind = 1:L
        outEq(:) = 0;
        % 每次丢弃一个符号（多个采样点）,number is the taps
        indIn = indTaps + (ind - 1) * SpS;

        % pass signal sequence through the equalizer
        for N = 1:nModes
            inEq = x(indIn, N);  % slice input coming from the Nth mode
            % every model corresponding two row, and results get the 2×1
            % number ,so the out is also the 2×1 vector, it is can be add
            % together
            outEq = outEq + H(indMode + (N - 1) * nModes, :) * inEq;  % add contribution from the Nth mode
        end
        %make out Transpose, so get 1 × model vector ， every column is one model， 
        % so we can use the row vector replace the Eq out.
        yEq(ind, :) = outEq.';

        % update equalizer taps according to the specified algorithm and save squared error
        if strcmp(alg, 'nlms')
            [H, errSq(:, ind)] = nlmsUp(x(indIn, :), dx(ind, :), outEq, mu, H, nModes);
        elseif strcmp(alg, 'cma')
            [H, errSq(:, ind)] = cmaUp(x(indIn, :), Rcma, outEq, mu, H, nModes);
        elseif strcmp(alg, 'dd-lms')
            [H, errSq(:, ind)] = ddlmsUp(x(indIn, :), constSymb, outEq, mu, H, nModes);
        elseif strcmp(alg, 'rde')
            [H, errSq(:, ind)] = rdeUp(x(indIn, :), Rrde, outEq, mu, H, nModes);
        elseif strcmp(alg, 'da-rde')
            [H, errSq(:, ind)] = dardeUp(x(indIn, :), dx(ind, :), outEq, mu, H, nModes);
        elseif strcmp(alg, 'rls')
            [H, Sd, errSq(:, ind)] = rlsUp(x(indIn, :), dx(ind, :), outEq, lambdaRLS, H, Sd, nModes);
        elseif strcmp(alg, 'dd-rls')
            [H, Sd, errSq(:, ind)] = ddrlsUp(x(indIn, :), constSymb, outEq, lambdaRLS, H, Sd, nModes);
        elseif strcmp(alg, 'static')
            errSq(:, ind) = errSq(:, ind - 1);
        else
            error('Equalization algorithm not specified (or incorrectly specified).');
        end
        
        % record
        if storeCoeff
            Hiter(:, :, ind) = H;
        else
            Hiter(:, :, 1) = H;
        end
    end
end
