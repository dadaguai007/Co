function y = quantizer(x, nBits, maxV, minV)
    % Quantize the input signal using a uniform quantizer with the specified precision.

    % Input arguments:
    % x: The input signal to be quantized.
    % nBits: Number of bits used for quantization. The quantizer will have 2^nBits levels.
    % maxV: Maximum value for the quantizer's full-scale range (default is 1).
    % minV: Minimum value for the quantizer's full-scale range (default is -1).

    % Output:
    % y: The quantized output signal with the same shape as 'x', quantized using 'nBits' levels.
    if size(x,1)<size(x,2)
        x=x';
    end
    detal = (maxV - minV) / (2^nBits - 1);

    d = minV:detal:maxV;

    y = zeros(size(x));

    for indMode = 1:size(x, 2)
        for idx = 1:length(x)
%             Find the nearest value
            [~, minIdx] = min(abs(x(idx, indMode) - d));
            % make the nearest value of d
            y(idx, indMode) = d(minIdx);
        end
    end
end
