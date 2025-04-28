function y = clockSamplingInterp(x, Fs_in, Fs_out, jitter_rms)
% ADC中添加频偏和抖动  
% Interpolate signal to a given sampling rate.

    % Input arguments:
    % x: Input signal.
    % Fs_in: Sampling frequency of the input signal.
    % Fs_out: Sampling frequency of the output signal.
    % jitter_rms: Standard deviation of the time jitter.

    % Output:
    % y: Resampled signal.
    
    % the input x should be N*nmodes
    if size(x,1)<size(x,2)
        x=x';
    end
    nModes = size(x, 2);
    % 输入时间
    inTs = 1 / Fs_in;
    % 输出时间，偏离
    outTs = 1 / Fs_out;
   
    tin = (0:length(x)-1) * inTs;

    tout = 0:outTs:(length(x)-1)*inTs;
    L=length(tout);
    jitter=jitter_rms.*randn(1,L);
    %
    tout=tout+jitter;


    y = zeros(length(tout), size(x, 2));
    % 使用插值算法进行数据插取
    for k = 1:nModes
        y(:, k) = interp1(tin, x(:, k), tout,"spline");
    end
end
