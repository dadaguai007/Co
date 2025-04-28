function Eo = adc(Ei, param)
    % Analog-to-digital converter (ADC) model.

    % Input arguments:
    % Ei: Input signal.
    % param: Resampling parameters (struct).
    %   - param.Fs_in  : sampling frequency of the input signal [default: 1 sample/s]
    %   - param.Fs_out : sampling frequency of the output signal [default: 1 sample/s]
    %   - param.jitter_rms : root mean square (RMS) value of the jitter in seconds [default: 0 s]
    %   - param.nBits : number of bits used for quantization [default: 8 bits]
    %   - param.Vmax : maximum value for the ADC's full-scale range [default: 1V]
    %   - param.Vmin : minimum value for the ADC's full-scale range [default: -1V]
    %   - param.AAF : flag indicating whether to use anti-aliasing filters [default: True]
    %   - param.N : number of taps of the anti-aliasing filters [default: 201]

    % Output:
    % Eo: Resampled and quantized signal.

    % Check and set default values for input parameters
%     param.Fs_in = getfield(param, 'Fs_in', 1);
%     param.Fs_out = getfield(param, 'Fs_out', 1);
%     param.jitter_rms = getfield(param, 'jitter_rms', 0);
%     param.nBits = getfield(param, 'nBits', 8);
%     param.Vmax = getfield(param, 'Vmax', 1);
%     param.Vmin = getfield(param, 'Vmin', -1);
%     param.AAF = getfield(param, 'AAF', true);
%     param.N = getfield(param, 'N', 201);

    % Extract individual parameters for ease of use
    Fs_in = param.Fs_in;
    Fs_out = param.Fs_out;
    jitter_rms = param.jitter_rms;
    nBits = param.nBits;
    Vmax = param.Vmax;
    Vmin = param.Vmin;
    AAF = param.AAF;
    % 混叠滤波器的抽头数
    N = param.N;

    % Reshape the input signal if needed to handle single-dimensional inputs
    if size(Ei, 1) < size(Ei,2)
        Ei = Ei.';
    end

    % Get the number of modes (columns) in the input signal
%     nModes = size(Ei, 2);

    % Apply anti-aliasing filters if AAF is enabled
    if strcmp(AAF,'on')
        % Anti-aliasing filters:
        Ntaps = min(size(Ei, 1), N);
        hi = lowpassFIR(Fs_out / 2, Fs_in, Ntaps, 'rect');
        ho = lowpassFIR(Fs_out / 2, Fs_out, Ntaps, 'rect');

        Ei = firFilter(hi, Ei);
    end

    if ~isreal(Ei)
        % Signal interpolation to the ADC's sampling frequency
        Eo = clockSamplingInterp(real(Ei), Fs_in, Fs_out, jitter_rms) + 1i * clockSamplingInterp(imag(Ei), Fs_in, Fs_out, jitter_rms);

        % Uniform quantization of the signal according to the number of bits of the ADC
        Eo = quantizer(real(Eo), nBits, Vmax, Vmin)+1i*quantizer(imag(Eo), nBits, Vmax, Vmin);
    else
        % Signal interpolation to the ADC's sampling frequency
        Eo = clockSamplingInterp(Ei, Fs_in, Fs_out, jitter_rms);

        % Uniform quantization of the signal according to the number of bits of the ADC
        Eo = quantizer(Eo, nBits, Vmax, Vmin);
    end

    % Apply anti-aliasing filters to the output if AAF is enabled
    if strcmp(AAF,'on')
        Eo = firFilter(ho, Eo);
    end
end
