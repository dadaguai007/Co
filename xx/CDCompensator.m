classdef CDCompensator < unit 

    properties
        nInputs = 1;
        nOutputs = 1;
        L;
        D;
        S;
        lambda;
    end

    methods
        function obj = CDCompensator(param)
            if isfield(param, 'L') && isscalar(param.L), obj.L = param.L; else, slog('Length of the optical fiber (L) must be entered.', 'ERR'); end 
            if isfield(param, 'D') && isscalar(param.D), obj.D = param.D; else, obj.D = 16e-6; end 
            if isfield(param, 'S') && isscalar(param.S), obj.S = param.S; else, obj.S = 0.08e3; end 
            if isfield(param, 'lambda') && isscalar(param.lambda), obj.lambda = param.lambda; else, obj.lambda = 1550e-9; end 
        end

        function out = generate(obj, in)
            f = (-round(in.L/2)+1:(in.L-round(in.L/2)))*in.Fs/in.L;
            w = 2*pi*f';
            
            beta2 = -obj.D * obj.lambda^2 / (2*pi*const.c);   
            beta3 = (obj.lambda / (2*pi*const.c))^2 * (obj.lambda^2*obj.S+2*obj.lambda*obj.D);

            op = beta2 / factorial(2)*w.^(2) + beta3/factorial(3)*w.^(3);
            DTF = exp(1j*op*obj.L);
            out = fun1(in, @(x) ifft(fft(x).*DTF));
        end

    end

end