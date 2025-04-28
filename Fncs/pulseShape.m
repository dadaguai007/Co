function filterCoeffs = pulseShape(pulseType, sps, N, alpha, Ts)
    % Generate a pulse shaping filter.
    if nargin < 2
        sps = 2;
    end
    if nargin < 3
        N = 1024;
    end
    if nargin < 4
        alpha = 0.1;
    end
    if nargin < 5
        Ts = 1;
    end
    
    fs = 1/Ts* sps;

    if strcmp(pulseType, 'rect')
        filterCoeffs = [zeros(1, sps/2), ones(1, sps), zeros(1, sps/2)];
    elseif strcmp(pulseType, 'nrz')
        t = linspace(-2, 2, sps);
        Te = 1;
        filterCoeffs = conv(ones(1,sps), 2/(sqrt(pi)*Te)*exp(-t.^2/Te),"full");
    elseif strcmp(pulseType, 'rrc')
        t = linspace(-N/2, N/2, N) * (1/fs);
        filterCoeffs = rrcFilterTaps(t, alpha, Ts);
    elseif strcmp(pulseType, 'rc')
        t = linspace(-N/2, N/2, N) * (1/fs);
        filterCoeffs = rcFilterTaps(t, alpha, Ts);
    end
    filterCoeffs = filterCoeffs / sqrt(sum(filterCoeffs.^2));
end
