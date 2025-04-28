function Theta = bps(Ei, N, constSymb, B)
    % Blind phase search (BPS) algorithm
    % N windows length/2
    % B number of test phi
    % Ei complex singal
    % constellation symbol

    % % The Ei must to be the N × (1,2)
    % Get the number of modes from the size of Ei
    nModes = size(Ei, 2);

    % Generate test phases
    phi_test = (0:(B-1)) * (pi/2) / B;

    % Initialize the theta array
    theta = zeros(size(Ei));

    % Create zero padding for the signal
    zeroPad = zeros(N, nModes);
    x = [zeroPad; Ei; zeroPad];  % Pad the signal with zeros

    %row
    L = size(x, 1);
    Theta=zeros(size(Ei));
    for n = 1:nModes

        dist = zeros(B, size(constSymb, 1));
        dmin = zeros(B, 2 * N + 1);

        for k = 1:L
            for indPhase = 1:B
                phi = phi_test(indPhase);
                %distance between test and constellation
                dist(indPhase, :) = abs(x(k, n) * exp(1i * phi) - constSymb).^2;
                dmin(indPhase, end) = min(dist(indPhase, :));
            end
            % get the sliding window length
            if k >= 2 * N
                % sum row
                sumDmin = sum(dmin, 2);
                [~, indRot] = min(sumDmin);
                % get the best theata
                theta(k-2 * N+1, n) = phi_test(indRot);
            end
              %shift symbol
            dmin = circshift(dmin, [0, -1]);
        end
        Theta(:,n)=theta(1:end-1,n);
    end
end
