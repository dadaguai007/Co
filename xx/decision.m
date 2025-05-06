function out = decision(in, constellation)
    out = zeros(size(in));
    % in = normalize_ichi(in, 'meanpwr');
    % constellation_ref = normalize_ichi(constellation, 'meanpwr');
    for idx = 1:size(in, 2)
        distance = abs(in(:, idx) - constellation(:).');
        [~, index] = min(distance');
        out(:, idx) = constellation(index);
    end
end