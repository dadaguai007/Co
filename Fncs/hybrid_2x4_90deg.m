function Eo = hybrid_2x4_90deg(Es, Elo)
% Optical 2 x 4 90° hybrid.
%
% Parameters
% ----------
% Es : Input signal optical field.
% Elo : Input LO optical field.
% Eo : Optical hybrid outputs.
if ~isequal(length(Es), length(Elo))
    error('E1 and E2 need to have the same shape');
end
if ~isrow(Es)
    Es=Es.';
end
if ~isrow(Elo)
    Elo=Elo.';
end

% Optical hybrid transfer matrix
T = [
    1/2, 1i/2, 1i/2, -1/2;
    1i/2, -1/2, 1/2, 1i/2;
    1i/2, 1/2, -1i/2, -1/2;
    -1/2, 1i/2, -1/2, 1i/2
    ];

Ei = [Es; zeros(size(Es)); zeros(size(Es)); Elo];

Eo = T * Ei;
end
