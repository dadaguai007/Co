function s = coherentReceiver(Es, Elo,param)
% Single polarization coherent optical front-end.
%
% Parameters
% ----------
% Es : Input signal optical field.
% Elo :  Input LO optical field.

% use row to calualator 
if ~isequal(size(Es), size(Elo))
    error('Es and Elo need to have the same shape');
end

% optical 2 x 4 90° hybrid
Eo = hybrid_2x4_90deg(Es, Elo);

% balanced photodetection
sI = bpd(Eo(2, :), Eo(1, :),param);
sQ = bpd(Eo(3, :), Eo(4, :),param);

s = sI + 1i * sQ;
end
