function S = pdmCoherentReceiver(Es, Elo, theta_sig, param)
% Polarization multiplexed coherent optical front-end.

if ~isequal(length(Es), length(Elo))
    error('Es and Elo need to have the same shape');
end
[Elox,Eloy] = pbs(Elo, pi/4);  % Split LO into two orthogonal polarizations
[Esx, Esy] = pbs(Es, theta_sig);  % Split signal into two orthogonal polarizations

Sx = coherentReceiver(Esx, Elox, param);  % Coherent detection of pol.X
Sy = coherentReceiver(Esy, Eloy, param);  % Coherent detection of pol.Y
Sx=Sx.';
Sy=Sy.';
S = [Sx, Sy];
end
