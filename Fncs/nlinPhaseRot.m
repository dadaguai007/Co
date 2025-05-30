function phi = nlinPhaseRot(Ex, Ey, Pch, gamma)
% Calculate nonlinear phase-shift per step for the Manakov SSFM.
%
% Ex : 
%     Input optical signal field of x-polarization.
% Ey : 
%     Input optical signal field of y-polarization.
% Pch : 
%     Input optical power.
% gama : 
%     fiber nonlinearity coefficient.
%     nonlinear phase-shift of each sample of the signal.
phi = real((8 / 9) * gamma * (Pch + Ex .* conj(Ex) + Ey .* conj(Ey)) / 2);
end

