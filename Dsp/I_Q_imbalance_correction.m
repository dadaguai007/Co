% 补偿IQ不平衡


I = real(Ein);
Q = imag(Ein);
%I和Q的相关性，计算I和Q的功率。然后，对Q进行校正，以减少I和Q之间的相关性
rho = mean(I.*Q);
P_I =pwr.meanpwr(I);
Q = Q - rho*I/P_I;
P_Q = pwr.meanpwr(Q);

Eout = pwr.normpwr(I/sqrt(P_I) + 1j*Q/sqrt(P_Q));

%save imbalance
idx = find(isnan(obj.results.rho), 1);
obj.results.rho(idx) = rho;