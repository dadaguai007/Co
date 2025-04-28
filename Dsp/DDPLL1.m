function [wout, esphase] = DDPLL1(win, tau1, tau2, k, Fs, phase, constellation)
	[L, N] = size(win);
	win = win ./ sqrt(bandpower(win)); 
	filtercoeff = [1, 1/(2*tau1*Fs)*(1+[-1, 1]/tan(1/(2*tau2*Fs)))];
	esphase = phase * ones(L+1, N);
	wout = zeros(L, N);
	for i = 1:N
		phaseout = 0;
		filterout = 0;
		for idx = 1:L
			prephaseout = phaseout;
			wout(idx, i) = win(idx, i) * exp(-1j*esphase(idx, i));
			[~, index] = min(abs(constellation-wout(idx, i)));
			phaseout = angle(wout(idx, i)*conj(constellation(index)));
			filterout = filtercoeff * [filterout; prephaseout; phaseout];
			esphase(idx+1, i) = esphase(idx, i) + k*filterout;
		end
	end
	esphase = esphase(1:end-1, :);
end