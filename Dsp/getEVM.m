function EVM = getEVM(sig, constellation)
sig = sig ./ sqrt(bandpower(sig));
constellation = constellation / sqrt(bandpower(constellation));
Pe = bandpower(sig - decision(sig, constellation));
Pref = bandpower(sig);
EVM = sqrt(Pe./Pref);
end

