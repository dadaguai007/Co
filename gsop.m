function [I, Q] = gsop(I, Q)
    I = I / sqrt(bandpower(I));
    Q = Q / sqrt(bandpower(Q));
    Q = Q - mean(I.*Q)*I/bandpower(I);
    Q = Q / sqrt(bandpower(Q));
end