function [I, Q] = gsopTest(I, Q)
    I = I / sqrt(bandpower(I));
    Q = Q / sqrt(bandpower(Q));
    Q = Q - mean(I.*Q)*I/bandpower(I);
    Q = Q / sqrt(bandpower(Q));
end