function [wout, error, w] = my_mma_lms(win, ML, step, trainlen, constellation, sps, iter)
ntaps = ML*sps;
w = zeros(ntaps, 1);
w(round(ntaps/2)) = 1;
L = floor((length(win)-ntaps)/sps+1);
win = win / sqrt(bandpower(win));
r = uniquetol(abs(constellation/sqrt(bandpower(constellation))), 1e-8);
r2 = r.^2;
wout = zeros(L, 1);
error = zeros(iter*trainlen, 1);
for i = 1:iter
    block = zeros(length(w), 1);
    for j = 1:trainlen
        for idx = 1:sps
            block = [win(idx+(j-1)*sps); block(1:end-1)];
        end
        y = w.' * block;
        amp2 = abs(y)^2;
        [~, index] = min(abs(r2 - amp2));
        error((i-1)*trainlen+j) = r2(index) - amp2;
        w = w + step * error((i-1)*trainlen+j) * y * conj(block);
    end
end

block = zeros(length(w), 1);
for i = 1:L
    for idx = 1:sps
        block = [win(idx+(i-1)*sps); block(1:end-1)];
    end
    wout(i) = w.' * block;
end
wout = wout(ML:end);
end