function [wout, error, w] = mimo2_2_mma(win, ML, step, trainlen, constellation, sps, iter)
ntaps = ML*sps;

w11 = zeros(ntaps, 1);
w12 = zeros(ntaps, 1);
w21 = zeros(ntaps, 1);
w22 = zeros(ntaps, 1);

w11(round(ntaps/2)) = 1;
w22(round(ntaps/2)) = 1;

L = floor((length(win(:, 1))-ntaps)/sps+1);

wout = zeros(L, 2);
error = zeros(iter*trainlen, 2);

r = uniquetol(abs(real(constellation)/sqrt(bandpower(constellation))), 1e-8);
r2 = r.^2;

for i = 1:iter
    block = zeros(length(w11), 2);

    for j = 1:trainlen
        for idx = 1:sps
            block = [win(idx+(j-1)*sps, :); block(1:end-1, :)];
        end

        y(1) = w11.' * block(:, 1) + w12.' * block(:, 2);
        y(2) = w21.' * block(:, 1) + w22.' * block(:, 2);

        amp2 = abs(y).^2;

        [~, index1(j)] = min(abs(r2 - amp2(1)));
        [~, index2(j)] = min(abs(r2 - amp2(2)));

        error((i-1)*trainlen+j, 1) = r2(index1(j)) - amp2(1);
        error((i-1)*trainlen+j, 2) = r2(index2(j)) - amp2(2);

        w11 = w11 + step(1) * error((i-1)*trainlen+j, 1) * y(1) * conj(block(:, 1));
        w12 = w12 + step(2) * error((i-1)*trainlen+j, 1) * y(1) * conj(block(:, 2));
        w21 = w21 + step(1) * error((i-1)*trainlen+j, 2) * y(2) * conj(block(:, 1));
        w22 = w22 + step(2) * error((i-1)*trainlen+j, 2) * y(2) * conj(block(:, 2));
    end

end

block = zeros(length(w11), 2);

for i = 1:L
    for idx = 1:sps
        block = [win(idx+(i-1)*sps, :); block(1:end-1, :)];
    end
    wout(i, 1) = w11.' * block(:, 1) + w12 .' * block(:, 2);
    wout(i, 2) = w21.' * block(:, 1) + w22 .' * block(:, 2);
end
wout = wout(ML:end, :);

w(:, 1) = w11;
w(:, 2) = w12;
w(:, 3) = w21;
w(:, 4) = w22;

end