function y = firFilter(h, x)
% Perform FIR filtering and compensate for filter delay.

[M, N] = size(x);
if M == 1
    x = reshape(x, N, []);
end

y = x;
%按列进行排列
nModes = size(x, 2);

for n = 1:nModes
    % 直接使用conv中same满足延时补偿
    y(:, n) = conv(x(:, n), h, 'same');
end
% vector
if size(y, 2) == 1
    y = y';
end
end
% 这段代码是一个MATLAB函数，名为firFilter，用于执行FIR（有限冲激响应）滤波并补偿滤波器延迟。