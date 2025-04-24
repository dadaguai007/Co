function [z] = slice(y, M)
% 符号判决（映射到最近星座点）
% 输入:
%   y - 接收符号（可能为复数）
%   M - 调制阶数（如16-QAM对应M=16）
% 输出:
%   z - 判决后的符号

if isreal(y) % PAM调制
    z_index = round((real(y) + M - 1)/2); % 量化到整数索引
    z_index = max(0, min(M-1, z_index));  % 边界截断
    z = z_index*2 - (M-1);                % 恢复为±1, ±3...
else % QAM调制
    M_bar = sqrt(M); % 每维度星座点数（如16-QAM为4）
    % 实部量化
    idx_re = round((real(y) + M_bar - 1)/2);
    idx_re = max(0, min(M_bar-1, idx_re));
    % 虚部量化
    idx_im = round((imag(y) + M_bar - 1)/2);
    idx_im = max(0, min(M_bar-1, idx_im));
    % 组合为复符号
    z = (idx_re*2 - (M_bar-1)) + 1j*(idx_im*2 - (M_bar-1));
end
end
