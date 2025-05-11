function b = index_to_binary(iTX)
% 将一个十进制索引向量转换为二进制向量。
% 功能描述
% 输入：iTX 是一个列向量，包含非负整数（十进制索引）。
% 输出：b 是一个二进制矩阵，每一行对应输入向量 iTX 中的一个十进制数的二进制表示。
% 目标：将每个十进制数转换为固定长度的二进制向量，长度为 N_bit，其中 N_bit 是能表示输入向量中最大值所需的二进制位数。

%% Convert
N_bit = ceil(log2(max(iTX)));
b = NaN(size(iTX,1),N_bit);
for i = 1:N_bit
    tmp = fix(iTX/2);
    b(:,i) = iTX - 2*tmp;
    iTX = tmp;
end
end