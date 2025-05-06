%MINIEXAMPLE
% Run this basic example to understand how to initialize, encode and decode
% a message using CCDM. You can adjust the choosen output distribution and
% the output length. 运行此基本示例以了解如何使用 CCDM 初始化、编码和解码消息。 您可以调整选择的输出分布和输出长度。
%
% See also CCDM.INITIALIZE, CCDM.ENCODE, CCDM.DECODE, INSTALL

% 4. In case results are published that rely on this source code, please cite
%    our paper entitled "Constant Composition Distribution Matching" [1]. 
% [1] http://arxiv.org/abs/1503.05133

% choose aribtray target distribution and output length 选择 aribtray 目标分布和输出长度
pOpt = [0.0,0.2,0.3,0.5];
n = 1000000;
% calculate  input length m, and the optimal n-type approximation  计算输入长度 m，以及最优 n 型近似
[p_quant,num_info_bits,n_i] = ccdm.initialize(pOpt,n);

% generate uniform bits of input length m  生成输入长度为 m 的统一位
src_symbols = randi(2,1,num_info_bits)-1;
% encode with distribution matcher 使用分布匹配器进行编码
tic;
code_symbols = ccdm.encode(src_symbols,n_i);
% decode with distribution matcher 使用分布匹配器解码
src_symbols_hat = ccdm.decode(code_symbols,n_i,num_info_bits);
toc;
% check equality 检查相等
display(sum(src_symbols_hat ~= src_symbols));
% check distribution 检查分配
hist(code_symbols,0:length(p_quant)-1);
