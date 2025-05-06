function [syms] = bit2sym(bits,nBpS)
%bit2sym    Transform bits into symbol indices 将比特转换为符号索引
%   This function takes a stream of transmitted/received bits and
%   calculates the corresponding symbol indices for a given number of bits
%   per symbol.
... 该函数采用传输/接收的比特流，并为每个符号的给定比特数计算相应的符号索引。
    
%   INPUTS:
%   bits    :=      array of bits [Nsig x Nbits]  位数组
%   nBps    :=      number of bits per symbol of a given constellation  给定星座的每个符号的位数
%
%   OUTPUTS:
%   syms    :=      array of symbol indices [Nsig x Nsyms], taking values
%                   in [0 2^nBps-1] 
%                   ... 符号索引数组    取 [0 2^nBps-1] 中的值
%   Fernando Guiomar
%   Last Update: 02/03/2017

%% Validate Input Arguments  验证输入参数
validateattributes(bits,{'logical'},{'binary'},'','bits',1);
validateattributes(nBpS,{'numeric'},{'scalar','positive','integer'},'','nBpS',2);

%% Input Parameters
[nSig,nBits] = size(bits);
nSyms = nBits / nBpS;

%% Bit-to-Symbol Assignment 比特到符号分配
syms = zeros(nSig,nSyms);
for k = 1:nSig
    for n = 1:nBpS
        syms(k,:) = syms(k,:) + bits(k,n:nBpS:end)*2^(nBpS-n);
    end
end
