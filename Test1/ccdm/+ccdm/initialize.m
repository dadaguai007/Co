% CCDM.INITIALIZE
% CCDM.INITIALIZE(p,n) 
% select a desired distribution vektor p, desired output length n.
... 选择一个期望的分布向量 p，期望的输出长度 n。

% [p_quant,num_info_bits,n_i] = CCDM.INITIALIZE(p,n)
... p_quant、num_info_bits、n_i三个为输出参数，p、n为输入参数

% p_quant is the n-type distribution approximating p.
...p_quant 是近似于 p 的 n 型分布。

% num_info_bits is the input length for the encoder.
...num_info_bits 是编码器的输入长度。
    
% n_i is the counting vector whose jth entry is the number of times the
% letter j will occure in the output sequences of ccdm.encode. Note that
% p_quant = n_i./n.
... n_i 是计数向量，其第 j 个条目是字母 j 在 ccdm.encode 的输出序列中出现的次数。 注意 p_quant = n_i./n。

% See also CCDM.ENCODE, CCDM.DECODE, MINIEXAMPLE

% Copyright (c) 2015, Patrick Schulte, Georg B枚cherer
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without modification,
% are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
%    this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation 
%    and/or other materials provided with the distribution.
%
% 3. Neither the name of the copyright holder nor the names of its contributors may
%    be used to endorse or promote products derived from this software without 
%    specific prior written permission.
%
% 4. In case results are published that rely on this source code, please cite
%    our paper entitled "Constant Composition Distribution Matching" [1]. 
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
% IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
% INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
% NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
% OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
% WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.
%
% [1] http://arxiv.org/abs/1503.05133


function [p_quant,num_info_bits,n_i] = initialize(p,n) % num_info_bits 是编码器的输入长度
% quantize to n-type distribution  量化为 n 型分布
[n_i,p_quant] = ccdm.idquant(p,n);% p_quant 是近似于 p 的 n 型分布     n_i 是计数向量

% calculate actual supported input length  计算实际支持的输入长度
num_info_bits = floor(ccdm.n_choose_ks_recursive_log2(n,n_i));   % y = floor(x) 函数将x中元素取整，值y为不大于本身的最小整数。对于复数，分别对实部和虚部取整
end