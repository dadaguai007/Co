%ENCODE
% code_symbols = ENCODE(src_symbols,n_i) 

% src_symbols is a binary vector of length num_info_bits (see
% ccdm.initialize). 
...src_symbols 是长度为 num_info_bits 的二进制向量（参见 ccdm.initialize）
    
%n_i is the counting vector as returned by ccdm.initialize. 
...n_i 是由ccdm.initialize 返回的计数向量。
% See also CCDM.INITIALIZE, CCDM.DECODE, MINIEXAMPLE

%CCDM.ENCODE
% encode a vektor of source symbols src_symbols into a Type(n_i) vector.
...将源符号 src_symbols 的向量编码为 Type(n_i) 向量。
% a basic example is provided in miniexample
% 
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



function code_symbols = encode(src_symbols,n_i)    %源符号 src_symbols    ...将源符号 src_symbols 的向量编码为 Type(n_i) 向量 
%encoding src_symbols to n_i type sequence. Consider ccdm.initialize as well
...将 src_symbols 编码为 n_i 类型序列。同时考虑 ccdm.initialize   
  code_symbols = double(ccdm.encodeCCADM(int32(src_symbols), int32(n_i))-1);
end