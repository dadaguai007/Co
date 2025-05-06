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


function [ out ] = n_choose_ks_recursive_log2( n,k )   %递归
%Calculates log_2(n!/prod_i k(i)!)
%   n is skalar,  % n是标量
%   k is vector with sum_i k(i) = n   % k 是 sum_i k(i) = n 的向量

out = 0;  %初始化
k = sort(k);  %当k是向量时，sort(k)对k的元素进行升序排序

for i = 1:length(k)-1
    out = out + ccdm.n_choose_k_iter_log(n,k(i));
    n = n - k(i);
end

end
