function [MI] = MI_eval(Srx,Stx,C,N0,symProb,useGPU)
%MI_eval    Estimate the mutual information (MI) in AWGN channels
%  估计 AWGN 通道中的互信息 (MI)
%   This function evaluates the MI between the received signal and the
%   transmitted signal, considering AWGN channel statistics.
%   This function also supports non-uniform modulation (probabilistic
%   shaping) in which the constellation points have different probabilities
%   of occurrence. 
%  该函数评估接收信号和发射信号之间的 MI，考虑 AWGN 信道统计。 此功能还支持星座点具有不同出现概率的非均匀调制（概率整形）。
%   INPUTS:
%   Srx     := received signal [1 x nSyms]
%   Stx     := transmitted signal [1 x nSyms]
%   C       := reference constellation points [M x 1]
%   N0      := noise variance of Srx [1 x 1]
%   symProb := symbol probabilities of each constellation point [M x 1]
%   useGPU  := flag to decide whether to use or not GPU-based processing
%
%   Note that Srx and Stx must be properly synchronized and scaled before
%   calling the MI_eval function. Srx and C must also be properly scaled
%   (e.g. normalizing their average power or minimizing their MMSE).
%
%   OUTPUTS:
%   MI      := mutual information [1 x 1]
%
%
%   Author: Fernando Guiomar
%   Last Update: 04/06/2019

%% Input Parser
if nargin < 5
    symProb = repmat(1/numel(C),numel(C),1);
end

%% Cast to Single (for increased speed)
Srx = single(Srx);
Stx = single(Stx);
C = single(C);
N0 = single(N0);
symProb = single(symProb);

%% Transfer Variable to GPU
if nargin < 6
    useGPU = false;
end
if useGPU
    Srx = gpuArray(Srx);
    Stx = gpuArray(Stx);
    C = gpuArray(C);
    N0 = gpuArray(N0);
    symProb = gpuArray(symProb);
end

%% MI Estimation Considering an AWGN Channel
qYonX = exp((-abs(Srx-Stx).^2)/N0);
qY = sum(repmat(symProb,1,numel(Srx)).*exp((-abs(Srx-C).^2)/N0));
MI = mean(log2(max(qYonX,realmin)./max(qY,realmin)));

%% Gather Variables from GPU
if useGPU
    MI = gather(MI);
end

