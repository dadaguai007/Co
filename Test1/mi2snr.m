function [SNR_dB] = mi2snr(MI,M,E)
%ngmi2snr   Evaluate the theoretical SNR that corresponds to given
%           measured MI, constellation template for probabilistic
%           shaping (M) and signal entropy (E)
%      评估与给定测量 MI、用于概率整形 (M) 和信号熵 (E) 的星座模板相对应的理论 SNR
%   INPUTS:
%   MI      :=  mutual information [1 x 1]
%   M       :=  QAM constellation sizes [1 x 1]
%   E       :=  entropy after constellation shaping [1 x 1]
%
%   OUTPUTS:
%   SNR_dB  :=  signal-to-noise ratio (dB) [1 x 1]
%
%
%   Examples:
%       [SNR_dB] = mi2snr(5.2,256,6);
%
%
%   Author: Fernando Guiomar
%   Last Update: 04/07/2019

%% Input Parser
if nargin < 3
    % If no entropy is indicated, then assume uniform modulation:
    ...如果没有指示熵，则假设均匀调制：
    E = log2(M);
end

%% Load SNR vs MI Table
fileName = [num2str(M) 'QAM_snr2gmi'];
SNR_vs_GMI = load(fileName);
E_LUT = SNR_vs_GMI.entropy;
SNR_LUT = SNR_vs_GMI.SNR_dB;
MI_LUT = SNR_vs_GMI.MI;

%% Find Closest Entropy Values in the LUT
[minErr,idx] = min(abs(E_LUT-E));
if minErr > 0
    idx = idx-1:idx+1;
end
E_LUT = E_LUT(idx);

%% Interpolate Over MI for Each Trial Entropy 为每个试验熵插值 MI
thisSNR_dB = NaN(1,numel(idx));
for n = 1:numel(idx)
    thisMI = MI_LUT(idx(n),:);
    idxNaN = isnan(thisMI);
    thisMI = thisMI(~idxNaN);
    thisSNR = SNR_LUT(~idxNaN);
    [thisMI,idxUnique] = unique(thisMI);
    thisSNR = thisSNR(idxUnique);
    thisSNR_dB(n) = interp1(thisMI,thisSNR,MI,'linear'); 
end

%% Interpolate Over Entropy  对熵进行插值
if numel(idx) > 1
    SNR_dB = interp1(E_LUT,thisSNR_dB,E,'linear');
else
    SNR_dB = thisSNR_dB;
end
