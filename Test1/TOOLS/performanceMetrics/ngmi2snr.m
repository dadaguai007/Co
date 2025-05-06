function [SNR_dB] = ngmi2snr(NGMI,M,E)
%ngmi2snr   Evaluate the theoretical SNR that corresponds to given
%           measured NGMI, constellation template for probabilistic
%           shaping (M) and signal entropy (E)   
%本程序用来评估与给定测量 NGMI、概率整形 (M) 和信号熵 (E) 的星座模板相对应的理论 SNR
%
%   INPUTS:
%   NGMI    :=  normalized generalized mutual information [1 x 1]   归一化广义互信息
%   M       :=  QAM constellation sizes [1 x 1]  QAM 星座大小
%   E       :=  entropy after constellation shaping [1 x 1]  星座整形后的熵
%
%   OUTPUTS:
%   SNR_dB  :=  signal-to-noise ratio (dB) [1 x 1]   信噪比
%
%
%   Examples:
%       [SNR_dB] = ngmi2snr(0.9,256,6);
%% Input Parser
if nargin < 3
    % If no entropy is indicated, then assume uniform modulation:
    E = log2(M);
end

%% Load SNR vs NGMI Table
fileName = [num2str(M) 'QAM_snr2gmi'];
SNR_vs_GMI = load(fileName);
E_LUT = SNR_vs_GMI.entropy;
SNR_LUT = SNR_vs_GMI.SNR_dB;
NGMI_LUT = SNR_vs_GMI.NGMI;

%% Find Closest Entropy Values in the LUT
[minErr,idx] = min(abs(E_LUT-E));
if minErr > 0
    idx = idx-1:idx+1;
end
E_LUT = E_LUT(idx);

%% Interpolate Over NGMI for Each Trial Entropy
thisSNR_dB = NaN(1,numel(idx));
for n = 1:numel(idx)
    thisNGMI = NGMI_LUT(idx(n),:);
    idxNaN = isnan(thisNGMI);
    thisNGMI = thisNGMI(~idxNaN);
    thisSNR = SNR_LUT(~idxNaN);
    [thisNGMI,idxUnique] = unique(thisNGMI);
    thisSNR = thisSNR(idxUnique);
    thisSNR_dB(n) = interp1(thisNGMI,thisSNR,NGMI,'linear'); 
end

%% Interpolate Over Entropy  对熵进行插值
if numel(idx) > 1
    SNR_dB = interp1(E_LUT,thisSNR_dB,E,'linear');
else
    SNR_dB = thisSNR_dB;
end
