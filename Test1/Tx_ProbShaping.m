function [Stx,txSyms,QAM] = Tx_ProbShaping(txBits,QAM,SIG,R_FEC)

% Last Update: 15/10/2019


%% Input Parser 输入解析器
% Set Default Shaping Method:    设置默认整形方法
if ~isfield(QAM,'PCS') || ~isfield(QAM.PCS,'method')%PCS是否为QAM的域 method是否为QAM.PCS的域 如果是 则CCDM
    QAM.PCS.method = 'CCDM';
end
% Set Default FEC Rate:          设置默认 FEC 速率
if nargin < 4 || isempty(R_FEC)  %nargin是用来判断输入变量个数的函数 若变量个数<4且 R_FEC域为空 则返回设置速率为1
    R_FEC = 1;
end

%% Impact of FEC on the PAS Scheme      FEC 对 PAS方案的影响
nBpS = SIG.nBpS;   %Nbps = 2; M = 2^Nbps;     %调制阶数 = 2/4/6; QPSK/16-QAM/64-QAM
C = QAM.IQmap;%#
M_PS = numel(C);% n = numel(C); % 返回数组C中元素个数
flag = false;%赋给flag“否”这么一个状态
if mod(log2(M_PS),1)
    M_PS = 2^nextpow2(M_PS);
    if mod(sqrt(M_PS),1)
        M_PS = 2^nextpow2(M_PS+1);
    end
    flag = true;
end
H_quad = (1-R_FEC)*log2(M_PS);
R_FEC_min = (log2(M_PS) - 2) / log2(M_PS);
if H_quad > 2
    error(['PAS scheme requires to allocate FEC bits to quadrant ',...
        'positions, whose maximum entropy is 2 bits/symbol. ',...
        'With the requested FEC rate of ',num2str(R_FEC,'%1.2f'),...
        ', the required entropy for FEC bits is ',...
        num2str(H_quad,'%1.2f'),', which exceeds the 2 bits/sym limit!',...
        ' Consider increasing the FEC rate. The minimum FEC rate '...
        'allowed for this system is (log2(M)-2)/log2(M)',...
        num2str(R_FEC_min,'%1.4f'),'.']);
end
H_PAS = nBpS * R_FEC + H_quad;
H_DM = H_PAS - 2;
H_PS = nBpS * R_FEC;

%% Apply Distribution Matcher   应用分布匹配器
[Stx,txSyms,QAM.R_CCDM] = Tx_PS_CCDM(C,H_PAS,txBits);

%%
if flag
    QAM.M = M_PS;
    QAM.IQmap = qammod(0:M_PS-1,M_PS,'gray').';
    txSyms = signal2symbol(Stx,QAM.IQmap);
    nBpS = log2(M_PS);
    bMap = false(M_PS,nBpS);
    for n = 0:M_PS-1
        [~,e] = log2(n);
        bMap(n+1,:) = rem(floor(n * pow2(1-max(e,nBpS):0)),2);
    end
    QAM.sym2bitMap = bMap;
    QAM.nBpS = log2(M_PS);
end

%% Set QAM Parameters
QAM.maxConstPower = max(abs(Stx).^2);
QAM.meanConstPower = mean(abs(Stx).^2);
QAM.maxConstPower = max(abs(Stx).^2);
edges = -0.5:1:QAM.M-0.5;
QAM.symProb = histcounts(txSyms(:),edges,'Normalization','prob').'; %histcounts 根据设置的区间，根据概念来计算pdf，受数据影响
tmp = log2(QAM.symProb);
tmp(isinf(tmp)) = 0;
QAM.entropy = -sum(QAM.symProb.*tmp);

% % Debug
% RGB = fancyColors;
% plotProbShaping_PDF_const('const',QAM.IQmap,...
%     'symbols',txSyms(:),'color',RGB.itred);
% plotConstMap(QAM.IQmap,QAM.sym2bitMap,QAM.symProb);

% Debug
RGB = fancyColors;
plotProbShaping_PDF_const('const',QAM.IQmap,...    %PCS概率密度函数图
    'symbols',txSyms(:),'color',RGB.itred);  %A(:)：将矩阵A中的每列合并成一个长的列向量
plotConstMap(QAM.IQmap,QAM.sym2bitMap,QAM.symProb);%常量映射

