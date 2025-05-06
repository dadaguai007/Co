function [Stx,txSyms,R_CCDM] = Tx_PS_CCDM(C,H,txBits)
% Last Update: 08/01/2019
% Apply Distribution Matcher   应用分布匹配器
% 用于Tx_ProbShaping
%% Input Parameters  输入参数
M = numel(C);
[nPol,nBits] = size(txBits);%npol 激光极化数
nSyms = ceil(nBits/log2(M));%nsyms 模拟符号总数  
%ceil:  w=ceil(z)函数将输入z中的元素取整，值w为不小于本身的最小整数。
...对于复数B，分别对实部和虚部取整。

%% Assign Symbol Probability According to Maxwell-Boltzman Distribution  根据麦克斯韦-玻尔兹曼分布分配符号概率
lambda = entropy2lambda(H,C);%第一个关键function:entropy2lambda
symProb = exp(-lambda*abs(C).^2);

%% Initialize CCDM  初始化 CCDM
[symProb,nBitsInfo,symFreq] = ccdm.initialize(symProb,nSyms);  %initialize量化 利用递归 n_choose_ks_recursive_log2( n,k )实际支持的输入长度来计算期望的输出长度
R_CCDM = nBitsInfo/nBits;  %%%%%？

%% Encode with Distribution Matcher 使用分布匹配器进行编码
[Stx,txSyms] = deal(NaN(nPol,nSyms));  %生成传输符号   deal：同时给多个变量赋值
for n = 1:nPol
    i_TX = ccdm.encode(txBits(n,1:nBitsInfo),symFreq).' + 1;  %ccdm.encode
    Stx(n,:) = C(i_TX).';
    txSyms(n,:) = i_TX.'-1;
end

