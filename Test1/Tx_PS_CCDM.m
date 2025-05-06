function [Stx,txSyms,R_CCDM] = Tx_PS_CCDM(C,H,txBits)
% Last Update: 08/01/2019
% Apply Distribution Matcher   Ӧ�÷ֲ�ƥ����
% ����Tx_ProbShaping
%% Input Parameters  �������
M = numel(C);
[nPol,nBits] = size(txBits);%npol ���⼫����
nSyms = ceil(nBits/log2(M));%nsyms ģ���������  
%ceil:  w=ceil(z)����������z�е�Ԫ��ȡ����ֵwΪ��С�ڱ������С������
...���ڸ���B���ֱ��ʵ�����鲿ȡ����

%% Assign Symbol Probability According to Maxwell-Boltzman Distribution  �������˹Τ-���������ֲ�������Ÿ���
lambda = entropy2lambda(H,C);%��һ���ؼ�function:entropy2lambda
symProb = exp(-lambda*abs(C).^2);

%% Initialize CCDM  ��ʼ�� CCDM
[symProb,nBitsInfo,symFreq] = ccdm.initialize(symProb,nSyms);  %initialize���� ���õݹ� n_choose_ks_recursive_log2( n,k )ʵ��֧�ֵ����볤���������������������
R_CCDM = nBitsInfo/nBits;  %%%%%��

%% Encode with Distribution Matcher ʹ�÷ֲ�ƥ�������б���
[Stx,txSyms] = deal(NaN(nPol,nSyms));  %���ɴ������   deal��ͬʱ�����������ֵ
for n = 1:nPol
    i_TX = ccdm.encode(txBits(n,1:nBitsInfo),symFreq).' + 1;  %ccdm.encode
    Stx(n,:) = C(i_TX).';
    txSyms(n,:) = i_TX.'-1;
end

