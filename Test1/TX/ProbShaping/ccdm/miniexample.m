%MINIEXAMPLE
% Run this basic example to understand how to initialize, encode and decode
% a message using CCDM. You can adjust the choosen output distribution and
% the output length. ���д˻���ʾ�����˽����ʹ�� CCDM ��ʼ��������ͽ�����Ϣ�� �����Ե���ѡ�������ֲ���������ȡ�
%
% See also CCDM.INITIALIZE, CCDM.ENCODE, CCDM.DECODE, INSTALL

% 4. In case results are published that rely on this source code, please cite
%    our paper entitled "Constant Composition Distribution Matching" [1]. 
% [1] http://arxiv.org/abs/1503.05133

% choose aribtray target distribution and output length ѡ�� aribtray Ŀ��ֲ����������
pOpt = [0.0,0.2,0.3,0.5];
n = 1000000;
% calculate  input length m, and the optimal n-type approximation  �������볤�� m���Լ����� n �ͽ���
[p_quant,num_info_bits,n_i] = ccdm.initialize(pOpt,n);

% generate uniform bits of input length m  �������볤��Ϊ m ��ͳһλ
src_symbols = randi(2,1,num_info_bits)-1;
% encode with distribution matcher ʹ�÷ֲ�ƥ�������б���
tic;
code_symbols = ccdm.encode(src_symbols,n_i);
% decode with distribution matcher ʹ�÷ֲ�ƥ��������
src_symbols_hat = ccdm.decode(code_symbols,n_i,num_info_bits);
toc;
% check equality ������
display(sum(src_symbols_hat ~= src_symbols));
% check distribution ������
hist(code_symbols,0:length(p_quant)-1);
