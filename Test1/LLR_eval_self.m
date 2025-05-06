clc;
clear all; 
%%
load('C:\Users\Administrator\Desktop\work\emulation\PS MATLAB CODE\16QAM');
Constellation
C=Constellation;

load('C:\Users\Administrator\Desktop\work\emulation\PS MATLAB CODE\receiveGG5_1');
receiveGG5_1;
y = receiveGG5_1;  %support_Yֻ��Ҫ�ź�������,����Ҫ���ʷֲ����ź� 
%% Input Parameters
nSyms = size(y,2);   %�Ѿ���y��������ֵ��nSyms
M = numel(C);          %M��������С
m = log2(M);        %ÿ���Ŷ���bit
bitMap = dec2bin(0:M-1,m);  
 
PX=[0.0247 0.0539 0.0539 0.0247 0.0539 0.1175  0.1175 0.0539  0.0539 0.1175 0.1175 0.0539 0.0247 0.0539 0.0539 0.0247 ];%������ʳ���
% PX = [0.0406,0.0601,0.0601,0.0406,0.0601,0.0892, 0.0892,0.0601, 0.0601,0.0892,0.0892,0.0601,0.0406,0.0601,0.0601,0.0406];%H=3.9451 
%----------------------------------------------------------------------------
load('C:\Users\Administrator\Desktop\work\emulation\PS MATLAB CODE\16QAM')
Constellation
X2D=Constellation;%��ά����  
X2D=X2D.';
d1=(real(X2D)).^2; % A������i�� ��ȡ����A�ĵ�i��  
d2=(imag(X2D)).^2;
snr_dB =2:20; 
snr_lin = 10.^(snr_dB/10);
m   = log2(M);              %ÿ�����ŵĵ�ַλ
Eav = sum(PX .*(d1+d2));    %��άƽ������ �ź�Xi�ܵ���ƽ������Լ�� ������ȷ���� y �� C�����磬��׼�����ǵ�ƽ�����ʻ���С�����ǵ� MMSE��
%% Calculate LLRs
% for isnr = 1:length(snr_lin) 
isnr =4;
    %---------------------------------------
    %��ά�źŵ��������ʣ������ŵ���
    noise_var     = Eav/snr_lin(isnr);  %noise_var=sigma^2
    noise_std_dev = sqrt(noise_var);
    inv2Var       = 1 / (2*noise_var);
    invVar        = 1 / (noise_var);
    invSq2piVar   = 1 / sqrt(2*pi*noise_var); 
    invSqpiVar    = 1 / sqrt(pi*noise_var);
    %----------------------------------------
LLRs = NaN(nSyms*m,1);
for n = 1:m                      % n��1��4��16QAM�� ÿ���Ŷ��ٱ���
    PB = zeros(m,2);            % bit-level distributions  bit���ֲ�
    idx0 = bitMap(:,n)=='0';  
    idx1 = bitMap(:,n)=='1';  
for ibit = 1:m
    PB(ibit,0+1) = sum(PX(idx0)) ;        % being zero
    PB(ibit,1+1) = sum(PX(idx1)) ;        % being one
end
    %��ȡ�Ӽ� Xk_b �� Pk_b��
    idx = bitMap(:,n)=='0';        
    Xk_0 = C(idx);                  %��ʾ��nλΪ 0 ��������ż���  Xk_0 = X2D(labeling(:,ibit)==0);% Bi
    Pk_0 =( PX(idx))/PB(ibit,0+1);
    idx = bitMap(:,n)=='1';         
    Xk_1 = C(idx);                  %��ʾ��nλΪ 1 ��������ż���  Xk_1 = X2D(labeling(:,ibit)==1);
    Pk_1 =( PX(idx))/PB(ibit,1+1);
 
    % ������Ȼ�ȵķ��Ӻͷ�ĸ��
    %-----------------
    A = zeros(1,nSyms);
    for k = 1:numel(Xk_0)
         A = A + exp((-abs(y-Xk_0(k)).^2)*invVar)*Pk_0(k); %���  ������ȷ���� y �� C�����磬��׼�����ǵ�ƽ�����ʻ���С�����ǵ� MMSE��
% p_ygb0(iy) = p_ygb0(iy) + PX0(ix) * invSqpiVar * exp( -(y-Xk_0(ix))^2 * invVar);  % invVar ��ʾ 1/sigma^2��1/N0 ...��ŷ�Ͼ���
    end
     %-----------------
    B = zeros(1,nSyms);
    for k = 1:numel(Xk_1)
        B = B + exp((-abs(y-Xk_1(k)).^2)*invVar)*Pk_1(k); %���  ������ȷ���� y �� C�����磬��׼�����ǵ�ƽ�����ʻ���С�����ǵ� MMSE��
% p_ygb1(iy) = p_ygb1(iy) + PX1(ix) * invSqpiVar * exp( -(d7+d8) * invVar);  %invSq2piVar * exp( -(y-X1D_b1(ix))^2 * inv2Var)Ϊ�������ת�Ƹ���
    end
     %-----------------
    %���������Ȼ��:
    LLRs(n:m:end) = -log(B./A);  % ע��������n:nBpS:end ��Ϊÿ�μ������ÿ�����ŵĵ�1�����أ���2������
end
% end