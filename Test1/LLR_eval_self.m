clc;
clear all; 
%%
load('C:\Users\Administrator\Desktop\work\emulation\PS MATLAB CODE\16QAM');
Constellation
C=Constellation;

load('C:\Users\Administrator\Desktop\work\emulation\PS MATLAB CODE\receiveGG5_1');
receiveGG5_1;
y = receiveGG5_1;  %support_Y只需要信号与噪声,不需要概率分布的信号 
%% Input Parameters
nSyms = size(y,2);   %把矩阵y的列数赋值给nSyms
M = numel(C);          %M：星座大小
m = log2(M);        %每符号多少bit
bitMap = dec2bin(0:M-1,m);  
 
PX=[0.0247 0.0539 0.0539 0.0247 0.0539 0.1175  0.1175 0.0539  0.0539 0.1175 0.1175 0.0539 0.0247 0.0539 0.0539 0.0247 ];%引入概率成形
% PX = [0.0406,0.0601,0.0601,0.0406,0.0601,0.0892, 0.0892,0.0601, 0.0601,0.0892,0.0892,0.0601,0.0406,0.0601,0.0601,0.0406];%H=3.9451 
%----------------------------------------------------------------------------
load('C:\Users\Administrator\Desktop\work\emulation\PS MATLAB CODE\16QAM')
Constellation
X2D=Constellation;%二维星座  
X2D=X2D.';
d1=(real(X2D)).^2; % A（：，i） 提取矩阵A的第i列  
d2=(imag(X2D)).^2;
snr_dB =2:20; 
snr_lin = 10.^(snr_dB/10);
m   = log2(M);              %每个符号的地址位
Eav = sum(PX .*(d1+d2));    %二维平均能量 信号Xi受到的平均功率约束 必须正确缩放 y 和 C（例如，标准化它们的平均功率或最小化它们的 MMSE）
%% Calculate LLRs
% for isnr = 1:length(snr_lin) 
isnr =4;
    %---------------------------------------
    %二维信号的条件概率（辅助信道）
    noise_var     = Eav/snr_lin(isnr);  %noise_var=sigma^2
    noise_std_dev = sqrt(noise_var);
    inv2Var       = 1 / (2*noise_var);
    invVar        = 1 / (noise_var);
    invSq2piVar   = 1 / sqrt(2*pi*noise_var); 
    invSqpiVar    = 1 / sqrt(pi*noise_var);
    %----------------------------------------
LLRs = NaN(nSyms*m,1);
for n = 1:m                      % n从1到4（16QAM） 每符号多少比特
    PB = zeros(m,2);            % bit-level distributions  bit级分布
    idx0 = bitMap(:,n)=='0';  
    idx1 = bitMap(:,n)=='1';  
for ibit = 1:m
    PB(ibit,0+1) = sum(PX(idx0)) ;        % being zero
    PB(ibit,1+1) = sum(PX(idx1)) ;        % being one
end
    %获取子集 Xk_b 和 Pk_b：
    idx = bitMap(:,n)=='0';        
    Xk_0 = C(idx);                  %表示第n位为 0 的输入符号集合  Xk_0 = X2D(labeling(:,ibit)==0);% Bi
    Pk_0 =( PX(idx))/PB(ibit,0+1);
    idx = bitMap(:,n)=='1';         
    Xk_1 = C(idx);                  %表示第n位为 1 的输入符号集合  Xk_1 = X2D(labeling(:,ibit)==1);
    Pk_1 =( PX(idx))/PB(ibit,1+1);
 
    % 计算似然比的分子和分母：
    %-----------------
    A = zeros(1,nSyms);
    for k = 1:numel(Xk_0)
         A = A + exp((-abs(y-Xk_0(k)).^2)*invVar)*Pk_0(k); %求和  必须正确缩放 y 和 C（例如，标准化它们的平均功率或最小化它们的 MMSE）
% p_ygb0(iy) = p_ygb0(iy) + PX0(ix) * invSqpiVar * exp( -(y-Xk_0(ix))^2 * invVar);  % invVar 表示 1/sigma^2或1/N0 ...用欧氏距离
    end
     %-----------------
    B = zeros(1,nSyms);
    for k = 1:numel(Xk_1)
        B = B + exp((-abs(y-Xk_1(k)).^2)*invVar)*Pk_1(k); %求和  必须正确缩放 y 和 C（例如，标准化它们的平均功率或最小化它们的 MMSE）
% p_ygb1(iy) = p_ygb1(iy) + PX1(ix) * invSqpiVar * exp( -(d7+d8) * invVar);  %invSq2piVar * exp( -(y-X1D_b1(ix))^2 * inv2Var)为输入输出转移概率
    end
     %-----------------
    %计算对数似然比:
    LLRs(n:m:end) = -log(B./A);  % 注意这里是n:nBpS:end 因为每次计算的是每个符号的第1个比特，第2个比特
end
% end