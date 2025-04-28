function [Frequency,Frequency_M,Y,DFT] = Get_tone_CLK(X,fs)
% 对单音信号的最大似然函数估计

%Coarch search
N = length(X);
M=4*N;
X_fft = fft(X,M);
[~,k] = max(abs(X_fft(1:M/2-1)));
wc = k*fs/M;
n=1:N;
%Fine search
% 定义左边和右边的边界
index=1e2;
bounds=2e4;
left_boundary = wc - bounds;
right_boundary = wc + bounds;

% 生成频率数组
fre_array = left_boundary:index:right_boundary;

% DFT_max=abs(sum(X.*exp(-1j*2*pi*wc*n/fs)')./N);

for j =1:1:length(fre_array)
    %似然函数
    DFT(j)=abs(sum(X.*exp(-1j*2*pi*fre_array(j)*n/fs)')./N);
end

[~,L]=min(DFT);
[~,M]=max(DFT);
%可根据情况选择拟合数据的长度
H=DFT(1:L-10);
Fre=fre_array(1:L-10);
%拟合，更精确估计
polyForWave = polyfit(Fre,H, 2);
%导数
cofficient=polyder(polyForWave);
Frequency = roots(cofficient);
Frequency_M=fre_array(M);
Y=abs(sum(X.*exp(-1j*2*pi*Frequency*n/fs)')./N);

end