% PMD 建模
clc;clear;close all;
N=12;%段数
tau_mean = 3; % average tau 3ps
sigma = 0.01; % variance
% 每个光纤段生成一个随机的初始相位，模拟光脉冲在光纤中的随机偏振态变化

Tau = sqrt(3*pi/(8*N))*(1+sigma*randn(N,1)).*tau_mean; %ti的取值,分为N段
fs=10e3;
% 时间轴
T_start = 0;
T_end = 5e-2-1/fs;
T = T_start:1/fs:T_end;
% RSOP 随时间变化，初始参数
w_arfa = 111e3; %rsop w
w_kappa = 61e3;
w_phi = 115e3;
kapa0 = pi+randn(N+1,1)*2*pi;
arfa0 = pi+randn(N+1,1)*2*pi;
phi0 = pi+randn(N+1,1)*2*pi;
% 频率轴
w_start = 1905;%1594nm            取中心波长的值
w_stop = 1932;%1530nm             取中心波长加上信号的频宽的值
w_spacing = 0.03; % angular freq spacing, unit: Trad/s       是不是应该取fs/length(signal)的长度
w = w_start:w_spacing:w_stop; % unit: T rad/s


% 随时间变化的RSOP矩阵，每个时间下，应为N+1 长度
% 横轴为时间T，纵轴为段数N
for i=1:length(T)

    arfa(:,i) = w_arfa.*T(i)+arfa0;
    phi(:,i) = w_phi.*T(i)+phi0;
    kapa(:,i) = w_kappa.*T(i)+kapa0;

    oo11(:,i)  = cos(kapa(:,i)).*exp(1j*(arfa(:,i)));
    oo12(:,i)  = -sin(kapa(:,i)).*exp(-1j*(phi(:,i)));
    oo21(:,i)  = -conj(oo12(:,i));
    oo22(:,i)  = conj(oo11(:,i));

end

% PMD矩阵，随频率变化
% 横轴为频率w，纵轴为光纤段数
for j=1:length(w)
    ee(:,j) = exp(1j*w(j).*Tau/2); % 群延时矩阵
    ff (:,j)= conj(ee(:,j));
end

% PMD动态建模
% U_i 为横轴为时间T，纵轴为频率w

omega=cell(1,length(w));
for j=1:length(w)
    for i=1:length(T)
        U=single(:,i);
        for K=1:N
            H=[oo11(K,i),oo12(K,i);oo21(K,i),oo22(K,i)];
            B=[ee(K,j),0;0,ff(K,j)];
           
            U1=H*U;
            U_f=fft(U1);
            U1=ifft(B*U_f);
            U=U1;
        end
        H_N1=[oo11(N+1,i),oo12(N+1,i);oo21(N+1,i),oo22(N+1,i)];
        single_out(:,i)=U*H_N1;
    end
    omega{j}=single_out;
end

