clear;
close all;
clc;

% constant
hp = 6.6260657e-34;                        % [J*s]
c0 = 299792458;                            % [m/s]
lambda = 1550*1e-9;                        % [m]
CenterFrequency = c0/lambda;               % [Hz]

M = 4;
lengthSequence = 10000;
sps =2 ;
roll_off = 0.1;

% DD-LMS
mu_DDLMS = 3e-3;
ntaps_DDLMS = 31;

% ==========================Tx=========================%
data = randi([0,M-1], lengthSequence,2);
sig_ori = qammod(data, M, 'UnitAveragePower', true);
%sig_ori = pskmod(data,M,pi/M);
sig_x = sig_ori(:,1);
sig_y = sig_ori(:,2);
% -----------------------------------------------------%
figure;
plot(real(sig_x), imag(sig_x),'.','MarkerSize',15);
hold on;
plot(real(sig_y), imag(sig_y),'.','MarkerSize',15);
title('QPSK');
legend('X-Pol', 'Y-Pol');
% ---------------------upsamping-----------------------%
sig_resample = upsample(sig_ori, sps);
% ---------------------pulse shaping-------------------%
filterSymbolLength = 32; % 2*10*16+1
filterCoeffs2 = rcosdesign(roll_off,filterSymbolLength, sps,'sqrt').';

sig_ps_x = conv(sig_resample(:,1), filterCoeffs2, 'same');
sig_ps_x = (sig_ps_x-mean(sig_ps_x))./sqrt(mean(abs(sig_ps_x).^2));
sig_ps_y = conv(sig_resample(:,2), filterCoeffs2, 'same');
sig_ps_y = (sig_ps_y-mean(sig_ps_y))./sqrt(mean(abs(sig_ps_y).^2));
sig_ps = [sig_ps_x, sig_ps_y];

% -----------------------------------------------------%
figure;
plot(real(sig_ps(:,1)), imag(sig_ps(:,1)),'.');
hold on;
plot(real(sig_ps(:,2)), imag(sig_ps(:,2)),'.');
title('rrc');
legend('X-Pol', 'Y-Pol');
fc = 193.5*1e12;
%% RSOP+pmd
fc = 193.5*1e12; %ceter frequency THz
tau_mean = 3*1e-12; % average tau 3ps
sigma = 0.01; % variance
B = 10e9;
fs = 2*B;
rr = 100;
t_spacing =1/fs*rr;
T = floor(length(sig_ps_x)/rr);
t = 0:t_spacing:T*t_spacing;
N =1;
w_alpha = 0; %krad/s---rad/s
w_phi = 0;
w_kappa = 0;
alpha= zeros(N,1);
phi = zeros(N,1);
kappa = zeros(N,1);
U = cell(1,T);
kappa0 = -pi+randn(N,1)*2*pi;
alpha0 = -pi+randn(N,1)*2*pi;
phi0 = -pi+randn(N,1)*2*pi;
Tau = sqrt(3*pi/(8*N))*(1+sigma*randn(N,1)).*tau_mean; %ti的取值

% 构造PMD矩阵B（频域，每个频率点对应一个2x2矩阵）
B = zeros(2, 2, N);

for a = 1:T
    alpha = w_alpha.*t(a)+alpha0;
    phi = w_phi.*t(a)+phi0;
    kappa = w_kappa.*t(a)+kappa0;
    oo11 = cos(kappa).*exp(1j*(alpha));
    oo12 = -sin(kappa).*exp(-1j*(phi));
    oo21 = -conj(oo12);
    oo22 = conj(oo11);
    x_rsop = sig_ps_x((a-1)*rr+1:a*rr);
    y_rsop = sig_ps_y((a-1)*rr+1:a*rr);
    for b = 1:N
        if(b == 1)
            X_rsop = oo11(b).*x_rsop+oo12(b).*y_rsop;
            Y_rsop = oo21(b).*x_rsop+oo22(b).*y_rsop;
        else
            X_rsop =  oo11(b).*X_pmd+oo12(b).*Y_pmd;
            Y_rsop =  oo21(b).*X_pmd+oo22(b).*Y_pmd;
        end
        X_rsop_f = fftshift(fft(X_rsop)); % FFT并调整频率顺序
        Y_rsop_f = fftshift(fft(Y_rsop));
        X_pmd_f = zeros(size(X_rsop_f));
        Y_pmd_f = zeros(size(Y_rsop_f));
        for k = 1:rr
            phase = 2*pi*Tau(b)*fc;
            B(:,:,k) = [exp(-1i*phase/2),0;0,exp(1i*phase/2)]; % PMD矩阵
        end

        for k = 1:rr
            Bk = B(:,:,k);
            sig_in = [X_rsop_f(k);Y_rsop_f(k)];
            sig_out = Bk * sig_in; % 矩阵乘法
            X_pmd_f(k) = sig_out(1);
            Y_pmd_f(k) = sig_out(2);
        end
        % 转换回时域
        X_pmd = ifft(ifftshift(X_pmd_f)); % 逆FFT并调整频率顺序
        Y_pmd = ifft(ifftshift(Y_pmd_f));
    end
    Eoutx((a-1)*rr+1:a*rr) = X_pmd;
    Eouty((a-1)*rr+1:a*rr) = Y_pmd;
end
Eoutx = Eoutx';
Eouty = Eouty';
Eout = [Eoutx Eouty];
figure;
plot(real(Eout(:,1)), imag(Eout(:,1)),'.');
hold on;
plot(real(Eout(:,2)), imag(Eout(:,2)),'.');
title('经过rsop+pmd后的星座图');
legend('X-Pol', 'Y-Pol');
axis equal;
axis([-2.1,2.1,-2.1,2.1]);



% %================ 匹配滤波 =======================%
sig_hps_x = conv(sig_ps_x, filterCoeffs2, 'same');
sig_hps_x = (sig_hps_x-mean(sig_hps_x))./sqrt(mean(abs(sig_hps_x).^2));
sig_hps_y = conv(sig_ps_y, filterCoeffs2, 'same');
sig_hps_y = (sig_hps_y-mean(sig_hps_y))./sqrt(mean(abs(sig_hps_y).^2));
sig_hps = [sig_hps_x, sig_hps_y];
%--------------------------------------------------------%
figure;
plot(real(sig_hps_x), imag(sig_hps_x),'.');
hold on;
plot(real(sig_hps_y), imag(sig_hps_y),'.');
title('RRC');
legend('X-Pol', 'Y-Pol');
%% downsample
signal_dp = downsample(sig_hps, sps);
signal_dp_x = (signal_dp(:,1)-mean(signal_dp(:,1)))./sqrt(mean(abs(signal_dp(:,1)).^2));
signal_dp_y = (signal_dp(:,2)-mean(signal_dp(:,2)))./sqrt(mean(abs(signal_dp(:,2)).^2));
% --------------------------------------------------------%
figure;
plot(real(signal_dp_x), imag(signal_dp_x),'.');
hold on;
plot(real(signal_dp_y), imag(signal_dp_y),'.');
legend('X-Pol', 'Y-Pol');
axis square;
axis([-1.1,1.1,-1.1,1.1]);
title('downsample');




%% 参数设定 %%

L = length(signal_dp_x);
tap = 41;
step = 4;
rx = signal_dp_x;
ry = signal_dp_y;
v = [1e-5 1e-5 1e-5 2e-5 2e-5 2e-5];
Q = diag(v);
H = zeros(3);
vv = [1e-2 1e-2];
R = diag(vv);
x = [0.0000001,0,0,0,0,0]';
pp = [1e-4 1e-4 1e-4 2e-5 2e-5 2e-5];
P = diag(pp);
I = eye(6);
z = [0,0]';
ux = zeros(tap,1);
uy = zeros(tap,1);
x1 = cell(L,1);
hk = cell(L,1);

d = 1;
%% 滤波算法 %%
for i = tap:step:L
    %% 预测
    x1{d} = x;
    P1 = P + Q;
    in_x = rx(i-tap+1:1:i);
    in_y = ry(i-tap+1:1:i);
    f_inx = fft(in_x);
    f_iny = fft(in_y);
    d_tau = sqrt(x(1)^2+x(2)^2+x(3)^2)*1e-12;
    t1 = x(1)*1e-12;
    t2 = x(2)*1e-12;
    t3 = x(3)*1e-12;
    ksi = x(4);
    eta = x(5);
    kappa = x(6);
    uu11 = cos(pi*fc*d_tau)-1j*t1*sin(pi*fc*d_tau)/d_tau;
    uu12 = -t3*sin(pi*fc*d_tau)/d_tau-1j*t2*sin(pi*fc*d_tau)/d_tau;
    uu21 = -conj(uu12);
    uu22 = conj(uu11);
    U_pmd = [uu11 uu12;uu21 uu22];
    f_pmdx = uu11.*f_inx+uu12.*f_iny;
    f_pmdy = uu21.*f_inx+uu22.*f_iny;
    outx = ifft(f_pmdx);
    outy = ifft(f_pmdy);
    J = [cos(kappa)*exp(-1j*ksi),sin(kappa)*exp(1j*eta);-sin(kappa)*exp(-1j*eta),cos(kappa)*exp(1j*ksi)];
    ux = J(1,1).*outx+J(1,2).*outy;
    uy = J(2,1).*outx+J(2,2).*outy;
    xout{d}=ux(1:step);
    yout{d}=uy(1:step);
    uxx = mean(ux);
    uyy = mean(uy);
    hk{d} = [1-uxx*conj(uxx);1-uyy*conj(uyy)];
    inxx = mean(in_x);
    inyy = mean(in_y);
    Eox = uu11*inxx+uu12*inyy;
    Eoy = uu21*inxx+uu22*inyy;

    cc1 = pi*fc*t1*sin(pi*fc*d_tau)/d_tau;
    cc2 = ((t2^2+t3^2)*sin(pi*fc*d_tau)+pi*fc*t1^2*d_tau*cos(pi*fc*d_tau))/d_tau^3;
    cc3 = (t1*pi*fc*d_tau*cos(pi*fc*d_tau)-t1*sin(pi*fc*d_tau))/d_tau^3;
    xt1_a = J(1,1)*(-cc1-1j*cc2)+J(2,1)*(-1j*t2-t3)*cc3;
    xt1_b = J(1,2)*(-cc1-1j*cc2)+J(2,2)*(-1j*t2-t3)*cc3;
    xt1 =(xt1_a*inxx+ xt1_b*inyy)*conj(uxx)+uxx*conj(xt1_a*inxx+ xt1_b*inyy);
    dd1 = pi*fc*t2*sin(pi*fc*d_tau)/d_tau;
    dd2 = t1*t2*(pi*fc*d_tau*cos(pi*fc*d_tau)-sin(pi*fc*d_tau))/d_tau^3;
    dd3 = ((t1^2+t3^2)*sin(pi*fc*d_tau)+pi*fc*t2^2*d_tau*cos(pi*fc*d_tau))/(d_tau)^3;
    dd4 = t3*t2*(pi*fc*d_tau*cos(pi*fc*d_tau)-sin(pi*fc*d_tau))/d_tau^3;
    xt2_a = J(1,1)*(-dd1-1j*dd2)+J(2,1)*(-1j*dd3-dd4);
    xt2_b = J(1,2)*(-dd1-1j*dd2)+J(2,2)*(-1j*dd3-dd4);
    xt2 = (xt2_a*inxx+ xt2_b*inyy)*conj(uxx)+uxx*conj(xt2_a*inxx+ xt2_b*inyy);
    ee1 = pi*fc*t3*sin(pi*fc*d_tau)/d_tau;
    ee2 = t1*t3*(pi*fc*d_tau*cos(pi*fc*d_tau)-sin(pi*fc*d_tau))/d_tau^3;
    ee3 = t3*t2*(pi*fc*d_tau*cos(pi*fc*d_tau)-sin(pi*fc*d_tau))/d_tau^3;
    ee4 = ((t1^2+t2^2)*sin(pi*fc*d_tau)+pi*fc*t3^2*d_tau*cos(pi*fc*d_tau))/(d_tau)^3;
    xt3_a = J(1,1)*(-ee1-1j*ee2)+J(2,1)*(-1j*ee3-ee4);
    xt3_b = J(1,2)*(-ee1-1j*ee2)+J(2,2)*(-1j*ee3-ee4);
    xt3 = (xt3_a*inxx+ xt3_b*inyy)*conj(uxx)+uxx*conj(xt3_a*inxx+ xt3_b*inyy);
    x_kappa1 = -sin(kappa)*exp(-1j*ksi)*uu11-cos(kappa)*exp(-1j*eta)*uu12;
    x_kappa2 = cos(kappa)*exp(1j*eta)*uu11-sin(kappa)*exp(1j*ksi)*uu12;
    x_kappa = (x_kappa1*inxx+x_kappa2*inyy)*conj(uxx)+conj(x_kappa1*inxx+x_kappa2*inyy)*uxx;

    x_ksi1 = -cos(kappa)*exp(-1j*ksi)*uu11;
    x_ksi2 = cos(kappa)*exp(1j*ksi)*uu12;
    x_ksi = (x_ksi1*inxx+x_ksi2*inyy)*conj(uxx)+conj(x_ksi1*inxx+x_ksi2*inyy)*uxx;

    x_eta1 = sin(kappa)*exp(-1j*eta)*uu11;
    x_eta2 = sin(kappa)*exp(1j*eta)*uu12;
    x_eta = (x_eta1*inxx+ x_eta2*inyy)*conj(uxx)+conj(x_eta1*inxx+ x_eta2*inyy)*uxx;


    %ff1 = (t1*pi*fc*d_tau*cos(pi*fc*d_tau)-t1*sin(pi*fc*d_tau))/d_tau^3;
    ff1 = cc3;
    %ff2 = pi*fc*t1*sin(pi*fc*d_tau)/d_tau;
    ff2 = cc1;
    ff3 = ((t2^2+t3^2)*sin(pi*fc*d_tau)+pi*fc*t1^2*d_tau*cos(pi*fc*d_tau))/(d_tau^3);

    %ff3 = cc2;
    yt1_a = J(1,1)*(-1j*t2+t3)*ff1+J(2,1)*(-ff2+1j*ff3);
    yt1_b = J(1,2)*(-1j*t2+t3)*ff1+J(2,2)*(-ff2+1j*ff3);
    yt1 = (yt1_a*inxx+yt1_b*inyy)*conj(uyy)+uyy*conj(yt1_a*inxx+yt1_b*inyy);
    gg1 = dd3;
    gg2 = dd4;
    gg3 = dd1;
    gg4 = dd2;
    yt2_a = J(1,1)*(-1j*gg1+gg2)+J(2,1)*(-gg3+1j*gg4);
    yt2_b = J(1,2)*(-1j*gg1+gg2)+J(2,2)*(-gg3+1j*gg4);
    yt2 = (yt2_a*inxx+yt2_b*inyy)*conj(uyy)+uyy*conj(yt2_a*inxx+yt2_b*inyy);
    hh1 = ee3;
    hh2 = ee4;
    hh3 = ee1;
    hh4 = ee2;
    yt3_a = J(1,1)*(-1j*hh1+hh2)+J(2,1)*(-hh3+1j*hh4);
    yt3_b = J(1,2)*(-1j*hh1+hh2)+J(2,2)*(-hh3+1j*hh4);
    yt3 = (yt3_a*inxx+yt3_b*inyy)*conj(uyy)+uyy*conj(yt3_a*inxx+yt3_b*inyy);

    y_kappa1 = -sin(kappa)*exp(-1j*ksi)*uu21-cos(kappa)*exp(-1j*eta)*uu22;
    y_kappa2 = cos(kappa)*exp(1j*eta)*uu21-sin(kappa)*exp(1j*ksi)*uu22;
    y_kappa = (x_kappa1*inxx+x_kappa2*inyy)*conj(uyy)+conj(x_kappa1*inxx+x_kappa2*inyy)*uyy;

    y_ksi1 = -cos(kappa)*exp(-1j*ksi)*uu21;
    y_ksi2 = cos(kappa)*exp(1j*ksi)*uu22;
    y_ksi = (x_ksi1*inxx+x_ksi2*inyy)*conj(uyy)+conj(x_ksi1*inxx+x_ksi2*inyy)*uyy;

    y_eta1 = sin(kappa)*exp(-1j*eta)*uu21;
    y_eta2 = sin(kappa)*exp(1j*eta)*uu22;
    y_eta = (x_eta1*inxx+ x_eta2*inyy)*conj(uyy)+conj(x_eta1*inxx+ x_eta2*inyy)*uyy;

    xt1 = xt1*1e-12;
    xt2 = xt2*1e-12;
    xt3 = xt3*1e-12;
    yt1 = yt1*1e-12;
    yt2 = yt2*1e-12;
    yt3 = yt3*1e-12;

    H = [xt1 xt2 xt3  x_ksi x_eta x_kappa ;yt1  yt2  yt3  y_ksi y_eta y_kappa];
    K = P1*H'*((H*P1*H' + R))^-1; 
    x = x1{d} + K* hk{d};
    P = (I - K*H)*P1;
    d = d+1;
end
matrix_x = [xout{:}];
col_xout = matrix_x(:);
matrix_y = [yout{:}];
col_yout = matrix_y(:);

figure;
plot(real(col_xout), imag(col_xout),'.');
hold on;
plot(real(col_yout), imag(col_yout),'.');
title('kalma');
legend('X-Pol', 'Y-Pol');
%%

P = qammod(0:M-1,M).';
% 归一化
P=pnorm(P);

sig_LMS_xi = (col_xout-mean(col_xout))./sqrt(mean(abs(col_xout).^2));
sig_LMS_yi = (col_yout-mean(col_yout))./sqrt(mean(abs(col_yout).^2));
[sig_LMS_xo, sig_LMS_yo] = DDLMS_v0(sig_LMS_xi,sig_LMS_yi,mu_DDLMS,ntaps_DDLMS,1,P);
figure;
plot(real(sig_LMS_xo), imag(sig_LMS_xo),'.');
hold on;
plot(real(sig_LMS_yo), imag(sig_LMS_yo),'.');
title('DD-LMS');
legend('X-Pol', 'Y-Pol')

%%
rxSignal = hard_decision_qam(M,sig_LMS_xo);
rxSignal = rxSignal(:);
rySignal = hard_decision_qam(M,sig_LMS_yo);
rySignal = rySignal(:);
%%
x = rxSignal;
y = rySignal;

ref_seq_1 = data(:,1);
ref_seq_2 = data(:,2);
swap_options = [false, true];
rot_angles = [0, 90, 180, 270];
conj_options = [false, true];
min_ber = Inf;
best_swap = false;
best_rot = 0;
best_conj = false;
m=0;

for swap = swap_options
    for rot_angle = rot_angles
        for conj_flag = conj_options
            % 应用交换
            if swap
                current_x = y;
                current_y = x;
            else
                current_x = x;
                current_y = y;
            end

            % 应用旋转因子
            switch rot_angle
                case 0
                    rot_factor = 1;
                case 90
                    rot_factor = 1i;
                case 180
                    rot_factor = -1;
                case 270
                    rot_factor = -1i;
            end
            current_x = current_x * rot_factor;
            current_y = current_y * rot_factor;

            % 应用复共轭
            if conj_flag
                current_x = conj(current_x);
                current_y = conj(current_y);
            end
            xxx = qamdemod(current_x, M , 'UnitAveragePower', true);
            xxx = xxx(:);
            yyy = qamdemod(current_y, M , 'UnitAveragePower', true);
            yyy = yyy(:);
            [xxx,ref_seq1] = sync(xxx,ref_seq_1);
            [yyy,ref_seq2] = sync(yyy,ref_seq_2);
            % 计算误码率（符号错误率）
            [errors_x,num1,~] = CalcBER(xxx,ref_seq1); %计算误码率
            [errors_y,num2,~] = CalcBER(yyy,ref_seq2); %计算误码率
            total_errors = errors_x + errors_y;
            ber = total_errors/2 ; % 总符号数是2N

            % 更新最佳方案
            if ber < min_ber
                min_ber = ber;
                best_swap = swap;
                best_rot = rot_angle;
                best_conj = conj_flag;
            end
            m=m+1;
        end
    end
end

fprintf('最小误码率: %.4f\n', min_ber);
fprintf('最佳变换方案:\n');
fprintf('交换x和y: %d\n', best_swap);
fprintf('旋转角度: %d°\n', best_rot);
fprintf('取复共轭: %d\n', best_conj);




%% phase compensation
% NN = 100;
% X = ux;
% Y = uy;
% MM = floor(length(X)/NN);
%
%
%
% for m= 1:MM
%    x_prime = X((m-1)*NN+1:m*NN);
%    V_pxr = real(x_prime);
%    V_pxi = imag(x_prime);
%    V_p2xr =  V_pxr.^2-V_pxi.^2;
%    V_p2xi =  2*(V_pxr.*V_pxi);
%    V_p4xr = V_p2xr.^2-V_p2xi.^2;
%    V_p4xi = 2*V_p2xr.*V_p2xi;
%    mean_V_p4xr = sum(V_p4xr);
%    mean_V_p4xi = sum(V_p4xi);
%  y
%    y_prime = Y((m-1)*NN+1:m*NN);
%    V_pyr = real(y_prime);
%    V_pyi = imag(y_prime);
%    V_p2yr =  V_pyr.^2-V_pyi.^2;
%    V_p2yi =  2*(V_pyr.*V_pyi);
%    V_p4yr = V_p2yr.^2-V_p2yi.^2;
%    V_p4yi = 2*V_p2yr.*V_p2yi;
%    mean_V_p4yr = sum(V_p4yr);
%    mean_V_p4yi = sum(V_p4yi);
%
%
% PEx14 = atan(mean_V_p4xi/mean_V_p4xr);
% PEy14 = atan(mean_V_p4yi/mean_V_p4yr);
%
% if (mean_V_p4xr <0 && mean_V_p4xi <0)
%    PEx4 = PEx14-pi;
% end
% if (mean_V_p4xr <0 && mean_V_p4xi >0)
%    PEx4 = PEx14+pi;
% end
% if (mean_V_p4xr >0 )
%    PEx4 = PEx14;
% end
% if (mean_V_p4xr == 0 )
%     if(mean_V_p4xi >0)
%     PEx4 = pi/2;
%     elseif(mean_V_p4xi <0)
%      PEx4 = -pi/2;
%     end
% end
%
% if (mean_V_p4yr <0 && mean_V_p4yi <0)
%    PEy4 = PEy14-pi;
% end
% if (mean_V_p4yr <0 && mean_V_p4yi >0)
%    PEy4 = PEy14+pi;
% end
% if (mean_V_p4yr >0 )
%    PEy4 = PEy14;
% end
% if (mean_V_p4yr == 0 )
%     if(mean_V_p4yi >0)
%     PEy4 = pi/2;
%     elseif(mean_V_p4yi <0)
%      PEy4 = -pi/2;
%     end
% end
%
% PEx = PEx4/4;
% PEy = PEy4/4;
%
%
% out_xr((m-1)*NN+1:m*NN) = -V_pxr.*cos(PEx)-V_pxi.*sin(PEx);
% out_xi((m-1)*NN+1:m*NN) = -V_pxi.*cos(PEx)+V_pxr.*sin(PEx);
% out_yr((m-1)*NN+1:m*NN) = -V_pyr.*cos(PEy)-V_pyi.*sin(PEy);
% out_yi((m-1)*NN+1:m*NN) = -V_pyi.*cos(PEy)+V_pyr.*sin(PEy);
% end
% figure;
% plot(out_xr,out_xi,'.');
% hold on;
% plot(out_yr,out_yi,'.');
% title('v-v相位补偿');
% legend('X-Pol', 'Y-Pol');
%
% eta = pi/4;
% eta1 = pi/4;
% out_xr1 = cos(eta).*out_xr-sin(eta).*out_xi;
% out_xi1 = sin(eta).*out_xr+cos(eta).*out_xi;
% out_yr1 = cos(eta).*out_yr-sin(eta).*out_yi;
% out_yi1 = sin(eta).*out_yr+cos(eta).*out_yi;
%
%
%
% figure;
% plot(out_xr1,out_xi1,'.');
% hold on;
% plot(out_yr1,out_yi1,'.');
% title('45度后');
% legend('X-Pol', 'Y-Pol');




