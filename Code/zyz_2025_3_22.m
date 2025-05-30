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
sps = 2;
roll_off = 0.1;
% % RSOP set
% r = 33e-3;
% psi = rand*pi/2;
% chi = rand*pi/2;
% CMA set
mu_CMA = 4e-3;
ntaps =9;
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
figure;
plot(real(sig_resample(:,1)), imag(sig_resample(:,1)),'.','MarkerSize',15);
hold on;
plot(real(sig_resample(:,2)), imag(sig_resample(:,2)),'.','MarkerSize',15);
title('上采样');
legend('X-Pol', 'Y-Pol');
% ---------------------pulse shaping-------------------%

filterSymbolLength = 31; % 2*10*16+1
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

% ---------------------RSOP-------------------%
t_spacing = 5e-11;
T = length(sig_ps);
t = 0:t_spacing:(T-1)*t_spacing;
N =12;
w_alpha = 5e3; %krad/s---rad/s
w_phi = 13e3;
w_kappa =11e3;
alpha= zeros(N,1);
phi = zeros(N,1);
kappa = zeros(N,1);
U = cell(1,T);

kappa0 = pi+randn(N,1)*2*pi;
alpha0 = pi+randn(N,1)*2*pi;
phi0 = pi+randn(N,1)*2*pi;

for a = 1:T
    alpha = w_alpha.*t(a)+alpha0;
    phi = w_phi.*t(a)+phi0;
    kappa = w_kappa.*t(a)+kappa0;
    oo11 = cos(kappa).*exp(1j*(alpha));
    oo12 = -sin(kappa).*exp(-1j*(phi));
    oo21 = -conj(oo12);
    oo22 = conj(oo11);
    y=1;
    for b = 1:N
        y = y*[oo11(b) oo12(b);oo21(b) oo22(b)];
    end
    U{a} = y;
sig_chx_out(a) = U{a}(1,1).*sig_ps_x(a)+U{a}(1,2).*sig_ps_y(a);
sig_chy_out(a) = U{a}(2,1).*sig_ps_x(a)+U{a}(2,2).*sig_ps_y(a);
end
sig_chx_out = sig_chx_out(:);

sig_chy_out = sig_chy_out(:);

% --------------------------------------------------------%
figure;
plot(real(sig_chx_out),imag(sig_chx_out),'.');
hold on;
plot(real(sig_chy_out),imag(sig_chy_out),'.');
title('RSOP');
legend('X-Pol', 'Y-Pol');
% %================ 匹配滤波 =======================%

sig_hps_x = conv(sig_chx_out, filterCoeffs2, 'same');
sig_hps_x=pnorm(sig_hps_x);
sig_hps_y = conv(sig_chy_out, filterCoeffs2, 'same');
sig_hps_y=pnorm(sig_hps_y);
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
% ==========================Rx=========================%
%%

L = length(signal_dp_x);
rx = signal_dp_x;
ry = signal_dp_y;
v = [1e-3 1e-3 1e-3 1e-3];
Q = diag(v);
H = zeros(3);
vv = [0.1 0.1];
R = diag(vv); 
x = [1,0,0,0]';
pp = [2e-5 2e-5 2e-5 2e-5];
P = diag(pp);
I = eye(4);
z = [0,0]';
ux = zeros(L,1);
uy = zeros(L,1);
x1 = cell(L,1);
hk = cell(L,1);

for m = 1:L
% 预测
x1{m} = x;

J = [x(1)-1j*x(2),x(3)+1j*x(4);-x(3)+1j*x(4),x(1)+1j*x(2)];
ux(m)=J(1,1)*rx(m)+J(1,2)*ry(m);
uy(m)=J(2,1)*rx(m)+J(2,2)*ry(m);
hk{m} = [1-ux(m)*conj(ux(m));1-uy(m)*conj(uy(m))];
%hk{m} = [ux(m)*conj(ux(m))-1;uy(m)*conj(uy(m))-1];
H11 = 2*real(rx(m)*conj(ux(m)));
H12 = 2*imag(rx(m)*conj(ux(m)));                           
H13 = 2*real(ry(m)*conj(ux(m)));
H14 = -2*imag(ry(m)*conj(ux(m)));
H21 = 2*real(ry(m)*conj(uy(m)));
H22 = -2*imag(ry(m)*conj(uy(m)));
H23 = -2*real(rx(m)*conj(uy(m)));
H24 = -2*imag(rx(m)*conj(uy(m)));
H = [H11 H12 H13 H14;H21 H22 H23 H24];
P1 = P + Q;
K = P1*H'*((H*P1*H' + R))^-1;
x = x1{m} + K* hk{m};
P = (I - K*H)*P1;
end
 ux = ux(end/2:end);
 uy = uy(end/2:end);
figure;
plot(real(ux), imag(ux),'.');
hold on;
plot(real(uy), imag(uy),'.');
title('kalma');
legend('X-Pol', 'Y-Pol');



% %================ DD-LMS ===========================%


% P = qammod(0:M-1,M).';
% % 归一化
% P=pnorm(P);
% 
% sig_LMS_xi = (ux-mean(ux))./sqrt(mean(abs(ux).^2));
% sig_LMS_yi = (uy-mean(uy))./sqrt(mean(abs(uy).^2));
% [sig_LMS_xo, sig_LMS_yo] = DDLMS_v0(sig_LMS_xi,sig_LMS_yi,mu_DDLMS,ntaps_DDLMS,1,P);
% figure;
% plot(real(sig_LMS_xo), imag(sig_LMS_xo),'.');
% hold on;
% plot(real(sig_LMS_yo), imag(sig_LMS_yo),'.');
% title('DD-LMS');
% legend('X-Pol', 'Y-Pol')
%%
% NN = 10;
% X = sig_cma_xo;
% Y = sig_cma_yo;
% MM = floor(length(X)/NN);
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
%  %y
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
% out_xr((m-1)*NN+1:m*NN) = V_pxr.*cos(PEx)+V_pxi.*sin(PEx);
% out_xi((m-1)*NN+1:m*NN) = V_pxi.*cos(PEx)-V_pxr.*sin(PEx);
% out_yr((m-1)*NN+1:m*NN) = V_pyr.*cos(PEy)+V_pyi.*sin(PEy);
% out_yi((m-1)*NN+1:m*NN) = V_pyi.*cos(PEy)-V_pyr.*sin(PEy);
% end
% figure;
% plot(out_xr(2000:end),out_xi(2000:end),'.');
% hold on;       
% plot(out_yr(2000:end),out_yi(2000:end),'.');
% title('v-v相位补偿');
% legend('X-Pol', 'Y-Pol');
% 
% eta = -pi/4;
% %eta1 = pi/4;
% out_xr1 = cos(eta).*out_xr-sin(eta).*out_xi;
% out_xi1 = sin(eta).*out_xr+cos(eta).*out_xi;
% out_yr1 = cos(eta).*out_yr-sin(eta).*out_yi;
% out_yi1 = sin(eta).*out_yr+cos(eta).*out_yi;
% out_xr1 = out_xr1(4000:end);
% out_xi1 = out_xi1(4000:end);
% out_yr1 = out_yr1(4000:end);
% out_yi1 = out_yi1(4000:end);
% figure;
% plot(out_xr1,out_xi1,'.');
% hold on;
% plot(out_yr1,out_yi1,'.');
% title('45度后');
% legend('X-Pol', 'Y-Pol');
% out_x = out_xr1+1j*out_xi1;
% out_y = out_yr1+1j*out_yi1;


%%
rxSignal = hard_decision_qam(M,ux);
rxSignal = rxSignal(:);
rySignal = hard_decision_qam(M,uy);
rySignal = rySignal(:);
%%
x = rxSignal;
y = rySignal;

ref_seq_1 = data(end/2:end,1);
ref_seq_2 = data(end/2:end,2);
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








%%





% 参考序列


% 全部的序列进行解码

% [ber1,num1,error_location] = CalcBER(xxx,ref_seq_1); %计算误码率

%fprintf('Num of Errors = %d, SER = %1.7f\n',num1,ber1);
%%












