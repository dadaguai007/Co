% 卡尔曼滤波
% 状态向量 x = [τ1, τ2, τ3, κ, η, ζ]^T
% 其中：
% τ1,τ2,τ3: PMD矢量参数(ps)
% κ: 偏振旋转角
% η: 偏振椭圆角
% ζ: 偏振方位角

%% 参数设定
L = length(signal_dp_x);          % 信号长度
tap = 41;                        % 滤波器抽头数
step = 4;                        % 处理步长
rx = signal_dp_x;                % 接收X偏振信号
ry = signal_dp_y;                % 接收Y偏振信号
v = [1e-5 1e-5 1e-5 2e-5 2e-5 2e-5]; % 过程噪声方差
Q = diag(v);                     % 过程噪声协方差矩阵
H = zeros(2,6);                    % 测量矩阵(后续动态计算)
vv = [1e-2 1e-2];                % 测量噪声方差
R = diag(vv);                    % 测量噪声协方差矩阵
x = [0.0000001,0,0,0,0,0]';      % 初始状态估计
pp = [1e-4 1e-4 1e-4 2e-5 2e-5 2e-5]; % 初始状态协方差
P = diag(pp);                    % 状态协方差矩阵
I = eye(6);                      % 6维单位矩阵
z = [0,0]';                      % 测量向量
ux = zeros(tap,1);               % X偏振补偿信号缓存
uy = zeros(tap,1);               % Y偏振补偿信号缓存
x1 = cell(L,1);                  % 状态预测值存储
hk = cell(L,1);                  % 测量残差存储

d = 1;                           % 数据块索引



for i = tap:step:L % 滑动窗口处理
    %% 预测步骤
    x1{d} = x;                   % 存储当前状态预测值
    P1 = P + Q;                  % 状态协方差预测

    %% 信号预处理
    in_x = rx(i-tap+1:i);        % 提取X偏振当前窗口信号
    in_y = ry(i-tap+1:i);        % 提取Y偏振当前窗口信号
    f_inx = fft(in_x);           % X信号FFT
    f_iny = fft(in_y);           % Y信号FFT

    %% 状态参数提取
    d_tau = sqrt(x(1)^2+x(2)^2+x(3)^2)*1e-12; % PMD矢量幅值
    t1 = x(1)*1e-12;             % τ1分量(转换为秒)
    t2 = x(2)*1e-12;             % τ2分量
    t3 = x(3)*1e-12;             % τ3分量
    ksi = x(4);                  % 偏振参数κ
    eta = x(5);                  % 偏振参数η
    kappa = x(6);                % 偏振参数ζ

    %% PMD补偿矩阵构建
    uu11 = cos(pi*fc*d_tau)-1j*t1*sin(pi*fc*d_tau)/d_tau;
    uu12 = -t3*sin(pi*fc*d_tau)/d_tau-1j*t2*sin(pi*fc*d_tau)/d_tau;
    uu21 = -conj(uu12);
    uu22 = conj(uu11);
    U_pmd = [uu11 uu12; uu21 uu22]; % PMD补偿矩阵

    %% 频域PMD补偿
    f_pmdx = uu11.*f_inx + uu12.*f_iny; % X偏振频域补偿
    f_pmdy = uu21.*f_inx + uu22.*f_iny; % Y偏振频域补偿
    outx = ifft(f_pmdx);          % X时域信号恢复
    outy = ifft(f_pmdy);          % Y时域信号恢复

    %% 偏振旋转补偿
    J = [cos(kappa)*exp(-1j*ksi), sin(kappa)*exp(1j*eta);
        -sin(kappa)*exp(-1j*eta), cos(kappa)*exp(1j*ksi)]; % 旋转矩阵
    ux = J(1,1).*outx + J(1,2).*outy; % X偏振最终输出
    uy = J(2,1).*outx + J(2,2).*outy; % Y偏振最终输出

    %% 测量残差计算
    uxx = mean(ux);              % X偏振均值
    uyy = mean(uy);              % Y偏振均值
    hk{d} = [1-uxx*conj(uxx);    % 测量残差计算(理想QPSK模长应为1)
        1-uyy*conj(uyy)];    % h(x) = 1-|ux|²
    % 输入
    inxx = mean(in_x);
    inyy = mean(in_y);
    %% 雅可比矩阵计算（线性化关键步骤）
    % 详细偏导数计算（此处约300行代码）
    % 主要计算H矩阵: H = dh/dx
    % 包含对τ1,τ2,τ3,κ,η,ζ的偏导
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
    y_kappa = (y_kappa1*inxx+y_kappa2*inyy)*conj(uyy)+conj(y_kappa1*inxx+y_kappa2*inyy)*uyy;

    y_ksi1 = -cos(kappa)*exp(-1j*ksi)*uu21;
    y_ksi2 = cos(kappa)*exp(1j*ksi)*uu22;
    y_ksi = (y_ksi1*inxx+y_ksi2*inyy)*conj(uyy)+conj(y_ksi1*inxx+y_ksi2*inyy)*uyy;

    y_eta1 = sin(kappa)*exp(-1j*eta)*uu21;
    y_eta2 = sin(kappa)*exp(1j*eta)*uu22;
    y_eta = (y_eta1*inxx+ y_eta2*inyy)*conj(uyy)+conj(y_eta1*inxx+ y_eta2*inyy)*uyy;

    xt1 = xt1*1e-12;
    xt2 = xt2*1e-12;
    xt3 = xt3*1e-12;
    yt1 = yt1*1e-12;
    yt2 = yt2*1e-12;
    yt3 = yt3*1e-12;
    % H 矩阵生成，对h向量求偏导
    H = [xt1 xt2 xt3  x_ksi x_eta x_kappa ;yt1  yt2  yt3  y_ksi y_eta y_kappa];

    %% 卡尔曼增益计算
    K = P1*H'/(H*P1*H' + R);     % 卡尔曼增益公式
%     K = P1*H'*inv(H*P1*H' + R);   % 卡尔曼增益公式
    %% 状态更新
    x = x1{d} + K*hk{d};         % 状态估计更新
    P = (I - K*H)*P1;            % 协方差更新

    d = d+1;                     % 更新数据块索引
end
