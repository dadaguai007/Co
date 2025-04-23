EQ=struct();
EQ.u=0.01;
EQ.k1=31;
EQ.k2=15;
EQ.ref=8;
EQ.sps=2;
EQ.lamda=0.9999;
EQ.delta=0.01;
% FFE_LMS
% [sigRx_E,en,w] = FFE_LMS(EQ, xn, label);
% FFE_RLS
% [sigRx_E,en,w] = FFE_RLS(EQ, xn, label);
% DFE_LMS
[sigRx_E,en,w] = DFE_LMS(EQ, xn, label);
% DFE_RLS
% [sigRx_E,en,w] = DFE_RLS(EQ, xn, label);
% APA
% [sigRx_E,en,w]=APA(EQ, xn, label);
% % Volterra
sps=2;
ref=8;
taps_list = [31 11 0 15 0 0];
taps_list_ABS = [31 0 0 15 5 3];
step_len = 0.001;
lamda=0.9999;
% volterra_lms
% [sigRx_E,en,w]=volterra_dfe_lms(xn,label,sps,ref,taps_list,step_len);
% volterra_rls
% [sigRx_E,en,w]=volterra_dfe_rls(xn,label,sps,ref,taps_list,lamda);
% abf_rls
% [sigRx_E,en,w]=abf_lms(xn,label,sps,ref,taps_list_ABS,step_len);

% pwl
% [sigRx_E,en,w] = pwl(xn,label,sps,ref,15,0.005,[-0.7661,0,0.7661]);

figure;
subplot(3,1,1)
hold on;
plot(xn(1:1e5),'.')
plot(sigRx_E(1:1e5),'k.')
plot(label(1:1e5),'.')
legend('接收信号','均衡后信号','发送信号')

subplot(3,1,2)
stem(w(:))


subplot(3,1,3)
semilogy(abs(en(1:1e4)).^2)
xlabel("迭代次数")
ylabel("误差")