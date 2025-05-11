function plotEquParam(xn,sigRx_E,label,w,en)

% 均衡后信号时域分布
figure;
subplot(3,1,1)
hold on;
plot(xn(1:1e5),'.')
plot(sigRx_E(1:1e5),'k.')
plot(label(1:1e5),'.')
legend('接收信号','均衡后信号','发送信号')

% 抽头数
subplot(3,1,2)
stem(w(:))

% 误差是否收敛
subplot(3,1,3)
semilogy(abs(en(1:1e4)).^2)
xlabel("迭代次数")
ylabel("误差")

end