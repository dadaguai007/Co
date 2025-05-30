function plot_transform_num_analog(Hs, Hz, f, fc)
% 数字和模拟滤波器的对比
figure(113)
subplot(211), hold on, box on
plot(f, 20*log10(abs(Hs)), f, 20*log10(abs(Hz)))
aa = axis;
plot([0 fc], [-3 -3], ':k')
plot(fc*[1 1], [aa(3) -3], ':k')
xlabel('Normalized frequency f/f_s')
ylabel('Magnitude (dB)')
legend('Analog filter', 'Bilinear transform')
axis([0 0.5 -50 10])

subplot(212), hold on, box on
plot(f, rad2deg(unwrap(angle(Hs))), f, rad2deg(unwrap(angle(Hz))))
xlabel('Normalized frequency f/f_s')
ylabel('Phase (deg)')
legend('Analog filter', 'Bilinear transform')
axis([0 0.5 -360 360])
end