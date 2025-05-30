function plot_filter(filt)
% filt is the 结构体
% 绘制数字滤波器的响应曲线
% 并给出相应的带宽对应点
f = linspace(0, 0.5);
Hz = filt.H(f);
figure(111)
subplot(211), hold on, box on
plot(f, 20*log10(abs(Hz)), 'DisplayName', filt.type)
aa = axis;
h = plot([0 filt.fcnorm/2], [-3 -3], ':k');
hasbehavior(h, 'legend', false);   % line will not be in legend
h = plot(filt.fcnorm/2*[1 1], [aa(3) -3], ':k');
hasbehavior(h, 'legend', false);   % line will not be in legend
xlabel('Normalized frequency f/f_s')
ylabel('Magnitude (dB)')
legend('-DynamicLegend')
axis([0 0.5 -50 10])

subplot(212), hold on, box on
plot(f, rad2deg(unwrap(angle(Hz))), 'DisplayName', filt.type)
xlabel('Normalized frequency f/f_s')
ylabel('Phase (deg)')
legend('-DynamicLegend')
axis([0 0.5 -360 360])
end
