function filt=Gaussian_filter(order,fcnorm,nlength,verbose)
% order=5;
% fs=20e9;
% f3dB=5e9;
% fcnorm = f3dB/(fs/2);
% verbose=1;
type='gaussian';
% nlength = 256;   % number of samples used for fitting.
%nlength 为信号的样本数

% Gaussian filter
df = 1/nlength;
f = -0.5:df:0.5-df;
f0 = (fcnorm)/(2*log(2)^(1/(2*order)));
H = exp(-0.5*(f/f0).^(2*order));
%从频域得到时域的抽头数
h = real(fftshift(ifft(ifftshift(H))));

% 配置为FIR滤波器
den = 1;
num = h;

% 一般不需要使用
% Noise bandwidth
nbw = @(fs) sum(abs(num).^2)*fs; % requires unit gain at DC

if verbose
    %     verbose = false;
    f = linspace(0, 1/2);
    % 模拟滤波器响应
    Hs = exp(-0.5*(f/f0).^(2*order));
    %数字滤波器响应
    Hz = freqz(num, den, f, 1);
    plot_approx(Hs, Hz, f, fcnorm/2)
end

filt.type = type;
filt.order = order;
filt.num = num;
filt.den = den;
filt.grpdelay = grpdelay(num, den, 1);
filt.fcnorm = fcnorm;
filt.H = @(f) freqz(num, den, 2*pi*f).*exp(1j*2*pi*f*filt.grpdelay);
filt.noisebw = nbw;


if den == 1 % FIR
    filt.h = num;
end

if verbose
    plot_filter(filt)
end

    function plot_approx(Hs, Hz, f, fc)
        figure(112)
        subplot(211), hold on, box on
        plot(f, 20*log10(abs(Hs)), f, 20*log10(abs(Hz)), '--')
        aa = axis;
        plot([0 fc], [-3 -3], ':k')
        plot(fc*[1 1], [aa(3) -3], ':k')
        xlabel('Normalized frequency f/f_s')
        ylabel('Magnitude (dB)')
        legend('Analog filter', 'Approximation')
        axis([0 0.5 -50 10])

        subplot(212), hold on, box on
        plot(f, rad2deg(unwrap(angle(Hs))), f, rad2deg(unwrap(angle(Hz))), '--')
        xlabel('Normalized frequency f/f_s')
        ylabel('Phase (deg)')
        legend('Analog filter', 'FIR approximation')
        axis([0 0.5 -360 360])
    end
end