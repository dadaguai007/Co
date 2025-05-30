% 如何调用
% example:fcnorm=0.7*Rs/(fs/2);
%filt.H(f/fs)
function filt=Besslf_filter(order,fcnorm,verbose)
% fs=20e9;
% f3dB=5e9;
% order=5;
% %归一化的截止频率
% fcnorm = f3dB/(fs/2);
% verbose=1;
type='Besslf';
% Variables used when calculating impulse response
maxMemoryLength = 2^10; % maximum memory length
threshold = 1-1e-6;      % energy threshold above which impulse respone is neglected

wc2wo = 1.621597623980423;
wow = wc2wo*(2*pi*fcnorm/2);
% anaglo
[nums, dens] = besself(order, wow); % design prototype in CT with frequency prewarped
%         [num, den] = bilinear(nums, dens, 1, fcnorm/2);
% numels
[num, den] = impinvar(nums, dens, 1);
if verbose
%     verbose = false;
    f = linspace(0, 1/2);
    Hs = freqs(nums, dens, 2*pi*f);
    Hz = freqz(num, den, f, 1);
    plot_transform_num_analog(Hs, Hz, f, fcnorm/2)
end
filt.type = type;
filt.order = order;
filt.num = num;
filt.den = den;
% 确定滤波器的延时值
filt.grpdelay = grpdelay(num, den, 1);
filt.fcnorm = fcnorm;
% 滤波器的频率响应，去掉了延迟
filt.H = @(f) freqz(num, den, 2*pi*f).*exp(1j*2*pi*f*filt.grpdelay);
filt.H_test = @(f) freqz(num, den, 2*pi*f);
% 噪声，一般情况下不使用
filt.noisebw = @(fs) noisebw(num, den, 2^15, fs); % equivalent two-sided noise bandwidth over larger number of samples (2^15) given a sampling frequency fs

% 因为是从频域定义滤波器的，可以将频域转换为时域的脉冲响应使用
% 也就是得到impulse response
if den == 1 % FIR
    filt.h = num;
else
    x = zeros(1, maxMemoryLength+1);
    x(1) = 1;
    y = filter(filt.num, filt.den, x);
    E = cumsum(abs(y).^2)/sum(abs(y).^2);
    y(E > threshold) = [];
    y = y/abs(sum(y)); % normalize to have unit gain at DC
    %滤波器的时域响应
    filt.h = y;
end

if verbose
    plot_filter(filt)
end


%     function plot_filter(filt)
%         f = linspace(0, 0.5);
%         Hz = filt.H(f);
%         figure(111)
%         subplot(211), hold on, box on
%         plot(f, 20*log10(abs(Hz)), 'DisplayName', filt.type)
%         aa = axis;
%         h = plot([0 filt.fcnorm/2], [-3 -3], ':k');
%         hasbehavior(h, 'legend', false);   % line will not be in legend
%         h = plot(filt.fcnorm/2*[1 1], [aa(3) -3], ':k');
%         hasbehavior(h, 'legend', false);   % line will not be in legend
%         xlabel('Normalized frequency f/f_s')
%         ylabel('Magnitude (dB)')
%         legend('-DynamicLegend')
%         axis([0 0.5 -50 10])
% 
%         subplot(212), hold on, box on
%         plot(f, rad2deg(unwrap(angle(Hz))), 'DisplayName', filt.type)
%         xlabel('Normalized frequency f/f_s')
%         ylabel('Phase (deg)')
%         legend('-DynamicLegend')
%         axis([0 0.5 -360 360])
%     end
% 
%     function plot_transform(Hs, Hz, f, fc)
%         figure(113)
%         subplot(211), hold on, box on
%         plot(f, 20*log10(abs(Hs)), f, 20*log10(abs(Hz)))
%         aa = axis;
%         plot([0 fc], [-3 -3], ':k')
%         plot(fc*[1 1], [aa(3) -3], ':k')
%         xlabel('Normalized frequency f/f_s')
%         ylabel('Magnitude (dB)')
%         legend('Analog filter', 'Bilinear transform')
%         axis([0 0.5 -50 10])
% 
%         subplot(212), hold on, box on
%         plot(f, rad2deg(unwrap(angle(Hs))), f, rad2deg(unwrap(angle(Hz))))
%         xlabel('Normalized frequency f/f_s')
%         ylabel('Phase (deg)')
%         legend('Analog filter', 'Bilinear transform')
%         axis([0 0.5 -360 360])
%     end

end
