function y = delay_signal(x, delay)
%% Delay signal by delay samples, which need not be integer

transpose = false;
if size(x, 1) > size(x, 2)
    x = x.';
    transpose = true;
end

if is_Integer(delay) % just circshift signal
    % time domain delay way
    delay = round(delay);
    y = circshift(x, [0 delay]);
else
    % frequence dommain delay way
    df = 1/length(x);
    f = ifftshift(-0.5:df:0.5-df);    
    y = ifft(fft(x).*exp(-1j*2*pi*f*delay));
end

if isreal(x)
    y = real(y);
end

if transpose
    y = y.';
end
    %整数进行circshishift，不是整数使用频域延迟方法