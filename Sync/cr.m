function [f_clk] = cr(y,BaudRate,os,fig)
% [f_clk] = ClockRecovery(y,t,BaudRate,os)
% <Principle> : HPF -> square in time domain -> BPF
% <OUTPUT>
% f_clk : recovered clock frequency
% <INPUT>
% y : input signal
% t : time signal corresponding to 'y'
% BaudRate : Baud rate
% os : signal('y')'s oversampling rate


% HPF
D = fdesign.highpass('Fst,Fp,Ast,Ap',0.5/os/4,0.5/os/4+0.05,90,0.5);    %Fstop, Fpass, Attenuation in the stop band (dB), ripple in the pass band (dB)
Hd = design(D, 'equiripple');
y = filter(Hd, y);

y_freq = fftshift(fft(y));
df = os*BaudRate/length(y_freq);
if mod(length(y_freq),2) == 0   %when the length of sequence is an even number
    freq = [-length(y_freq)/2*df:df:-df, 0:df:df*(length(y_freq)/2-1)];
else                            %when the length of sequence is odd number
    freq = [-(length(y_freq)-1)/2*df:df:-df, 0:df:df*(length(y_freq)-1)/2];
end

if fig == 1
	plot(freq,abs(y_freq));
    title('Spectrum after HPF');
end



% Square in time domain
% format long
y=y-mean(y);        %%%%%%%%%%%%%%%%%%%
y = abs(y).^2;      %%%%%%%%%%%%%%%%%%%
y = y-mean(y);      %%%%%%%%%%%%%%%%%%%
% y = y(1:2:end);     %%%%%%%%%%%%%%%%%%%
y_freq = fftshift(fft(y));          %Frequency domain
df = os*BaudRate/length(y_freq);

if mod(length(y_freq),2) == 0   %when the length of sequence is an even number
    freq = [-length(y_freq)/2*df:df:-df, 0:df:df*(length(y_freq)/2-1)];
else                            %when the length of sequence is odd number
    freq = [-(length(y_freq)-1)/2*df:df:-df, 0:df:df*(length(y_freq)-1)/2];
end

if fig == 1
    figure
    plot(freq, abs(y_freq),'.-')
    title('Spectrum after square in time domain');
    % figure
    % plot(freq, angle(y_freq),'.-');
end



% BPF
a = find(y_freq == max(y_freq(find(freq>-BaudRate*1.002 & freq<-BaudRate*0.998))));   %-BaudRate근처 peak freq. 위치
b = find(y_freq == max(y_freq(find(freq>BaudRate*0.998 & freq<BaudRate*1.002))));     %+BaudRate근처 peak freq. 위치
f_clk_peak = freq(b);

% Gaussian fitting
num = 21;   %fitting length*0.5

temp = abs(y_freq(b-num:b+num));
f=fit([-num:num].',temp.','gauss2');

if fig == 1
    figure
    subplot(2,1,1)
    plot(f_clk_peak+df*[-num:num],temp,'.-')
    subplot(2,1,2)
    plot(f,[-num:num],temp)
end

temp = (f.a1*f.b1 + f.a2*f.b2)/(f.a1+f.a2); 
f_clk=f_clk_peak+temp*df;
% 
% fprintf('# Clock recovery\n');
% fprintf('Recovered clock frequency = %d GHz\n\n',double(f_clk)/1e9);

