function [out,e,debugInfo] = mycma(in,Ntaps,mu)
%% validate input
if size(in,2) ~= 2
    error('Input data should be N*2 matrix!');
end
if ~mod(Ntaps,2)
    Ntaps = Ntaps + 1;
    fprintf('Ntaps should be odd number, thus new Ntaps is %d\n',Ntaps);
end

%% initialize the filters
[hxx,hxy,hyx,hyy] = deal(zeros(Ntaps,1));
hxx(ceil(Ntaps/2)) = 1;
hyy(ceil(Ntaps/2)) = 1;
window = (1:Ntaps) - ceil(Ntaps/2);

%% equalization start here
for k = ceil(Ntaps/2):2:size(in,1)-ceil(Ntaps/2)+1
    % get the symbols
    xin = in(k + window, 1);
    yin = in(k + window, 2);
    % get the output
    xout = hxx'*xin + hxy'*yin; % Hermittian conjugate
    yout = hyx'*xin + hyy'*yin; % Hermittian conjugate
    % calculate the errors
    ex = 1 - abs(xout)^2;
    ey = 1 - abs(yout)^2;
    % update the equalizers
    hxx = hxx + mu*ex*conj(xout).*xin;
    hxy = hxy + mu*ex*conj(xout).*yin;
    hyx = hyx + mu*ey*conj(yout).*xin;
    hyy = hyy + mu*ey*conj(yout).*yin;
    % save the output
    out((k+1)/2,:) = [xout,yout]; % k is odd number
    e(k,:) = [ex,ey];
end

debugInfo.h = [hxx,hxy,hyx,hyy];
end
