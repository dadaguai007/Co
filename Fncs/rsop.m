function [rx,ry]=rsop(inx,iny)

t_spacing = 5e-11;
T = length(inx);
t = 0:t_spacing:T*t_spacing;
N =12;
w_alpha = 5e3; %krad/s---rad/s
w_phi = 13e3;
w_kappa =11e3;

U = cell(1,T);

kappa0 = pi+randn(N,1)*2*pi;
alpha0 = pi+randn(N,1)*2*pi;
phi0 = pi+randn(N,1)*2*pi;

for a = 1:T
    alpha = w_alpha.*t(a)+alpha0;
    phi = w_phi.*t(a)+phi0;
    kappa = w_kappa.*t(a)+kappa0;
    oo11 = cos(kappa).*exp(1j*(alpha));
    oo12 = -sin(kappa).*exp(-1j*(phi));
    oo21 = -conj(oo12);
    oo22 = conj(oo11);
    y=1;
    for b = 1:N
        y = y*[oo11(b) oo12(b);oo21(b) oo22(b)];
    end
    U{a} = y;
rx(a) = U{a}(1,1).*inx(a)+U{a}(1,2).*iny(a);
ry(a) = U{a}(2,1).*inx(a)+U{a}(2,2).*iny(a);
end
