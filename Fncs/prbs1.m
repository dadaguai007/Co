function prbs_seq = prbs1(n,N,h)
% h=1,为正常；h=0,为全取
if nargin < 2
    N = 2^n - 1;
end
switch n
    % 常用的特征多项式
    case 7
        g = [7 6];
    case 9
        g = [9 5];
    case 11
        g = [11 9];
    case 15
        g = [15 14];
    case 23
        g  =[23 18];
    case 31
        g = [31 28];
    otherwise
        error('Allowed lengths: 2^{7|9|13|15|23|31}-1');
end

% PRBS Generation
z = zeros(1,2^n-1);
%初始seed 100……
z(1)=1;
for i=(n+1):(2^n-1)
    q=z(i-g(1));
    for j=2:length(g)
        q=xor(q,z(i-g(j)));
    end
    z(i) = q;
end
z = [z 0];
z=z(:);
if h>0
    x=logical(z(1:2^n-1));
else
    x=logical(z);
end
% repeat to have N bits
nRepeat = ceil(N/length(x));
x_tmp = repmat(x,nRepeat,1);
prbs_seq = x_tmp(1:N);
end