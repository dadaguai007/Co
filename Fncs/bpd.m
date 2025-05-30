function  i = bpd(E1,E2,param)
%BPD
if nargin < 2
    param = struct();
end
R=1;
if isfield(param, 'R')
    R = param.R;
end
if ~isequal(size(E1), size(E2))
    error('E1 and E2 need to have the same shape');
end

i1 = pd(E1,param);
i2 = pd(E2,param);

i= R*(i1-i2);

end