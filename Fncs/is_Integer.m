function out = is_Integer(v)
%% Check if v is integer with tolerance tol
%判断是否为整数
tol = 1e-3;
err = abs(v - round(v));
out = (err < tol);