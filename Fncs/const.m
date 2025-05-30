% class const, filename 'const.m'
%   This class stores four constants used by the framework.
% 
% Examples:
%   @code
%   lambda = 1550e-9;
%   fc = const.c / lambda;
%   @endcode
%
% from @OCG-LAB OPC Framework
%
% 简体中文版注释：
% 类 const，文件名为 'const.m'
%   该类存储了框架所用到的四个常量。
% 
% 示例：
%   @code
%   lambda = 1550e-9;
%   fc = const.c / lambda;
%   @endcode
%
% 来自 @OCG-LAB OPC 框架

classdef const
    properties (Constant)
        % Boltzmann constant (J/K)
        kB = 1.3806488e-23;  
        % Elementary charge (C)
        q = 1.602176565e-19; 
        % Planck constant (J*s)
        h = 6.62606957e-34;  
        % Speed of light in vacuum (m/s)
        c = 299792458;       
    end
end