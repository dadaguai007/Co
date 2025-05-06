function lambda = entropy2lambda(H,C)
%entropy2lambda     Convert entropy values to required shaping parameter,
%                   lambda using a pre-stored look-up table (LUT)
...使用预存储的查找表 (LUT) 将熵值转换为所需的整形参数 lambda
%   This function performs the conversion of query entropy points, H, into
%   the correspondent shaping parameter, lambda, required by probabilistic
%   shaping. The probability distribution function is assumed to be 
%   exp(-lambda*abs(C).^2).
...此函数执行将查询熵点 H 转换为概率整形所需的相应整形参数 lambda。 假设概率分布函数为 exp(-lambda*abs(C).^2)。
%   INPUTS:
%       H := entropy values [1 x nEntropy]
%       M := size of the QAM constellation [scalar]  %scalar：标量
%
%   OUTPUTS:
%       lambda := shaping parameter [1 x nEntropy]  整形参数
%
%
%   Examples:
%       lambda = entropy2lambda(3.5:0.1:4,64);
%
%
% Authors: Fernando Guiomar
% Last Update: 23/10/2017

%% Optimize lambda to Achieve Desired Entropy  优化 lambda 以实现所需的熵
options = optimset('MaxFunEvals',1e4,'TolX',1e-4,'TolFun',1e-4,...
    'Display','none','PlotFcns',[]);  %optimset：优化集  TolX、TolFun：修改精度限制 皆为优化函数中设置参数的函数
[lambda,err] = fminsearch(@(lambda) fitMaxwellBoltzman(lambda,C,H),...
    0,options);             %函数命令拟合
...                          最常用的函数拟合命令为fit，语法为
...                          [拟合结果 拟合精度]＝fit（X数据，Y数据，‘拟合类型’）

end

%% Maxwell-Boltzman Fitting Function  Maxwell-Boltzmann 拟合函数
function err = fitMaxwellBoltzman(lambda,C,H)
    symProb = exp(-lambda*abs(C).^2);
    symProb = symProb/sum(symProb);
    entropy = -sum(symProb.*log2(symProb));%MB分布的Px用symProb描述    
    err = abs(H-entropy);  %err为拟合精度
end
