function [Eo, theta] = cpr(Ei, symbTx, paramCPR)
% Carrier phase recovery function (CPR)


if nargin < 2
    symbTx = [];
end

if nargin < 3
    paramCPR = [];
end

% M表示星座的阶数，默认为4。
% constType表示星座的类型，可选为'qam'或'psk'，默认为'qam'。
% B表示BPS算法中的测试相位的个数，默认为64。
% N表示BPS算法中的移动平均窗口的长度，默认为35。
% Kv，表示DDPLL算法中的环路滤波器的增益，默认为0.1。
% tau1表示DDPLL算法中的环路滤波器的参数1，默认为1 / (2 * np.pi * 10e6)。
% tau2表示DDPLL算法中的环路滤波器的参数2，默认为1 / (2 * np.pi * 10e6)。
% Ts表示符号周期，默认为1 / 32e9。
% pilotInd表示导频符号的位置的索引，默认为np.array([len(Ei) + 1])，表示没有导频符号。

% check input parameters

if isfield(paramCPR, 'alg')
    alg = paramCPR.alg;
end

if isfield(paramCPR, 'B')
    B = paramCPR.B;
end
if isfield(paramCPR, 'N')
    N = paramCPR.N;
end
if isfield(paramCPR, 'Kv')
    Kv = paramCPR.Kv;
end
if isfield(paramCPR, 'tau1')
    tau1 = paramCPR.tau1;
end
if isfield(paramCPR, 'tau2')
    tau2 = paramCPR.tau2;
end
if isfield(paramCPR, 'Ts')
    Ts = paramCPR.Ts;
end
if isfield(paramCPR, 'pilotInd')
    pilotInd = paramCPR.pilotInd;
end
if isfield(paramCPR, 'constSymb')
    constSymb = paramCPR.constSymb;
end
% alg = paramCPR.alg;
% % M = paramCPR.M;
% % constType = paramCPR.constType;
% B = paramCPR.B;
% N = paramCPR.N;
% Kv = paramCPR.Kv;
% tau1 = paramCPR.tau1;
% tau2 = paramCPR.tau2;
% Ts = paramCPR.Ts;
% pilotInd = paramCPR.pilotInd;
% constSymb=paramCPR.constSymb;


if isempty(Ei)
    error('Input signal Ei cannot be empty.');
end

% Reshape Ei if it has only one column
if size(Ei, 1) == 1
    Ei = reshape(Ei, length(Ei), 1);
end



% 4th power frequency offset estimation/compensation
[Ei, ~] = fourthPowerFOE(Ei, 1/Ts);
%norm
Ei = Ei./sqrt(mean(abs(Ei).^2));

if strcmp(alg, 'ddpll')
    theta = ddpll(Ei, Ts, Kv, tau1, tau2, constSymb, symbTx, pilotInd);
elseif strcmp(alg, 'bps')
    theta = bps(Ei, floor(N / 2), constSymb, B);
elseif strcmp(alg,'viterbi')
    theta = viterbi(Ei, N);
else
    error('CPR algorithm incorrectly specified.');
end

% 相位估计值的数组θ进行展开和缩放处理，使其在[-pi, pi]的范围内。
theta = unwrap(4 * theta, 1) / 4;

%compensation
Eo = Ei .* exp(1i * theta);
%如果输出信号的模式数为1，转换为一维
if size(Eo, 2) == 1
    Eo = Eo(:);
    theta = theta(:);
end
end
