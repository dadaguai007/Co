function [yEq, H, errSq, Hiter] = mimoAdaptEqualizer(x, dx, param)
% N-by-N MIMO adaptive equalizer.

% Algorithms available: 'cma', 'rde', 'nlms', 'dd-lms', 'da-rde', 'rls', 'dd-rls', 'static'.

% check input parameters
if nargin < 2
    dx = [];
end
if nargin < 3
    param = struct();
end

% Pre-convergence time 预先均衡的次数
numIter =1;

% Eq parameters
nTaps=15;
mu=1e-3;
%RLS forgetting factor
lambdaRLS=0.99;
SpS=2;
% coefficient taps
H=[];
% Eq length
L=[];
Hiter=[];
storeCoeff='true';
% alg is the type
alg='nlms';

%constSymb should be get from the   external，don't forget 

if isfield(param, 'numIter')
    numIter = param.numIter;
end
if isfield(param, 'nTaps')
    nTaps = param.nTaps;
end
if isfield(param, 'mu')
    mu = param.mu;
end
if isfield(param, 'lambdaRLS')
    lambdaRLS = param.lambdaRLS;
end
if isfield(param, 'SpS')
    SpS = param.SpS;
end
if isfield(param, 'H')
    H = param.H;
end
if isfield(param, 'L')
    L = param.L;
end
if isfield(param, 'Hiter')
    Hiter = param.Hiter;
end
if isfield(param, 'storeCoeff')
    storeCoeff = param.storeCoeff;
end
if isfield(param, 'alg')
    alg = param.alg;
end
if isfield(param, 'constSymb')
    constSymb = param.constSymb;
end


%  the all signal sequences to be disposed in columns,
% it is conventence for the modes
if isempty(dx)
    dx = x;
end

% data should be the L * Nmodel
if size(x, 2) > size(x, 1)
    x = x.';
end

% Reference Sequence should be the L * Nmodel
if size(dx, 2) > size(dx, 1)
    dx = dx.';
end

nModes = size(x, 2);
Lpad = floor(nTaps / 2);
zeroPad = zeros(Lpad, nModes);
x = [zeroPad; x; zeroPad];  % pad start and end of the signal with zeros


% 除去滤波器本身所占用的一组采样点，最后得到输出信号的符号数
%将 X 的每个元素朝零方向四舍五入为最近的整数。
%此操作实际上是通过删除 X 中每个数的小数部分，将它们截断为整数。
totalNumSymb = fix((length(x) - nTaps) / SpS + 1);

if isempty(L)  % if L is not defined
    % L 为输出的符号数
    L = totalNumSymb;  % Length of the output (1 sample/symbol) of the training section
end

% MIMO equize system
if isempty(H)  % if H is not defined
    H = zeros(nModes.^2, nTaps);
    
    for initH = 1:nModes  % initialize filters' taps
        H(initH + (initH-1 ) * nModes, floor(nTaps / 2) + 1) = 1;  % Central spike initialization
    end
    %         for initH = 1:nModes  % initialize filters' taps
    %             H(initH , floor(nTaps / 2) + 1) = 1;  % Central spike initialization
    %         end
end

% Equalizer training:
% alg 可以 写成{'cma'，'cma'}形式，相对应 L，mu
if iscell(alg)
    %以下过程包括预均衡 和 均衡
    % 预均衡就是取 前一部分 信号进行 均衡 ， 均衡之后 再对之后的代码进行均衡，能加快进度是真的
    yEq = zeros(totalNumSymb, nModes);
    errSq = zeros(nModes, totalNumSymb);
    nStart = 0;
    for indstage = 1:length(alg)
        runAlg = alg{indstage};
        nEnd = nStart + L(indstage);
        if indstage == 1
            for indIter = 1:numIter
                disp([runAlg, ' pre-convergence training iteration #', num2str(indIter)]);
                [yEq(nStart + 1:nEnd, :), H, errSq(:, (nStart + 1):nEnd), Hiter] = coreAdaptEq(...
                    x( (nStart*SpS+1) : (nEnd*SpS+2*Lpad), :),...
                    dx((nStart+1):nEnd, :),...
                    SpS, H, L(indstage), mu(indstage), lambdaRLS, nTaps,...
                    storeCoeff, runAlg, constSymb);

                disp([runAlg, ' MSE = ', num2str(mean(errSq(:, nStart + 1:nEnd)))]);
            end
        else
            % Formal use
            [yEq(nStart + 1:nEnd, :), H, errSq(:, nStart + 1:nEnd), Hiter] = coreAdaptEq(...
                x((nStart*SpS+1) : (nEnd*SpS+2*Lpad), :),...
                dx(nStart+1:nEnd, :),...
                SpS, H, L(indstage), mu(indstage), lambdaRLS, nTaps,...
                storeCoeff, runAlg, constSymb);

            disp([runAlg, ' MSE = ', num2str(mean(errSq(:, nStart + 1:nEnd)))]);
        end
        nStart = nEnd;
    end
else
    [yEq, H, errSq, Hiter] = coreAdaptEq(x(Lpad+1:end-Lpad), dx, SpS, H, L, mu,lambdaRLS, nTaps, storeCoeff, alg, constSymb);
    disp([alg, ' MSE = ', num2str(mean(errSq))]);
end
end
