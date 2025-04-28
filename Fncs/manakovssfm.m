function [Ech, param] = manakovssfm(Ei, param)
% 双偏的SSFM
% Default parameters
Ltotal = 400; %km
Lspan = 80;
hz= 0.5;
alpha=0.2;
D = 16;
gamma = 1.3;
Fc = 193.1e12;
NF = 4.5;
amp='edfa';
maxIter = 10;
tol=1e-5;
nlprMethod='true';
maxNlinPhaseRot=2e-2;

% check input parameters
if isfield(param, 'Fs')
    Fs = param.Fs;
end
if isfield(param, 'Ltotal')
    Ltotal = param.Ltotal;
end
if isfield(param, 'Lspan')
    Lspan = param.Lspan;
end
if isfield(param, 'hz')
    hz = param.hz;
end
if isfield(param, 'alpha')
    alpha = param.alpha;
end
if isfield(param, 'D')
    D = param.D;
end
if isfield(param, 'gamma')
    gamma = param.gamma;
end
if isfield(param, 'Fc')
    Fc = param.Fc;
end
if isfield(param, 'NF')
    NF = param.NF;
end
if isfield(param, 'amp')
    amp = param.amp;
end
if isfield(param, 'maxIter')
    %积分中的最大迭代次数
    maxIter = param.maxIter;
end

if isfield(param, 'tol')
    %积分中的收敛公差
    tol = param.tol;
end
if isfield(param, 'nlprMethod')
    % 自适应步长基于非线性相位旋转的开关
    nlprMethod = param.nlprMethod;
end
if isfield(param, 'maxNlinPhaseRot')
    %最大非线性相位旋转容差的阈值
    maxNlinPhaseRot = param.maxNlinPhaseRot;
end
%输出的光纤跨度的索引列表
saveSpanN=Ltotal/Lspan;



% channel parameters
c = 299792458; % speed of light (vacuum) in m/s
c_kms = c / 1e3;
lambda = c_kms / Fc;
Alpha = alpha / (10 * log(exp(1)));
beta2 = -(D * lambda^2) / (2 * pi * c_kms);


% edfa parameters
paramAmp = struct();
paramAmp.G = alpha * Lspan;
paramAmp.NF = NF;
paramAmp.Fc = Fc;
paramAmp.Fs = Fs;
paramAmp.type = 'noise';
% generate frequency axis
Nfft = length(Ei);
omega = 2 * pi * Fs * (-0.5:1/Nfft:0.5-1/Nfft);
omega=ifftshift(omega);

Nspans = floor(Ltotal / Lspan);

% Polarization vectors
Ech_x = Ei(:, 1)';
Ech_y = Ei(:, 2)';

% define static part of the linear operator
argLimOp = (-(Alpha / 2) - 1i * (beta2 / 2) * (omega.^2));

% row Ech_x 
if size(Ech_x, 1) > 1
    % change to martix size(Ech_x.row, argLimOp.column)
    argLimOp = repmat(argLimOp, [size(Ech_x, 1), 1]);
    warning('Please cheak the input');
else
    argLimOp = reshape(argLimOp, [1, length(argLimOp)]);
end

% Stores optical signals within a specified span
if ~isempty(saveSpanN)
    % (Ei.row, Ei.column * len(saveSpanN))的全零数组.
    % 有多少段 ， 就有多少列进行存储 
    % saveSpaN 为一个数
    % 将变量indRecSpan初始化为0，用于后续记录当前存储的位置。
    Ech_spans = zeros([size(Ei, 1), size(Ei, 2) * saveSpanN]);
    indRecSpan = 1;
end

for spanN = 1:Nspans
    Ex_conv = Ech_x;
    Ey_conv = Ech_y;

    z_current = 0;

    % fiber propagation steps
    while z_current < Lspan
        Pch = Ech_x .* conj(Ech_x) + Ech_y .* conj(Ech_y);

        phiRot = nlinPhaseRot(Ex_conv, Ey_conv, Pch, gamma);

        if strcmp(nlprMethod, 'true')
            hz_ = maxNlinPhaseRot / max(phiRot);
            if Lspan - z_current < maxNlinPhaseRot / max(phiRot)
                hz_ = Lspan - z_current;
            end
        elseif Lspan - z_current < hz
            hz_ = Lspan - z_current;
        else
            hz_ = hz;
        end

        % define the linear operator
        linOperator = exp(argLimOp * (hz_ / 2));

        % First linear step (frequency domain)
        Ex_hd = ifft(fft(Ech_x) .* linOperator);
        Ey_hd = ifft(fft(Ech_y) .* linOperator);

        % Nonlinear step (time domain)
        for nIter = 1:maxIter
            % take note , the nonlinear should be test
            rotOperator = exp(-1i * phiRot * hz_);

            Ech_x_fd = Ex_hd .* rotOperator;
            Ech_y_fd = Ey_hd .* rotOperator;

            % Second linear step (frequency domain)
            Ech_x_fd = ifft(fft(Ech_x_fd) .* linOperator);
            Ech_y_fd = ifft(fft(Ech_y_fd) .* linOperator);

            % check convergence o trapezoidal integration in phiRot
            lim = convergenceCondition(Ech_x_fd, Ech_y_fd, Ex_conv, Ey_conv);

            Ex_conv = Ech_x_fd;
            Ey_conv = Ech_y_fd;

            if lim < tol
                break;
            elseif nIter == maxIter
                warning(['Warning: SSFM error tolerance was not achieved in ' num2str(maxIter) ' iterations']);
            end

            phiRot = nlinPhaseRot(Ex_conv, Ey_conv, Pch, gamma);
        end

        Ech_x = Ech_x_fd;
        Ech_y = Ech_y_fd;

        z_current = z_current + hz_;  % update propagated distance
    end

    % amplification step
    if strcmp(amp, 'edfa')
        Ech_x = edfa(Ech_x, paramAmp);
        Ech_y = edfa(Ech_y, paramAmp);
        fprintf('the edfa %.2f span \n',spanN);
    elseif strcmp(amp, 'ideal')
        Ech_x = Ech_x .* exp(Alpha/2 * Lspan);
        Ech_y = Ech_y .* exp(Alpha/2 * Lspan);
        fprintf('the ideal %.2f span \n',spanN)
    elseif strcmp(amp, 'none')
        Ech_x = Ech_x .* exp(0);
        Ech_y = Ech_y .* exp(0);
    end

    if ~isempty(saveSpanN) 
        Ech_spans(:, indRecSpan:(indRecSpan + 1)) = [Ech_x(:) Ech_y(:)];
        indRecSpan = indRecSpan + 2;
    end
end

if ~isempty(saveSpanN)
    Ech = Ech_spans(:,indRecSpan-2:(indRecSpan-2+1));
else
    Ech = Ei;
    Ech(:, 1) = Ech_x.';
    Ech(:, 2) = Ech_y.';
end


end
