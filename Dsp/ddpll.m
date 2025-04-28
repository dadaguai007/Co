function theta = ddpll(Ei, Ts, Kv, tau1, tau2, constSymb, symbTx, pilotInd)
    % Decision-directed Phase-locked Loop (DDPLL) algorithm
%parameter
% Ei : 接收到的星座符号。
% Ts :  符号周期。
% Kv : 环路滤波器增益。
% tau1 : 环路滤波器参数1
% tau2 : 环路滤波器参数2。
% constSymb : 复数数组 理想的星座符号。
% symbTx : 复数数组 发送的符号序列。
% pilotInd : 整数数组 导频符号的位置索引。
% θ : 实数数组 时变的相位偏移估计值。
    % Get the number of symbols and modes from the size of Ei
    % Ei should be the N × （1,2）
    [M,N]=size(Ei);
    if M < N
       error('the Ei should be the N × models')
    end
    
    [M,N]=size(symbTx);
    if M < N
       error('the symbTx should be the N × models')
    end
    [nSymbols, nModes] = size(Ei);

    % Initialize the theta array
    theta = zeros(nSymbols, nModes);

    % Loop filter coefficients 一个三阶的一阶全通滤波器，可以实现相位跟踪和相位噪声抑制
    a1b = [
        1;
        Ts / (2 * tau1) * (1 - 1 / tan(Ts / (2 * tau2)));
        Ts / (2 * tau1) * (1 + 1 / tan(Ts / (2 * tau2)))
    ];
% 存储相位检测器和环路滤波器的输出，u[1]为环路滤波器的输出，u[2]和u[3]为相位检测器的前后两个输出，初始值为零
    u = zeros(1, 3);  % [u_f, u_d1, u_d]

    for n = 1:nModes
        u(3) = 0;  % Output of phase detector (residual phase error)
        u(1) = 0;  % Output of loop filter

        for k = 1:nSymbols
            u(2) = u(3);
            % Remove estimate of phase error from input symbol
            %得到去除相位偏移的符号Eo
            Eo = Ei(k, n) * exp(1i * theta(k, n));

            % Slicer (perform hard decision on symbol)
            %即根据导频符号或最近的星座符号，计算相位误差信号u[3]，它是Eo和发送的符号或星座符号的虚部
            if ismember(k, pilotInd)
                % phase estimation with pilot symbol
                % Generate phase error signal (also called x_n (Meyer))
                u(3) = imag(Eo * conj(symbTx(k, n)));
            else
                % find the closest constellation symbol
                [~, decided] = min(abs(Eo - constSymb));
                % Generate phase error signal (also called x_n (Meyer))
                u(3) = imag(Eo * conj(constSymb(decided)));
            end

            % Pass phase error signal through Loop Filter  get the filter
            % output
            u(1) = sum(a1b .* u');

            % Estimate the phase error for the next symbol
            %当前相位偏移θ[k, n]减去环路滤波器增益Kv和环路滤波器输出u[1]的乘积
            if k < nSymbols
                theta(k + 1, n) = theta(k, n) - Kv * u(1);
            end
        end
    end
end
