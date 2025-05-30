function [Eout, phase_est] = viterbiviterbi(E, N, M)
% Viterbi-Viterbi blind phase recovery for  M-PSK signal
% N是指要对信号进行平均的样本数，而不是信号的总样本数。
% 通过计算E信号的N个相邻样本的平均值来实现平均操作，以减小相位恢复的噪声影响，从而提高相位恢复的准确性。
% 确保输入向量是row向量,或者多维行向量，多维向量当然是2*N形式

%针对 QPSK 进行相偏估计
sizeVector = size(E);

if length(sizeVector) > 1 && sizeVector(2) >= 1
   disp('This is a mutal-dimensional vector')
   E2d=E;
elseif length(sizeVector) > 1
    if isrow(E)
        E2d=E;
    else
        E2d=E.';
    end
end
Eout = zeros(size(E2d));
    
    for i = 1:size(E2d, 1)
        Et = E2d(i, :);
        L = length(Et);
        phi = angle(Et);
        E_raised = exp(1i * phi).^M;


        %E_raised：将信号提升到M次幂。
        % sa：将提升后的信号划分为长度为N、步长为N-1的重叠段。
        % phase_est：通过对每个段求和来估计相位。
        % 估计的相位然后进行展开，并应用了校正。
        sa = buffer(E_raised, N, N-1, 'nodelay');
        phase_est = sum(sa, 1);
        phase_est = unwrap(angle(phase_est));
        phase_est = (phase_est - pi) / M;
        
        % 奇数范围和偶数范围
        if mod(N, 2)
            Eout(i, (N - 1) / 2 + 1:L - (N - 1) / 2) = ...
                E2d(i, (N - 1) / 2 + 1:L - (N - 1) / 2) .* exp(-1i * phase_est);
        else
            Eout(i, N / 2:L - (N / 2 - 1)) = ...
                E2d(i, N / 2:L - (N / 2 - 1)) .* exp(-1i * phase_est);
        end
    end
    
    % Handle QPSK case (uncomment if needed)
%     if M == 4
%         % QPSK needs pi/4 shift
%         phase_est = phase_est + pi/4;
%     end
%     
%     if isrow(E)
%         Eout = Eout(:).';
%         phase_est = phase_est(:).';
%     end
end
% 该函数实现了Viterbi-Viterbi肓相位恢复算法，用于对M-PSK信号进行相位恢复。该算法通过对输入的电场信号E进行处理，恢复出信号中的相位信息，并返回恢复后的电场信号Eout和相位估计结果phase_est。该函数的输入参数包括电场信号E、样本数量N和M-PSK的阶数M。函数首先将信号E转换为2维数组E2d，并创建一个与E2d相同大小的全零数组Eout用于存储恢复后的信号。然后，对于E2d中的每一行信号Et，算法执行以下操作：
% 
% 计算信号Et的相位phi。
% 将信号Et的相位phi提升到M次方，并得到E_raised。
% 将E_raised按长度为N的子序列进行划分，得到子序列数组sa。
% 对sa按行求和，并得到phase_est。
% 对phase_est进行unwrap操作，消除相位的折叠现象。
% 将phase_est除以M，并减去pi，得到归一化的相位值。
% 根据N的奇偶性，将E2d中的相应部分与np.exp(-1.j*phase_est)相乘，然后将结果赋给Eout。
% 最后，根据输入信号E的维度，将恢复后的电场信号Eout和相位估计结果phase_est返回。如果输入信号E是一维的，则返回的Eout和phase_est也是扁平化的一维数组；如果输入信号E是二维的，则返回的Eout和phase_est是二维数组。