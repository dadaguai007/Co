function lim = convergenceCondition(Ex_fd, Ey_fd, Ex_conv, Ey_conv)
    % Verify the convergence condition for the trapezoidal integration.
    % Calculate the squared root of the MSE normalized by the power of the fields.
    lim = sqrt(norm(Ex_fd - Ex_conv)^2 + norm(Ey_fd - Ey_conv)^2) / sqrt(norm(Ex_conv).^2 + norm(Ey_conv).^2);
end


%这段代码是一个名为`convergenceCondition`的MATLAB函数，用于验证梯形积分的收敛条件。函数接收四个输入参数：`Ex_fd`和`Ey_fd`分别是x方向和y方向电场在细网格上的数值解，`Ex_conv`和`Ey_conv`分别是x方向和y方向电场在粗网格上的数值解。函数的目标是比较细网格和粗网格解之间的差异，以确定积分是否收敛。
% 函数的工作原理如下：
% 1. `norm(Ex_fd - Ex_conv)^2`计算细网格和粗网格在x方向电场解之间的均方误差（Mean Squared Error, MSE）的平方。
% 2. `norm(Ey_fd - Ey_conv)^2`计算细网格和粗网格在y方向电场解之间的均方误差的平方。
% 3. 这两个误差的平方相加，得到总的均方误差的平方。
% 4. `sqrt(norm(Ex_fd - Ex_conv)^2 + norm(Ey_fd - Ey_conv)^2)`计算总的均方误差的平方根，这是细网格和粗网格解之间差异的一种度量。
% 5. `norm(Ex_conv).^2`计算粗网格x方向电场解的功率。
% 6. `norm(Ey_conv).^2`计算粗网格y方向电场解的功率。
% 7. 这两个功率相加，然后取平方根，得到粗网格解的归一化功率。
% 8. 最后，将步骤4中得到的总的均方误差的平方根除以步骤7中得到的粗网格解的归一化功率，得到收敛条件的值`lim`。
% 简而言之，`lim`表示细网格和粗网格解之间的相对误差，用于判断梯形积分是否收敛。如果`lim`的值很小，则可以认为积分是收敛的；如果`lim`的值较大，则可能需要进一步细化网格或检查其他数值方法以提高精度。
