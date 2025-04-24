classdef SymbolSynchronizer < handle
    %     loop_gain = 0.3
    % 只能实现 采样相位的偏移，还需要对采样时钟（频率）进行相应的纠正
    % 流程：接收信号 → 载波同步（频率补偿） → 符号同步（相位补偿） → 解调
    properties
        sync_type     % 同步算法类型: 'gardner', 'MM', 'EL'
        loop_gain     % 环路增益 控制环路收敛速度和稳定性的增益系数
    end

    methods
        function obj = SymbolSynchronizer(sync_type, loop_gain)
            % 构造函数
            obj.sync_type = sync_type;
            obj.loop_gain = loop_gain;
        end

        function [output_x, output_list] = process(obj, samples, samples_per_symbol)
            % 输入参数:
            %   samples: 复数采样数组 (1xN)
            %   samples_per_symbol: 每符号采样数 (sps)
            % 输出:
            %   output_x: 同步点标记向量 （用于调试）
            %   output_list: 同步后的符号序列

            in_index = 1;      % MATLAB索引从1开始
            sps = samples_per_symbol;
            tau = 0;           % 相位误差初始化

            output_x = zeros(size(samples));
            output_list = [];

            % 假设sps 为2 ， x12 和 x_pre 为相同的采样点
            % 假设sps 为4， x1代表第一个（index）， x12代表 index+2， x_pre 代表 index+3


            % 采取两个 符号 中的所有 点数 进行始终恢复  即： x1，x1,1，x2，x2,1，x3 （ 要到达第三个符号的第一个采样点）
            while (in_index + 2*sps - 1) <= length(samples) % 防止越界
                % 提取关键采样点 (注意MATLAB索引从1开始)
                x_1    = samples(in_index + 0*sps + 0*floor(sps/2));  %符号1的起始点，可用于MM算法
                x_12   = samples(in_index + 0*sps + 1*floor(sps/2));  %符号1与符号2之间的过渡点，用于Gardner算法
                x_pre2 = samples(in_index + 1*sps + 0*floor(sps/2) - 1); % 符号2起始点 前的一个采样点，可用于 早迟门算法

                % 用于所有算法：【x2】
                % Gardner算法：作为当前符号的参考点。
                % MM算法：与前一个符号比较。
                % EL算法：作为当前符号的“准时”门。

                x_2    = samples(in_index + 1*sps + 0*floor(sps/2));
                x_post2= samples(in_index + 1*sps + 0*floor(sps/2) + 1); % 符号2起始点 后一个采样点 ，用于 早迟门算法
                x_23   = samples(in_index + 1*sps + 1*floor(sps/2)); %  符号2与符号3之间的过渡点， 用于Gardner算法

                x_3    = samples(in_index + 2*sps + 0*floor(sps/2)); % 没有用处

                % 记录当前最佳采样点（因为按照顺序处理，最开始的信号会排序到末尾）
                output_list = [x_2; output_list]; % 头部插入
                output_x(in_index + 1*sps + 0*floor(sps/2)) = 1; % x2 相应的位置插入 参考点

                % 计算定时误差 （不同定时算法）
                switch obj.sync_type
                    case 'gardner'
                        % 使用 x_12（符号1中间点）和 x_23（符号2中间点）与 x_2（符号2起始点）计算误差
                        % Gardner算法
                        timing_error = real(x_2) * (real(x_12) - real(x_23));
                    case 'MM'
                        % Mueller & Müller算法
                        %  比较当前符号（x_2）与前一个符号（x_1）的符号方向
                        sng_x1 = sign(real(x_1));
                        sng_x2 = sign(real(x_2));
                        timing_error = real(x_2)*sng_x1 - real(x_1)*sng_x2;
                    case 'EL'
                        % 早迟门算法
                        sng_x2 = sign(real(x_2)); % 假设使用实部符号
                        timing_error = -sps * sng_x2 * (real(x_post2) - real(x_pre2))/2;
                    otherwise
                        error('Unknown synchronization type');
                end

                % 新相位计算 （pre相位减去 误差与环路增益的乘积）
                new_tau = tau - timing_error * obj.loop_gain;
                % 相位整数部分的插值（补偿量）
                index_step = floor(new_tau) - floor(tau);
                in_index = in_index + floor(sps + index_step); % 基础步进（一个符号周期）加上相位补偿量
                % 更新相位
                tau = new_tau ; % 保留小数部分
            end

            % 调整输出顺序
            output_list = flipud(output_list); % 翻转恢复正确顺序
        end
    end
end
