classdef DspSyncDecoding < handle
    % 设计一个 专门 为 DSP 流程中  同步，解码 的类
    % 流程：接收信号 → 载波同步（频率补偿） → 符号同步（相位补偿） → 解调
    properties
        % 为了方便后续添加变量，采用结构体方式进行变量命名
        signalPHY; % 传入的信号参数
        Nr%无量纲参数
        Implementation;% 参考信号实施参数
    end

    methods
        function obj = DspSyncDecoding(varargin)
            %初始化类的属性
            if numel(varargin) == 7
                obj.signalPHY.fs            = varargin{1} ;% 接收信号的采样率
                obj.signalPHY.fb            = varargin{2} ;% 接收信号的波特率
                obj.signalPHY.M             = varargin{3} ;% 接收信号的格式
                obj.Nr.sps                  = varargin{4} ;% 上采样率
                obj.Nr.time_osr             = varargin{5} ;% 时钟信号的上采样率
                obj.Nr.ncut_index           = varargin{6} ;% 误码计算起始位置
                obj.Implementation.ref      = varargin{7} ;% 参考信号
            end
            obj.Implementation.sync_type = 'gardner';      % 时钟误差计算方式
            obj.Implementation.loop_gain = 0.03;           % PLL的环路增益, 控制环路收敛速度和稳定性的增益系数
            obj.Implementation.TED       ="MLTED" ;        % 时钟误差计算方式，另一种实现方式
            obj.Implementation.eta       = 1;              % 环路阻尼因子（稳定性控制）
            obj.Implementation.Bn_Ts     = 0.01;           % 环路带宽 × 符号周期（控制同步速度）
        end


        function clk_output=ClockRecovery(obj,input_signal)
            %数字时钟恢复（载波频率）
            f_up=obj.Nr.time_osr *obj.signalPHY.fb; % 时钟恢复所需要的采样率
            data_up = resample(input_signal,f_up,obj.signalPHY.fs);  %一般采用四倍上采样，满足后续时钟恢复频率
            % 使用时钟恢复函数
            f_clk_rec = cr(data_up.',obj.signalPHY.fb,obj.Nr.time_osr,0); % 注意：f_clk_rec是对fb的估计，不是对fb*osr的估计！
            % 创建时间轴
            [~, t_cur] = freq_time_set(length(data_up), f_up);
            %t_cur = 1/obj.Nr.time_osr/obj.signalPHY.fb*(0:length(data_up)-1); % 注意：插值前data1的采样率是osr*fb
            % 正确时间信号，重建时间轴
            t_new = 0:1/f_clk_rec/obj.Nr.time_osr:t_cur(end); % 注意：插值后data2的采样率是osr*f_clk_rec
            % 插值
            data1 = interp1(t_cur,data_up,t_new,'spline');
            % 转换到所需要的采样率
            clk_output = resample(data1,obj.Nr.sps ,obj.Nr.time_osr);
            % 归一化
            clk_output=pnorm(clk_output);
        end


        function [syncedSignal,sync_position] = Synchronization(obj, clk_output)
            % 执行同步操作
            % 参考序列 达到 需要的采样率
            ref_sps = repelem(obj.Implementation.ref,obj.Nr.sps);
            [syncedSignal,~,sync_position] = sync(clk_output,ref_sps);
        end


        function syncedSignal=Time_recovery(obj,input_signal)
            % 对输入信号执行时间（频率）恢复 和 同步
            clk_output=obj.ClockRecovery(input_signal);
            [syncedSignal,~] = obj.Synchronization(clk_output);
        end

        % 采样相位恢复, 并进行序列的同步 (实现最佳相位采样)
        function   [DeWaveform,P,OptSampPhase,MaxCorrIndex]=time_phase_Recovery(obj,input_signal)
            % 参考信号
            label=obj.Implementation.ref;
            % 对输入信号执行时间（频率）恢复 和 同步
            clk_output=obj.ClockRecovery(input_signal);
            % 数字采样相位恢复
            fs_up=obj.Nr.sps*obj.signalPHY.fb; % 信号所采取的时钟
            [DeWaveform,P,OptSampPhase,MaxCorrIndex] = Quick_Syn_Vec(clk_output,label,1/fs_up,1/obj.signalPHY.fb);

        end



        % 采样相位恢复,（实现最佳相位采样）
        function [output_list,output_x] = recoverOptimalSamplingPoints(obj,samples)
            % 实现恢复最佳采样点的算法
            % 输入参数:
            %   samples: 复数采样数组 (1xN)
            %   samples_per_symbol: 每符号采样数 (sps)
            % 输出:
            %   output_x: 同步点标记向量 （用于调试）
            %   output_list: 同步后的符号序列

            in_index = 1;      % MATLAB索引从1开始
            sps = obj.Nr.sps;
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
                switch obj.Implementation.sync_type
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
                new_tau = tau - timing_error * obj.Implementation.loop_gain;
                % 相位整数部分的插值（补偿量）
                index_step = floor(new_tau) - floor(tau);
                in_index = in_index + floor(sps + index_step); % 基础步进（一个符号周期）加上相位补偿量
                % 更新相位
                tau = new_tau ;
            end

            % 调整输出顺序
            output_list = flipud(output_list); % 翻转恢复正确顺序

        end


        % 最佳相位恢复（使用解析法创造 TED）
        function [rxSync,Kp] = TED_recoverOptimalSamplingPoints(obj,rxSeq,mfOut,rollOff,rcDelay,const)
            % 计算环路增益(计算定时误差检测器（TED）的增益 )
            Kp  = calcTedKp(obj.ImplementationTED, rollOff);
            % 两种方法计算增益：'analytic'（解析法）或 'simulated'（仿真法），默认为解析法

            % 计数器增益（类似VCO增益）
            K0 = -1;                            % 控制环路响应速度
            % PI : PI环路的控制增益
            [ K1, K2 ] = piLoopConstants(Kp, K0, obj.Implementation.eta, obj.Implementation.Bn_Ts, obj.Nr.sps);

            fprintf("Loop constants:\n");
            fprintf("K1 = %g; K2 = %g; Kp = %g\n", K1, K2, Kp);

            % mfOut为匹配滤波后的信号，rxSeq为接收信号

            % debug 监视器
            debug_tl_static  = 0; % Show static debug plots after sync processing
            debug_tl_runtime = 0; % Open scope for debugging of sync loop iterations

            % 对参考信号进行归一化
            Ksym = modnorm(const, 'avpow', 1);
            const = Ksym * const;

            %symbol timing recovery implementation
            [ rxSync ] = symbolTimingSync(obj.Implementation.TED, obj.intpl, obj.Nr.sps, rxSeq, mfOut, K1, K2, ...
                const, Ksym, rollOff, rcDelay, debug_tl_static, debug_tl_runtime);

        end





        % 系统自带的最佳相位采样
        function rxSync = systemRecoverOptimalSamplingPoints(obj,mfOut,rollOff)

            % 计算环路增益(计算定时误差检测器（TED）的增益 )
            Kp  = calcTedKp(obj.ImplementationTED, rollOff);
            % 使用函数将缩写与全称 相对应
            tedMap = containers.Map({'ELTED', 'ZCTED', 'GTED', 'MMTED'}, ...
                {'Early-Late (non-data-aided)', ...
                'Zero-Crossing (decision-directed)', ...
                'Gardner (non-data-aided)', ...
                'Mueller-Muller (decision-directed)'
                });
            if strcmp(obj.Implementation.TED, "MLTED")
                warning("MLTED not supported by MATLAB's synchronizer - using ZCTED");
                matlabTed = "ZCTED";
            else
                matlabTed = obj.Implementation.TED;
            end

            SYMSYNC = comm.SymbolSynchronizer(...
                'TimingErrorDetector', tedMap(matlabTed), ...
                'SamplesPerSymbol', obj.Nr.sps, ...
                'NormalizedLoopBandwidth', obj.Implementation.Bn_Ts, ...
                'DampingFactor', obj.Implementation.eta, ...
                'DetectorGain', Kp);

            % mfOut为 匹配滤波输出 或者 为 接收信号
            rxSync = step(SYMSYNC, mfOut);


        end



        % 系统自带的载波偏移、相位补偿
        function  [recovered_carrier,CFO_error,mean_phase_error]=systemCarrierRecover(obj,input_signal)

            % 载波提取（四次方+滤波）
            carrier = input_signal.^4;              % 四次方运算产生4倍频分量
            carrier = carrier - mean(carrier);    % 去除直流偏移

            % 设计Chebyshev I型带通滤波器
            f_p4 = 4*obj.signalPHY.fb; % 预期载波频率
            [b, a] = cheby1(4, 0.5, [f_p4*0.98, f_p4*1.02]/(fs/2), 'bandpass');
            carrier_iso = filter(b, a, carrier);  % 滤波提取载波

            % 限幅处理（转换为方波）
            limited_carrier = sign(carrier_iso);  % 二值化限幅
            limited_carrier = limited_carrier(:); % 确保列向量

            % PLL载波恢复
            pll = comm.CarrierSynchronizer(...    % 创建锁相环对象
                'Modulation', 'QAM',...           % 设置调制类型
                'SamplesPerSymbol', 1,...         % 每符号采样数
                'DampingFactor', 0.707,...        % 阻尼系数（临界阻尼）
                'NormalizedLoopBandwidth', 0.01); % 归一化环路带宽

            % 执行载波恢复,输出经过补偿的信号
            [recovered_carrier, ~] = pll(limited_carrier);

            % 性能分析
            % 载波频率偏移计算
            estimated_freq = obj.signalPHY.fb; % 理论载波频率
            recovered_freq = abs(mean(diff(angle(recovered_carrier))))*obj.signalPHY.fs/(2*pi);
            %  载波偏移
            CFO_error = abs(estimated_freq - recovered_freq);
            fprintf(' CFO is : %.2f Hz\n', CFO_error);

            % 相位误差分析
            ideal_carrier = cos(2*pi*estimated_freq*t); % 理想参考载波
            % 相位误差均值
            mean_phase_error = mean(abs(angle(ideal_carrier) - angle(recovered_carrier)));

        end


        function [decodedData,ber] = NRZ_ExecuteDecoding(obj, eq_signal)
            % NRZ信号执行解码操作
            offset=0;
            % normalize
            sigRx_E = eq_signal - mean(eq_signal);
            sigRx_E = sigRx_E./max(sigRx_E);
            % make decision and convert into 0,1 sequence
            out = sign(sigRx_E-offset);
            decodedData = 0.5*(out+1);
            % 参考信号  重复 一定数量 ，满足解码 数量
            ref_seq=repmat(obj.Implementation.ref,1000,1);
            ref_seq=ref_seq(:);
            % NRZ 为 0,1 码元， 故解码为0,1码
            ref_seq=0.5*(ref_seq+1);

            % 解码
            [ber,num,~] = CalcBER(decodedData(obj.Nr.ncut_index:end),ref_seq(obj.Nr.ncut_index :end)); %计算误码率
            fprintf('Num of Errors = %d, BER = %1.7f\n',num,ber);
        end

        function [decodedData,ber] = PAM_ExecuteDecoding(obj, eq_signal)
            % 针对 PAM4 解码，计算误码
            % 量化区间
            A=[-2 0 2];
            % 参考信号  重复 一定数量 ，满足解码 数量
            ref_seq=repmat(obj.Implementation.ref,1000,1);
            ref_seq=ref_seq(:);
            % 参考序列
            [~,label] = quantiz(ref_seq,A,[-3,-1,1,3]);
            label_bit=obj.pam4demod(label);
            % 接收序列
            [~,I] = quantiz(eq_signal,A,[-3,-1,1,3]);
            decodedData=obj.pam4demod(I);
            % 解码
            [ber,num,~] = CalcBER(decodedData(obj.Nr.ncut_index:end),label_bit(obj.Nr.ncut_index:end)); %计算误码率
            fprintf('Num of Errors = %d, BER = %1.7f\n',num,ber);
        end

        function y= pam4demod(obj,sig)
            % 调制格式
            M= obj.signalPHY.M;
            % 解码类型
            y=zeros(1,length(sig)*2);
            %%PAM4
            for i = 1:length(sig)
                if sig(i) == -3
                    y(i*2-1) = 0;
                    y(i*2) = 0;
                elseif sig(i) == -1
                    y(i*2-1) = 1;
                    y(i*2) = 0;
                    %原始对应的是01
                elseif sig(i) == 1
                    y(i*2-1) = 1;
                    y(i*2) = 1;
                elseif sig(i) == 3
                    y(i*2-1) = 0;
                    y(i*2) = 1;
                    %原始对应的是10
                end
            end
            y = y(:);
        end


        function pam8demod= pam8demod(obj,sig)
            % 调制格式
            M= obj.signalPHY.M;
            % 解码类型
            pam8demod=zeros(1,length(sig)*3);
            k=1;
            %%PAM4
            for i = 1:length(sig)
                if sig(i) == -7
                    pam8demod(k) =0;
                    pam8demod(k+1)=0;
                    pam8demod(k+2)=0;
                elseif sig(i) == -5
                    pam8demod(k) =0;
                    pam8demod(k+1)=0;
                    pam8demod(k+2)=1;
                elseif sig(i) == -3
                    pam8demod(k) =0;
                    pam8demod(k+1)=1;
                    pam8demod(k+2)=1;
                elseif sig(i) == -1
                    pam8demod(k) =0;
                    pam8demod(k+1)=1;
                    pam8demod(k+2)=0;
                elseif sig(i) == 1
                    pam8demod(k) =1;
                    pam8demod(k+1)=1;
                    pam8demod(k+2)=1;
                elseif sig(i) == 3
                    pam8demod(k) =1;
                    pam8demod(k+1)=1;
                    pam8demod(k+2)=0;
                elseif sig(i) == 5
                    pam8demod(k) =1;
                    pam8demod(k+1)=0;
                    pam8demod(k+2)=0;
                elseif sig(i) == 7
                    pam8demod(k) =1;
                    pam8demod(k+1)=0;
                    pam8demod(k+2)=1;
                end
                k=k+3;
            end
        end

    end
end