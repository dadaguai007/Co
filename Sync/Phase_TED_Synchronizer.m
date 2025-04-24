classdef Phase_TED_Synchronizer < handle
    %     loop_gain = 0.3
    % rollOff  = 0.2;
    %     根据TED类型（如ZCTED、GTED等）计算误差信号增益Kp。
    % 自定义函数calcTedKp和piLoopConstants用于生成环路滤波器参数K1/K2。
    %     发送信号 → 重采样（txResamp） → 信道延迟 → 添加噪声（rxSeq） → 匹配滤波（mfOut）

    properties
        sync_type     % 同步算法类型: 'gardner', 'MM', 'EL' / % TED (MLTED, ELTED, ZCTED, GTED, or MMTED)
        loop_gain     % 环路增益 控制环路收敛速度和稳定性的增益系数
        rollOff      % Pulse shaping roll-off factor （影响脉冲成形和S曲线形状）
        TED
        sps
        eta
        Bn_Ts
        intpl   % 插值方法
        rcDelay
        Ex
        kp
        k1
        k2
        systemParam
    end

    methods
        function obj = Phase_TED_Synchronizer(sync_type, loop_gain,rollOff,Ted,sps,intpl)
            % 构造函数
            obj.sync_type = sync_type;
            obj.loop_gain = loop_gain;
            obj.rollOff=rollOff;
            obj.TED=Ted;
            obj.sps=sps;
            obj.eta      = 1;          % 环路阻尼因子（稳定性控制）
            obj.Bn_Ts    = 0.01;       % 环路带宽 × 符号周期（控制同步速度）
            obj.intpl = intpl;
            obj.rcDelay  = 10;         % Raised cosine (combined Tx/Rx) delay
            obj.Ex       = 1;          % Average symbol energy
        end

        function [K1, K2] = process(obj)
            % 计算环路增益(计算定时误差检测器（TED）的增益 )
            Kp  = calcTedKp(obj.TED, obj.rollOff);
            % 两种方法计算增益：'analytic'（解析法）或 'simulated'（仿真法），默认为解析法

            % 计数器增益（类似VCO增益）
            K0 = -1;                            % 控制环路响应速度
            % PI : PI环路的控制增益
            [ K1, K2 ] = piLoopConstants(Kp, K0, obj.eta, obj.Bn_Ts, obj.sps);

            fprintf("Loop constants:\n");
            fprintf("K1 = %g; K2 = %g; Kp = %g\n", K1, K2, Kp);

            obj.kp=Kp;
            obj.k1=K1;
            obj.k2=K2;
        end

        function  SYMSYNC=systemProcess(obj,Kp)

            % MATLAB's symbol synchronizer for comparison
            tedMap = containers.Map({'ELTED', 'ZCTED', 'GTED', 'MMTED'}, ...
                {'Early-Late (non-data-aided)', ...
                'Zero-Crossing (decision-directed)', ...
                'Gardner (non-data-aided)', ...
                'Mueller-Muller (decision-directed)'
                });
            if strcmp(obj.TED, "MLTED")
                warning("MLTED not supported by MATLAB's synchronizer - using ZCTED");
                matlabTed = "ZCTED";
            else
                matlabTed = obj.TED;
            end
            SYMSYNC = comm.SymbolSynchronizer(...
                'TimingErrorDetector', tedMap(matlabTed), ...
                'SamplesPerSymbol', obj.sps, ...
                'NormalizedLoopBandwidth', obj.Bn_Ts, ...
                'DampingFactor', obj.eta, ...
                'DetectorGain', Kp);

            obj.systemParam=SYMSYNC;
        end


        % 相位同步
        function rxSync=symbolPhaseSync(obj,rxSeq,mfOut,const)
            % mfOut为匹配滤波后的信号，rxSeq为接收信号

            % debug 监视器
            debug_tl_static  = 0; % Show static debug plots after sync processing
            debug_tl_runtime = 0; % Open scope for debugging of sync loop iterations

            % 对参考信号进行归一化
            Ksym = modnorm(const, 'avpow', obj.Ex);
            const = Ksym * const;


            if strcmp(method,'process')
                %symbol timing recovery implementation
                [ rxSync ] = symbolTimingSync(obj.TED, obj.intpl, obj.sps, rxSeq, mfOut, obj.k1, obj.k2, ...
                    const, Ksym, obj.rollOff, obj.rcDelay, debug_tl_static, debug_tl_runtime);
            elseif strcmp(method,'system')
                % MATLAB's implementation
                rxSync = step(obj.systemParam, mfOut);
            end
        end


    end



end
