classdef CoherentTx < handle
    % To Do:
    % 添加延迟函数

    % 定义类的属性
    properties
        % 为了方便后续添加变量，采用结构体方式进行变量命名
        TxPHY % 发射机参数
        Nr %无量纲参数
        Button % 开关讯号
    end

    % 定义类的方法
    methods
        % 构造函数
        function obj = CoherentTx(varargin)
            %初始化类的属性
            if numel(varargin) == 14
                obj.TxPHY.fs                           = varargin{1} ;% 发射信号的采样率
                obj.TxPHY.fb                           = varargin{2} ;% 发射信号的波特率
                obj.TxPHY.order                        = varargin{3} ;% 随机信号的阶数
                obj.TxPHY.prbsOrder                    = varargin{4} ;% prbs码的阶数
                obj.TxPHY.M                            = varargin{5} ;% 调制格式
                obj.TxPHY.sps                          = varargin{6} ;% 每符号采样点
                obj.TxPHY.NSym                         = varargin{7} ;% 码元数目
                obj.TxPHY.Nmodes                       = varargin{8} ;% 偏振状态
                obj.Nr.psfShape                        = varargin{9} ;% 脉冲形式
                obj.Nr.psfRollOff                      = varargin{10} ;% 滚降系数
                obj.Nr.psfLength                       = varargin{11};% 影响长度
                obj.Nr.userShape                       = varargin{12};% 用户的成型滤波器
                obj.Button.DataType                    = varargin{13};% 选择模式
                obj.Button.shapingFilter               = varargin{14};% 成型滤波器的生成方式
            end
        end

        % 输出信号
        function [signalDualPol,qamDualPol]=dataOutput(obj)

            % 成型滤波器
            if strcmp(obj.Button.shapingFilter,'system')
                pulse = obj.systemHsqrt();
            elseif strcmp(obj.Button.shapingFilter,'user')
                pulse = obj.getPulse(obj.Nr.userShape);
            end

            % 是否为双偏振信号
            for indMode = 1:obj.TxPHY.Nmodes
                %生成bit数
                switch lower(obj.Button.DataType)
                    case 'prbs'
                        [~,symbols]=obj.prbs_bits();
                    case 'rand'
                        symbols=obj.rand_bits();
                end
                %调制
                qam_signal=obj.qam(symbols);
                % 装载信号
                qamDualPol(:, indMode) = qam_signal;
                % 脉冲成型
                filteredSignal=obj.applyShapingFilter(qamDualPol(:, indMode),pulse);
                % 装载成型后的信号
                signalDualPol(:,indMode)=filteredSignal;
            end

        end

        % 生成比特数
        function [data,symbols_prbs]=prbs_bits(obj)
            %参数：obj.prbsOrder，NSym，M
            %采用prbs码生成基本数据
            data = prbs1(obj.TxPHY.prbsOrder,obj.TxPHY.NSym*log2(obj.TxPHY.M),0);
            data_2bit=reshape(data,log2(obj.TxPHY.M),[]);
            %             symbols_prbs = 2.^(0:log2(obj.M)-1)*data_2bit;
            symbols_prbs=double(data_2bit);
        end
        % 生成比特数
        function symbols_rand=rand_bits(obj)
            rng(obj.TxPHY.order);
            %参数：obj.prbsOrder，NSym，M
            symbols_rand=randi([0,1],log2(obj.TxPHY.M),obj.TxPHY.NSym);
        end

        % 调制生成信号
        function qam_signal=qam(obj,symbols)
            qam_signal=qammod(symbols,obj.TxPHY.M,'InputType','bit','UnitAveragePower',1) ;
        end

        %设计根升余弦脉冲成型滤波器
        function hsqrt = systemHsqrt(obj)
            hsqrt = rcosdesign(obj.Nr.psfRollOff,obj.Nr.psfLength,obj.TxPHY.sps,obj.Nr.psfShape);
        end

        % 手动设计成型滤波器
        function pulse=getPulse(obj,pulseshape)
            % pulse shape
            if strcmp(pulseshape, 'nrz')
                pulse = pulseShape('nrz', obj.TxPHY.sps);
            elseif strcmp(pulseshape, 'rrc')
                pulse = pulseShape('rrc', obj.TxPHY.sps, 4096, obj.Nr.psfRollOff, 1/obj.TxPHY.fs);
            end
            pulse = pulse / max(abs(pulse));
        end


        % 应用脉冲成型滤波器
        function    filteredSignal=applyShapingFilter(obj,symbTx,pulse)
            % Upsampling
            symbolsUp = upsample(symbTx, obj.TxPHY.sps);
            % Pulse shaping
            if strcmp(obj.Button.shapingFilter,'system')
                filteredSignal=conv(symbolsUp,pulse,'same');
            elseif strcmp(obj.Button.shapingFilter,'user')
                filteredSignal = firFilter(pulse, symbolsUp);
            end
        end

        % 相噪建模
        function   Pin=phaseNoise(obj,sigTx,lw)
            % 相噪建模
            phi_pn_lo = phaseNoise(lw, length(sigTx), 1/obj.TxPHY.fs);
            sigLO = exp(1i * phi_pn_lo);
            Pin=sigLO;
        end

        % 激光器建模
        function   sigLO = laserMode(obj,inputSignal,lw,RIN,Plo_dBm)
            % 转换为W
            Plo=obj.dBTow(Plo_dBm);
            % LO
            paramLO=struct();
            paramLO.P = Plo;
            paramLO.lw = lw;          % laser linewidth
            paramLO.RIN_var = RIN; % 一般设置为0
            paramLO.Fs = obj.TxPHY.fs;
            paramLO.N = length(inputSignal);
            sigLO = basicLaserModel(paramLO);
        end

        % 添加频偏
        function signalWithOffset = addFrequencyOffset(obj, sigLO,freqOffset)
            % 创建时间轴
            [~,t_up]=freq_time_set(obj.TxPHY.NSym*obj.TxPHY.sps,obj.TxPHY.fs);
            signalWithOffset = sigLO.*exp(1j*2*pi*freqOffset*t_up); % add frequency offset
        end

        % dBm to W (功率转换)
        function Pch = dBTow(~,Pch_dBm)
            Pch = (10 .^ (Pch_dBm / 10)) * 1e-3;
        end

        % 调制后的信号，进行功率转换
        function  signal=setSignalPower(obj,input,Channel_power_type,Pin_dBm)
            % 转换为W
            Pin=obj.dBTow(Pin_dBm);
            % 每个通道的发射功率是相对于所有偏振模式总功率，不是相对于单个偏振模式的功率。
            % 将发射功率除以偏振模式的数量可以得到每个通道在单个偏振模式下的功率。
            if strcmp(Channel_power_type,'output')
                % 这句的作用是，输出的光功率 就是 设置的功率，不用进行放大了，也即前面设置的功率为输出的光功率
                signal = sqrt(Pin / obj.TxPHY.Nmodes) * pnorm(input);
            elseif strcmp(Channel_power_type,'input')
                % 这里的作用是， 前面设置的即为输入光功率为多少，
                signal = sqrt(Pin / obj.TxPHY.Nmodes)* (input) ;
            end
        end


        
    end
end