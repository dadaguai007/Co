classdef CoherentRx < handle

    % 定义类的属性
    properties
        % 为了方便后续添加变量，采用结构体方式进行变量命名
        TxPHY % 发射机参数
        signalPHY; % 接收机的信号参数
        Nr %无量纲参数
        Implementation
        Button % 开关讯号
    end

    methods
        % 构造函数
        function obj = CoherentRx(varargin)
            %初始化类的属性
            if numel(varargin) == 7
                obj.TxPHY                              = varargin{1} ;% 发射机参数
                obj.signalPHY.fs                       = varargin{2} ;% 接收机采样率
                obj.signalPHY.sps                      = varargin{3} ;% Dsp的上采样因子
                obj.signalPHY.hsqrt                    = varargin{4} ;% 匹配滤波器
                obj.signalPHY.photocurrentSignal       = varargin{5} ;% 接收的光电流信号
                obj.Implementation.qam_signal          = varargin{6} ;% 参考信号
                obj.Button.selectAll                   = varargin{7}; % 是否选取全部信号 或者 分段选取
            end
            obj.Button.delayRealSignal = 'on';     % 延迟实部信号
        end

        % 参考信号(解码使用)
        function createReferenceSignal(obj)
            reff=obj.Implementation.qam_signal;
            if obj.TxPHY.Nmodes == 2
                for inMode=1:obj.TxPHY.Nmodes
                    % 函数体，用于创建参考信号
                    qam_signal_ref(:,inMode)=reff(:,inMode);
                    % 得到用于解码的参考信号
                    reeff(:,inMode) = repmat(qam_signal_ref(:,inMode),1,100);
                end
                obj.Implementation.ref=reeff;
            else
                % 函数体，用于创建参考信号
                qam_signal_ref=obj.Implementation.qam_signal;
                % 得到用于解码的参考信号
                obj.Implementation.ref = repmat(qam_signal_ref,1,100);
            end
        end

        % 创建参考星座图
        function creatReferenceConstellation(obj)
            % Constellation
            z = (0:obj.TxPHY.M-1)';
            obj.Implementation.constellation = qammod(z,obj.TxPHY.M);
            % 是否需要进行
        end


        % 延迟函数：
        function outSignal=delaySignal(obj,input,delay_sample)

            if strcmp(obj.Button.delayRealSignal,'on')
                signalReal=real(input);
                outSignalImag=imag(input);
                % 延迟样本数
                outSignalReal = delay_signal(signalReal, delay_sample);
                outSignal=outSignalReal+1j*outSignalImag;
            else
                signalImag=imag(input);
                outSignalReal=real(input);
                % 延迟样本数
                outSignalImag = delay_signal(signalImag, delay_sample);
                outSignal=outSignalReal+1j*outSignalImag;
            end
        end


        % 创建RSOP
        function outSignal=addRsop(obj,input)

            num_signal=length(input);
            [~, t] = freq_time_set(num_signal, obj.signalPHY.fs);

            % x,y偏振赋值
            inx=input(1,:);
            iny=input(2,:);

            % 初始参数
            N =12;
            w_alpha = 5e3; %krad/s---rad/s
            w_phi = 13e3;
            w_kappa =11e3;

            U = cell(1,num_signal);

            % 三参量随机值，模拟rsop
            kappa0 = pi+randn(N,1)*2*pi;
            alpha0 = pi+randn(N,1)*2*pi;
            phi0 = pi+randn(N,1)*2*pi;

            for a = 1:num_signal
                alpha = w_alpha.*t(a)+alpha0;
                phi = w_phi.*t(a)+phi0;
                kappa = w_kappa.*t(a)+kappa0;
                oo11 = cos(kappa).*exp(1j*(alpha));
                oo12 = -sin(kappa).*exp(-1j*(phi));
                oo21 = -conj(oo12);
                oo22 = conj(oo11);
                y=1;
                for b = 1:N
                    y = y*[oo11(b) oo12(b);oo21(b) oo22(b)];
                end
                U{a} = y;
                rx(a) = U{a}(1,1).*inx(a)+U{a}(1,2).*iny(a);
                ry(a) = U{a}(2,1).*inx(a)+U{a}(2,2).*iny(a);
            end
            % 输出信号赋值
            outSignal(:,1)=rx;
            outSignal(:,2)=ry;
        end

        % 创建PMD模型
        function outSignal=addPMD(obj,input,tau_mean,sigma,w)
            % PMD在每个频率段有不同的延时，一般情况下，w只取一个值
            % 需要判断输入信号是否为2*N的形式

            N=12;%建模段数
            % tau_mean = 3; % average tau 3ps
            % sigma = 0.01; % variance
            % 每个光纤段生成一个随机的初始相位，模拟光脉冲在光纤中的随机偏振态变化
            % PMD的延迟
            Tau = sqrt(3*pi/(8*N))*(1+sigma*randn(N,1)).*tau_mean; %ti的取值,分为N段

            % 输入信号长度和时间轴
            num_signal=length(input);
            [~, t] = freq_time_set(num_signal, obj.signalPHY.fs);


            % RSOP 随时间变化，初始参数
            w_arfa = 111e3; %rsop w
            w_kappa = 61e3;
            w_phi = 115e3;
            kapa0 = pi+randn(N+1,1)*2*pi;
            arfa0 = pi+randn(N+1,1)*2*pi;
            phi0 = pi+randn(N+1,1)*2*pi;
            % 频率轴
            %             w_start = 1905;%1594nm            取中心波长的值
            %             w_stop = 1932;%1530nm             取中心波长加上信号的频宽的值
            %             w_spacing = 0.03; % angular freq spacing, unit: Trad/s       是不是应该取fs/length(signal)的长度
            %             w = w_start:w_spacing:w_stop; % unit: T rad/s


            % 随时间变化的RSOP矩阵，每个时间下，应为N+1 长度
            % 横轴为时间T，纵轴为段数N
            for i=1:length(t)

                arfa(:,i) = w_arfa.*t(i)+arfa0;
                phi(:,i) = w_phi.*t(i)+phi0;
                kapa(:,i) = w_kappa.*t(i)+kapa0;

                oo11(:,i)  = cos(kapa(:,i)).*exp(1j*(arfa(:,i)));
                oo12(:,i)  = -sin(kapa(:,i)).*exp(-1j*(phi(:,i)));
                oo21(:,i)  = -conj(oo12(:,i));
                oo22(:,i)  = conj(oo11(:,i));

            end

            % PMD矩阵，随频率变化
            % 横轴为频率w，纵轴为光纤段数
            for j=1:length(w)
                ee(:,j) = exp(1j*w(j).*Tau/2); % 群延时矩阵
                ff (:,j)= conj(ee(:,j));
            end

            % PMD动态建模
            % U_i 为横轴为时间T，纵轴为频率w
            omega=cell(1,length(w));
            for j=1:length(w)
                for K=1:N
                    for i=1:length(t)
                        if K==1
                            % 输入信号
                            U=input(:,i);
                            H=[oo11(K,i),oo12(K,i);oo21(K,i),oo22(K,i)];
                            % RSOP在时域上处理
                            U1(:,i)=H*U;
                        else
                            % 输入信号
                            U=UU(:,i);
                            H=[oo11(K,i),oo12(K,i);oo21(K,i),oo22(K,i)];
                            % RSOP在时域上处理
                            U1(:,i)=H*U;
                        end
                    end
                    % PMD矩阵
                    B=[ee(K,j),0;0,ff(K,j)];
                    % PMD在频率上处理,分偏振方向进行时频转换
                    for index=1:2
                        U_f(index,:)=fft(U1(index,:));
                    end
                    % PMD在频率上处理
                    U1=ifft(B*U_f);
                    % 传递变量
                    UU=U1;
                end
                % 多一段的rsop
                for i=1:length(t)
                    % 输入信号
                    U=UU(:,i);
                    H_N1=[oo11(N+1,i),oo12(N+1,i);oo21(N+1,i),oo22(N+1,i)];
                    outSignal(:,i)=U*H_N1;
                end
                omega{j}=outSignal;
            end

        end





        % 应用skew
        function s_skewed=addSkew(obj,input,delay)
            if strcmp(obj.Button.delayRealSignal,'on')
                % 应用I路延迟（Q路不变）
                % delay = 50e-12;    % 50ps延迟
                s_skewed = iqdelay(input, obj.signalPHY.fs, delay);
            end
        end


        %应用jitter,频偏,量化噪声
        function outSignal=applyADC(obj,input,ppm)

            % adc的采样频偏
            Fs_adc = 2*obj.TxPHY.fb*(1 + ppm/1e6);
            ppm_meas = (Fs_adc-2*obj.TxPHY.fb)/(2*obj.TxPHY.fb)*1e6;
            fprintf('ADC sampling rate = %.5f GS/s\n',Fs_adc/1e9);
            fprintf('ADC sampling clock drift = %.2f ppm\n',ppm_meas);
            paramADC = struct();
            paramADC.Fs_in =  obj.signalPHY.fs  ;
            paramADC.Fs_out = Fs_adc;
            paramADC.jitter_rms = ppm_meas; %400e-15 % 代表抖动
            paramADC.nBits =  8;  % 8比特量化
            paramADC.Vmax = max(real(input));
            paramADC.Vmin = min(real(input));
            paramADC.AAF = 'on';
            paramADC.N = 1001;
            % adc
            outSignal = adc(input, paramADC);

        end


        % 信号衰减
        function outSignal=Attenuator(~,input,Att_dB)
            % 衰减，以dB为单位值
            outSignal=input*10.^(-Att_dB/20);
        end


        % 信号传输
        function outSignal=signalTran(obj,input,param)

            if obj.TxPHY.Nmodes == 2
                % 双偏振传输
                outSignal = manakovssfm(input, param);
            else
                % 单偏振传输
                outSignal=ssfm(input,param);
            end

        end



        % 相干接收
        function outSignal=coherentReceive(obj,input,LO,theta,paramPD)
            % 判断信号是双偏振 or 单偏振
            if obj.TxPHY.Nmodes == 2
                outSignal=pdmCoherentReceiver(input, LO, theta, paramPD);
            else
                outSignal=coherentReceiver(input,LO,paramPD);
            end

        end


        % 匹配滤波
        function outSignal=matchFiltering(obj,input)
            % 创建变量
            outSignal=zeros(size(input));
            for inMode=1:length(obj.TxPHY.Nmodes)
                input(:,indMode)=input(:,indMode)-mean(input(:,indMode));
                outSignal(:,indMode)=conv(input(:,indMode),obj.signalPHY.hsqrt ,'same');
                outSignal(:,indMode)=pnorm(outSignal(:,indMode));
            end
        end

        % 色散补偿
        function outSignal=cdCompensation(obj,input,paramEDC)
            paramEDC.Fs =obj.signalPHY.fs;
            outSignal = cdc(input, paramEDC);
        end

        % 硬判决
        function data_qam=hard_decision(obj,Receiver_data)
            data_qam = qamdemod(Receiver_data,obj.TxPHY.M,'OutputType','bit','UnitAveragePower',1);
            data_qam=data_qam(:);
        end
        % 误码计算
        function [ber,num,L]=Direcct_Cal_BER(obj,ReceivedSignal)
            % 信号解码
            qam_bit=obj.hard_decision(ReceivedSignal);
            % label信号
            ref_seq =qamdemod(obj.Implementation.ref ,obj.TxPHY.M,'OutputType','bit','UnitAveragePower',1);
            ref_seq=ref_seq(:);
            % 计算误码率
            [ber,num,~] = CalcBER(qam_bit,ref_seq);
            % 比特数量
            L = min(length(qam_bit),length(ref_seq));
            fprintf('Num of Errors = %d, BER = %1.7f\n',num,ber);
        end

        % EVM
        function EVM = getEVM(obj,input)
            input = input ./ sqrt(bandpower(input));
            % 归一化
            constellation = obj.Implementation.constellation / sqrt(bandpower(obj.Implementation.constellation));
            Pe = bandpower(input - decision(input, constellation));
            Pref = bandpower(input);
            EVM = sqrt(Pe./Pref);
        end




    end
end