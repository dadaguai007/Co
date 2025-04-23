classdef DataProcessor < handle
    % 定义类的属性
    properties
        % 为了方便后续添加变量，采用结构体方式进行变量命名
        TxPHY; % 发射机参数
        signalPHY; % 接收机的信号参数
        Nr%无量纲参数
        Implementation;% 参考信号实施参数
        Button; % 开关讯号
    end

    % 定义类的方法
    methods
        % 构造函数
        function obj = DataProcessor(varargin)
            %初始化类的属性
            if numel(varargin) == 15
                obj.TxPHY                              = varargin{1} ;% 传递而来的发射机参数
                obj.signalPHY.fs                       = varargin{2} ;% 接收信号的采样率
                obj.signalPHY.fb                       = varargin{3} ;% 接收信号的波特率
                obj.signalPHY.fUp                      = varargin{4} ;% KK恢复算法的采样率
                obj.signalPHY.photocurrentSignal       = varargin{5} ;% 接收的光电流信号
                obj.Nr.squ_num                         = varargin{6} ;% 选取第 x 段信号
                obj.Nr.nTrainSym                       = varargin{7} ;% 训练序列长度
                obj.Nr.pilotIndex                      = varargin{8} ;% 相噪估计——导频位置
                obj.Nr.nTrainCarrier                   = varargin{9} ;% 频偏估计——导频数量，一般对全体载波进行估计
                obj.Implementation.qam_signal          = varargin{10} ;% 调制信号参考矩阵
                obj.Implementation.label               = varargin{11} ;% 同步参考信号
                obj.Button.CPE_Status                  = varargin{12};% 默认 关闭 CPE
                obj.Button.Fre_Offset_Status           = varargin{13}; % 频偏补偿 默认 打开
                obj.Button.selectAll                   = varargin{14}; % 是否选取全部信号 或者 分段选取
                obj.Button.select_photocurrentSignal   = varargin{15}; % 是否选择光电流信号进行处理
            end
            % 发射机自身存在缺陷，造成一些频率点（即载波）有固定差值，影响性能
            obj.signalPHY.errSub = [8,16,125,157,253];
            obj.Nr.k=1;% 默认为1 ，不变， 选取全部信号时，在函数中更改,此项控制信号的解码，均衡范围
            obj.Nr.CL=0.2;% 削波比例
            obj.Button.Clipping='off';% 默认为不削波
        end

        % 参考信号(解码使用)
        function createReferenceSignal(obj)
            % 函数体，用于创建参考信号
            qam_signal_ref=obj.Implementation.qam_signal;
            qam_signal_ref(obj.signalPHY.errSub,:)=[];
            % 并串转化
            ref_seq=reshape(qam_signal_ref,1,[]);
            % 得到用于解码的参考信号
            obj.Implementation.ref = repmat(ref_seq,1,100);
        end


        % 参考矩阵，用于信号的DSP处理中
        % 信号同步后，再进行使用
        function createReferenceSignal_matrix(obj)
            obj.Implementation.qam_signal_mat=repmat(obj.Implementation.qam_signal,1,obj.Nr.k);
        end


        % 对接收信号进行预滤波（低通滤波）
        function filteredData = preFilter(obj, input_data, cutoffFreq)
            % 截止频率应设置为22e9
            % 低通滤波
            signal_orgin = LPF(input_data,obj.signalPHY.fs,cutoffFreq);
            % 取实部
            filteredData = real(signal_orgin(1:2*floor(length(signal_orgin)/2)));
        end

        % 选取某段信号
        function selectedPortion=selectSignal(obj,Index_P,DataGroup)
            if obj.Nr.squ_num > length(Index_P)
                warning('选取序列超出范围')
                return;
            else
                % 选取数据段
            selectedPortion= DataGroup{obj.Nr.squ_num};
        
            end
        end

        % 信号同步
        function [DataGroup,Index_P,selectedPortion]=Synchronization(obj,Rxsig)
            % 同步算法好像需要先进行KK恢复，再进行同步
            %  参考
            label=obj.Implementation.label;
            % 信号与相位 ———— 同步
            [DeWaveform,P,~,~] = Quick_Syn_Vec(Rxsig,label,1/obj.signalPHY.fUp,1/obj.signalPHY.fb);
            % 已知P为每个同步点的位置，但是可能会出现偏差
            % 修正同步点坐标
            Index_P=obj.getSyncPoint(P,DeWaveform);

            % 对光电流信号进行同步
            if strcmp(obj.Button.select_photocurrentSignal,'on')
                % 赋值光电流信号
                photocurrent_data=obj.signalPHY.photocurrentSignal;
                [DeWaveform,~,~,~] = Quick_Syn_Vec(photocurrent_data,real(label),1/obj.signalPHY.fs,1/obj.signalPHY.fb);
            end


            % 数据进行存储
            DataGroup = cell(1, length(Index_P));
            for Idx=1:length(Index_P)
                Data = DeWaveform(Index_P(Idx):Index_P(Idx)+obj.TxPHY.len-1);
                DataGroup{Idx} = Data;
            end
            % 是否选取全部信号
            if strcmp(obj.Button.selectAll,'on')
                % 选取全部信号
                selectedPortion=DeWaveform(Index_P(1):Index_P(length(Index_P))+obj.TxPHY.len-1);
                % 后续 DSP 使用的数组长度进行变化
                obj.Nr.k=length(Index_P);
            else
                % 选取数据段
                selectedPortion=DeWaveform(Index_P(obj.Nr.squ_num):Index_P(obj.Nr.squ_num)+obj.TxPHY.len-1);
            end

            % 创建DSP所需的信号矩阵
            obj.createReferenceSignal_matrix();
        end


        % 纠正同步点位置
        function Index_P=getSyncPoint(obj,P,DeWaveform)

            for i=1:length(P)
                Index_P(i)=P(1)+(i-1)*obj.TxPHY.len;
                % 判断每组数据是否超出范围
                if Index_P(i)+obj.TxPHY.len-1>length(DeWaveform)
                    warning('某一组数据长度超出接收范围');
                    % 移除最后一个变量
                    Index_P(i:end)=[];
                    break;
                end
            end

        end

        % KK接收机
        function Rxsig=KK_receiver(obj,rxSig)

            % 对光电流信号或者对场信号进行处理
            if strcmp(obj.Button.select_photocurrentSignal,'on')
                Rxsig_up = KK_New(rxSig,obj.signalPHY.fb ,obj.signalPHY.fUp );

                Rxsig=resample(Rxsig_up,obj.signalPHY.fb,obj.signalPHY.fUp);
            else
                % KK 算法
                Rxsig = KK_New(rxSig,obj.signalPHY.fs ,obj.signalPHY.fUp );

            end
        end


        % 信号预处理
        function  [ReceivedSignal,Dc]=Preprocessed_signal(obj,rxsig)
            
            % 直流
            Dc=mean(rxsig);
            % 削波
            if strcmp(obj.Button.Clipping,'on')
                [rxsig,~]=clipping_mean(rxsig,obj.Nr.CL);
            end           
            % KK
            ReceivedSignal=obj.KK_receiver(rxsig);

        end

        % 输入QAM信号进行解码
        function [ber,num,L]=Direcct_Cal_BER(obj,ReceivedSignal)
            % 信号解码
            qam_bit=obj.hard_decision(ReceivedSignal);
            % label信号
            ref_seq =qamdemod(obj.Implementation.ref ,obj.TxPHY.M,'OutputType','bit','UnitAveragePower',1);
            ref_seq=ref_seq(:);
            % 计算误码率
            [ber,num,~] = Calc_BER(qam_bit,ref_seq);
            % 比特数量
            L = min(length(qam_bit),length(ref_seq));
            fprintf('Num of Errors = %d, BER = %1.7f\n',num,ber);
        end

        % 信号接收
        function [signal_ofdm_martix,data_ofdm_martix,Hf,data_ofdm] = OFDM_ExecuteDecoding(obj, ReceivedSignal)

            % 解OFDM
            signal_ofdm = reshape(ReceivedSignal,obj.TxPHY.fft_size+obj.TxPHY.nCP,[]); % 转换为矩阵形式
            signal_ofdm(1: obj.TxPHY.nCP,:) = [];   %去除CP
            signal_ofdm = fft(signal_ofdm);
            % 存储矩阵形式的OFDM 信号
            signal_ofdm_martix=signal_ofdm;
            % get the modulated carriers
            data_ofdm = signal_ofdm(obj.TxPHY.dataCarrierIndex,:);
            % 信道均衡
            [data_ofdm,Hf]=obj.one_tap_equalization(data_ofdm);

            % 频偏问题
            if strcmp(obj.Button.Fre_Offset_Status,'on')
                % 计算频偏
                freqOffset = obj.calculateFrequencyOffset(data_ofdm);
                %频偏补偿
                data_ofdm =obj.compensateFrequencyOffset(data_ofdm,freqOffset);
            end

            % 相位偏差消除
            if strcmp(obj.Button.CPE_Status,'on')
                % CPE compensation
                data_ofdm=obj.CPE_Eliminate(data_ofdm);
            end

            %保留信号矩阵
            data_ofdm_martix=data_ofdm;

            % 去除发射机异常频点（载波）
            data_ofdm(obj.signalPHY.errSub,:)=[];

            % 并串转换
            data_ofdm=data_ofdm(:); % 输出为QAM信号向量，不进行判决
            %data_ofdm = data_ofdm./sqrt(mean(abs(data_ofdm(:)).^2));
        end

        % ZF信道均衡
        function [data_kk,Hf]=one_tap_equalization(obj,data_ofdm)
            % channel estimation
            rxTrainSymbol = data_ofdm(:,1:1:obj.Nr.nTrainSym);
            qam_signal_mat=obj.Implementation.qam_signal_mat; % 参考信号矩阵
            refTrainSymbol = qam_signal_mat(:,1:1:obj.Nr.nTrainSym);
            % 信道响应
            Hf = mean(rxTrainSymbol./refTrainSymbol,2);
            % channel equalization
            data_kk = data_ofdm.*repmat(1./Hf,1,obj.TxPHY.nPkts*obj.Nr.k);
        end

        % ZF信道均衡-平滑操作
        function [data_kk,Hf]=smooth_one_tap_equalization(obj,data_ofdm)
            % channel estimation
            rxTrainSymbol = data_ofdm(:,1:1:obj.Nr.nTrainSym);
            qam_signal_mat=obj.Implementation.qam_signal_mat; % 参考信号矩阵
            refTrainSymbol = qam_signal_mat(:,1:1:obj.Nr.nTrainSym);
            % 信道响应
            Hf = mean(rxTrainSymbol./refTrainSymbol,2);
            Hf = smooth(Hf,5);  % 对信道进行平滑，增加估计准确度
            % channel equalization
            data_kk = data_ofdm.*repmat(1./Hf,1,obj.TxPHY.nPkts*obj.Nr.k);
        end

        % 相位均衡 CPE
        function data=CPE_Eliminate(obj,data_ofdm)
            % 提取导频信号
            phi_mean=angle(mean(data_ofdm(obj.Nr.pilotIndex,:)./...
                obj.Implementation.qam_signal_mat(obj.Nr.pilotIndex,:),1));

            % 补偿相噪
            data=data_ofdm.*...
                repmat(exp(-1j.*phi_mean),size(data_ofdm,1),1);
        end

        % 计算信号频偏
        function freqOffset = calculateFrequencyOffset(obj,data_ofdm)
            % 参考矩阵生成
            qam_signal_mat=obj.Implementation.qam_signal_mat;
            % 接收信号
            rxTrainSymbol_phase = data_ofdm(1:obj.Nr.nTrainCarrier,:);
            %参考信号
            refTrainSymbol_phase = qam_signal_mat(1:obj.Nr.nTrainCarrier,:);
            % 计算角度
            phi=angle(rxTrainSymbol_phase./refTrainSymbol_phase);
            % 对OFDM符号的相位偏差估计
            freqOffset=zeros(size(qam_signal_mat,2),1);
            for j=1:size(qam_signal_mat,2)
                K=zeros((obj.Nr.nTrainCarrier),1);
                M=zeros((obj.Nr.nTrainCarrier),1);
                for i=1:obj.Nr.nTrainCarrier
                    K(j)=i*phi(i,j);
                    M(j)=i.^2;
                end
                % 升成每个符号的频偏估计量
                freqOffset(j)=sum(K)./sum(M);
            end
        end

        % 频偏补偿
        function compensatedSignal =compensateFrequencyOffset(obj,data_ofdm,freqOffset)
            % 参考矩阵生成
            qam_signal_mat=obj.Implementation.qam_signal_mat;

            % 对每一个符号进行补偿
            for k=1:size(qam_signal_mat,2)
                for m=1:size(qam_signal_mat,1)
                    compensatedSignal(m,k)=data_ofdm(m,k)*exp(-1j*m*freqOffset(k));
                end
            end

        end


        % 硬判决
        function data_qam=hard_decision(obj,Receiver_data)
            data_qam = qamdemod(Receiver_data,obj.TxPHY.M,'OutputType','bit','UnitAveragePower',1);
            data_qam=data_qam(:);
        end

        % EVM测算
        function [rmsEVM_symbol,rmsEVM_subcarrier,rmsEVM_martix]=EVM_Measure_martix(obj,data_ofdm_martix)
            % 排除异常频点（载波）
            data_ofdm_martix(obj.signalPHY.errSub,:)=[];

            qam_signal_mat=obj.Implementation.qam_signal_mat(obj.signalPHY.errSub,:);
            qam_signal_mat(obj.signalPHY.errSub,:)=[];

            %EVM measure
            evm = comm.EVM(AveragingDimensions=1);
            rmsEVM_symbol = evm(data_ofdm_martix,qam_signal_mat);


            evm = comm.EVM(AveragingDimensions=2);
            rmsEVM_subcarrier = evm(data_ofdm_martix,qam_signal_mat);


            evm = comm.EVM(AveragingDimensions=[1 2]);
            rmsEVM_martix = evm(data_ofdm_martix,qam_signal_mat);

            fprintf('EVM = %1.7f\n',rmsEVM_martix);
        end

    end
end