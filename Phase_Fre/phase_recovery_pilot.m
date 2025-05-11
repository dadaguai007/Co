function [y,ph,ph_pilot] = phase_recovery_pilot(x,a,Parameters)
%PHASE_RECOVERY_PILOT 导频辅助的相位恢复算法
% 本函数通过盲相位搜索(BPS)结合最大似然(ML)估计，并利用周期性导频符号消除循环滑动，
% 实现高效相位恢复。适用于高阶QAM调制系统（16-QAM及以上）
% 参考论文：
% [1] 10.1109/LPT.2010.2049644 (BPS算法)
% [2] 10.1364/OFC.2014.Th4D.1 (导频辅助方案)

%% 参数提取 ===============================================================
N = Parameters.cpe_memory;        % CPE记忆深度：滑动平均窗口长度（符号数）
stp = Parameters.bps_anglestep;   % BPS角度搜索步长（弧度）
C = Parameters.C;                 % 星座图（列向量，如16-QAM复数序列）
pL = Parameters.cpe_pilotlen;     % 导频块长度（每个导频块包含的符号数）
bS = Parameters.cpe_blocksize;    % 导频间隔（两个导频块之间的符号数）
sps = Parameters.cpe_sps;         % 输入信号的过采样率（处理后下采样至1SPS）
ml_dec = Parameters.cpe_ml_dec;   % 是否使用ML判决（true-启用，false-硬判决）
Ntx = Parameters.Ntx;             % 发射序列数量（多输入场景）
sym = Parameters.cpe_sym;         % 星座对称性（2π, π, π/2等）
ph = 1;                           % 采样相位（用于2SPS下采样）


% 星座有效性检查（高阶调制要求）
if size(C,1)<=4
    error('导频相位恢复不适用于QPSK/BPSK，请改用Viterbi-Viterbi算法');
end

%% 导频模式检测 ===========================================================
if pL==0||bS==0           % 如果导频长度或间隔为0
    pilot = false;         % 禁用导频模式
else
    pilot = true;          % 启用导频模式
end

%% 信号预处理 =============================================================
x = x(ph:sps:end,:);       % 下采样至1SPS（选择指定采样相位）

%% 导频序列准备 ===========================================================
a = repmat(a,ceil(size(x,1)/size(a,1)),1); % 重复发送序列以匹配接收长度
a = a(1:size(x,1),:);      % 截断至接收信号长度
a_P = add_pilot(a,pL,bS);  % 在发送序列中插入导频符号（NaN填充非导频位置）

%% 初始化变量 =============================================================
M = round(sym/stp);        % 计算BPS测试角度数量（根据对称性和步长）
pt = (0:M-1).'/M*sym-sym/2;% 生成测试相位向量（[-sym/2, sym/2)范围）
ph = NaN(size(x));         % 预分配相位估计矩阵
ph_pilot = NaN(size(x,1),1);% 预分配导频相位向量

%% 主处理循环（逐通道处理）===============================================
for n = 1:size(x,2)        % 遍历每个接收通道（如X/Y偏振）
    xn = x(:,n);           % 提取当前通道信号

    %% 盲相位搜索(BPS)阶段 ==============================================
    % 生成候选相位旋转信号矩阵（每列对应不同相位补偿）
    x_rot = xn.*exp(1j*pt).';

    % 误差计算（滑动平均窗口）
    if ml_dec
        % ML判决模式：计算与最近星座点的平方距离
        er = movmean(abs(x_rot - ML_dec(x_rot,C)).^2,N);
    else
        % 硬判决模式：快速近似计算
        er = movmean(abs(x_rot - hard_dec(x_rot,C)).^2,N);
    end

    % 寻找最小误差对应的相位索引
    [~,idx] = min(er,[],2);

    %% 导频辅助循环滑动消除 =============================================
    if pilot
        % 获取当前通道对应的发射导频序列索引
        a_idx = mod(n-1,Ntx)+1;

        % 提取导频位置相位信息（基于已知导频符号）

        ph_pilot = angle(movmean(reshape(xn,size(a_P(:,a_idx),1),[]).*...
            conj(a_P(:,a_idx)),pL,'Endpoints','fill'));
        % 定位导频块位置
        pilot_pos = find(~isnan(ph_pilot(:,1)));

        % 相位解绕（消除2π跳变）
        ph_pilot(pilot_pos,:) = reshape(unwrap(...
            reshape(ph_pilot(pilot_pos,:),[],1)),size(pilot_pos,1),[]);

        % 线性插值填补非导频区域相位
        ph_pilot = interp_pilot_phase(ph_pilot(:));

        % BPS相位初步估计
        ph_bps = -pt(idx);

        % 循环滑动校正：将相位偏差对齐到导频参考相位
        ph_bps = ph_bps - sym*round((ph_bps - ph_pilot)/sym);
    else
        % 无导频模式：直接解绕BPS相位
        ph_bps = -unwrap(pt(idx)*2*pi/sym)/(2*pi/sym);
    end

    %% 最大似然(ML)相位细化 =============================================
    % 初步相位补偿信号
    y1 = xn.*exp(-1j*ph_bps);

    % ML相位估计（滑动平均提升稳定性）
    ph(:,n) = unwrap(angle(movmean(conj(y1).*xn,N)));
end

%% 应用最终相位补偿 ======================================================
y = x.*exp(-1j*ph);
end

%% 辅助函数1：硬判决 =====================================================
function c = hard_dec(x,C)
% 硬判决函数（针对标准QAM优化）
% 输入：
%   x - 待判决信号
%   C - 完整星座图
% 输出：
%   c - 判决后符号
switch log2(size(C,1)) % 根据星座大小优化判决
    case 2  % QPSK（实际代码禁用，此处保留结构）
        c = sign(real(x)) + 1j*sign(imag(x));
    case 4  % 16-QAM
        c = 2*(sign(real(x)) + 1j*sign(imag(x)));
        c = c + sign(real(x-c)) + 1j*sign(imag(x-c));
    case 6  % 64-QAM
        c = 4*(sign(real(x)) + 1j*sign(imag(x)));
        c = c + 2*(sign(real(x-c)) + 1j*sign(imag(x-c)));
        c = c + sign(real(x-c)) + 1j*sign(imag(x-c));
    case 8  % 256-QAM
        c = 8*(sign(real(x)) + 1j*sign(imag(x)));
        c = c + 4*(sign(real(x-c)) + 1j*sign(imag(x-c)));
        c = c + 2*(sign(real(x-c)) + 1j*sign(imag(x-c)));
        c = c + sign(real(x-c)) + 1j*sign(imag(x-c));
    otherwise % 通用ML判决
        c = ML_dec(x,C);
end
end

%% 辅助函数2：ML判决 =====================================================
function c = ML_dec(x,C)
% 最大似然判决（全搜索）
% 输入：
%   x - 待判决信号（可多列）
%   C - 星座图（列向量）
% 输出：
%   c - 判决符号（最小欧氏距离）
[~,idx] = min(abs(x.'-C)); % 矩阵化快速计算
c = C(idx);
c = reshape(c,size(x));    % 保持原始维度
end

%% 辅助函数3：导频相位插值 ===============================================
function x = interp_pilot_phase(x)
% 导频相位线性插值（处理NaN值）
% 算法来源：MATLAB Central #8225
x(1) = x(find(~isnan(x),1));    % 首值填充
x(end) = x(find(~isnan(x),1,'last')); % 末值填充
nan_pos = isnan(x);             % 定位缺失值
x(nan_pos) = interp1(find(~nan_pos),x(~nan_pos),find(nan_pos),'linear');
end

%% 辅助函数4：导频插入 ===================================================
function a_P = add_pilot(a,cpe_pilotlen,blocksize)
% 在发送序列中插入周期性导频
% 输入：
%   a - 原始发送序列
%   cpe_pilotlen - 导频块长度
%   blocksize - 导频间隔
% 输出：
%   a_P - 插入导频后的序列（非导频位置为NaN）

a_P = NaN(size(a)); % 初始化导频矩阵
% 生成导频位置索引（每个blocksize插入cpe_pilotlen个导频）
pilot_idx = reshape((0:blocksize:size(a,1)-blocksize)+(1:cpe_pilotlen)',[],1);
% 填充导频符号
for n = 1:size(a,2)
    a_P(pilot_idx,n) = a(pilot_idx,n);
end
end
