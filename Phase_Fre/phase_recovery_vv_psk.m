function [pr,errphi_filt,ph_pilot] = phase_recovery_vv_psk(x,a,Parameters)
% PSK信号的Viterbi&Viterbi相位恢复算法
% 本函数实现PSK信号的相位恢复，结合导频符号消除循环滑动现象
% 主要功能：
%   1. M次方相位估计（经典V&V算法）
%   2. 导频辅助的循环滑动消除
%   3. 相位误差平滑与补偿
% 参考论文：
% [1] 10.1109/LPT.2010.2049644 (V&V算法)
% [2] 10.1364/OFC.2014.Th4D.1 (导频辅助方案)

%% ========================== 参数提取 ================================
N = Parameters.cpe_memory;        % 相位估计记忆深度（滑动窗口长度）
pL = Parameters.cpe_pilotlen;     % 导频块长度（每个导频块包含的符号数）
bS = Parameters.cpe_blocksize;    % 导频间隔（正常数据块长度）
sps = Parameters.cpe_sps;         % 输入信号过采样率（1或2）
Ntx = Parameters.Ntx;             % 发射序列数量（多通道场景）
C = Parameters.C;                 % 星座图（PSK调制）

%% ========================== 衍生参数计算 ============================
Mapsk = 2*pi/Parameters.cpe_sym;  % M-PSK的M值（根据星座对称性计算）
ph = 1;                           % 下采样相位选择（2SPS时使用）
phi0 = min(abs(angle(C)));        % 获取星座最小相位角（用于相位对齐）

%% ========================== 导频模式检测 ============================
Ntx = size(a,2); % 实际发射通道数（覆盖参数中的设置）
if pL==0||bS==0  % 检测是否需要导频
    pilot = false; 
else
    pilot = true;
end

%% ========================== 信号预处理 ==============================
x = x(ph:sps:end,:); % 下采样至1SPS（符号率采样）
% 示例：当sps=2时，保留第1,3,5...个样本

%% ========================== 导频序列准备 ============================
a = repmat(a,ceil(size(x,1)/size(a,1)),1); % 循环扩展发送序列
a = a(1:size(x,1),:);         % 精确截断至接收信号长度
a_P = add_pilot(a,pL,bS);      % 插入导频符号（非导频位置置NaN）

%% ========================== M次方相位估计 ==========================
yM = sign(x.^Mapsk);  % M次方运算+符号提取（消除PSK调制）
% 原理：对于M-PSK，x^M将星座点映射到单位圆上同一角度

%% ========================== 相位误差估计 ============================
yM = movmean(yM,N);   % 滑动平均滤波（抑制噪声）
errphi_filt = angle(yM); % 提取相位误差（未解绕）

%% ========================== 导频辅助处理 ============================
if pilot
    % 相位预处理：移除初始相位偏移
    errphi_filt = errphi_filt - phi0*Mapsk;
    
    % 逐通道处理
    for n = 1:Ntx
        % 获取当前通道的导频位置相位
        a_idx = mod(n-1,Ntx)+1;  % 导频序列索引计算
        xn = x(:,n);             % 当前通道接收信号
        
        % 导频相位计算（基于已知符号）
        ph_pilot = angle(movmean(reshape(xn,size(a_P(:,a_idx),1),[]).*...
            conj(a_P(:,a_idx)),pL,'Endpoints','fill'));  
        
        % 定位有效导频块
        pilot_pos = find(~isnan(ph_pilot(:,1)));
        
        % 相位解绕（消除2π跳变）
        ph_pilot(pilot_pos,:) = reshape(unwrap(...
            reshape(ph_pilot(pilot_pos,:),[],1)),size(pilot_pos,1),[]);
        
        % 线性插值得到连续相位参考
        ph_pilot = interp_pilot_phase(ph_pilot(:));
        
        % 循环滑动校正（2π整数倍校正）
        errphi_filt(:,n) = errphi_filt(:,n) - ...
            2*pi*round((errphi_filt(:,n)-Mapsk*ph_pilot)/2/pi);
    end
else
    % 无导频模式：直接解绕相位
    errphi_filt = unwrap(errphi_filt)-phi0*Mapsk;
    ph_pilot = NaN(size(errphi_filt)); % 导频相位无效
end

%% ========================== 相位归一化 ==============================
errphi_filt = errphi_filt/Mapsk; % 恢复实际相位偏差（除以M次方系数）

%% ========================== 相位补偿 ================================
pr = x .* exp(-1j*errphi_filt); % 应用相位校正
end

%% ====================== 辅助函数1：导频插入 ========================
function a_P = add_pilot(a,cpe_pilotlen,blocksize)
% 在发送序列中插入周期性导频
% 输入：
%   a - 原始序列
%   cpe_pilotlen - 导频块长度
%   blocksize - 导频间隔
% 输出：
%   a_P - 带导频的序列（非导频位置为NaN）

a_P = NaN(size(a)); % 初始化导频矩阵
% 生成导频索引（每个blocksize插入cpe_pilotlen个导频）
pilot_idx = reshape((0:blocksize:size(a,1)-blocksize)+(1:cpe_pilotlen)',[],1);
% 填充导频符号
for n = 1:size(a,2)
    a_P(pilot_idx,n) = a(pilot_idx,n);
end
end

%% ================== 辅助函数2：导频相位插值 ========================
function x = interp_pilot_phase(x)
% 导频相位线性插值（处理NaN值）
% 来源：MATLAB Central文件交换#8225

% 边界处理
x(1) = x(find(~isnan(x),1));    % 首值填充（避免起始NaN）
x(end) = x(find(~isnan(x),1,'last')); % 末值填充

% 线性插值
nan_pos = isnan(x);             % 定位缺失位置
x(nan_pos) = interp1(...
    find(~nan_pos),x(~nan_pos),find(nan_pos),'linear');
end
