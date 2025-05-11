function [SNR_dB,Ps,Pn] = SNR_estimateFromSpectrum(PSD,Fs,Rs,rollOff,fNoise,fCenter,debugPlots)
%SNR_estimateFromSpectrum 从功率谱估计信噪比(SNR)
% 本函数通过分析信号的功率谱密度(PSD)，计算信号与噪声功率比
% 核心功能：
%   1. 信号功率计算（积分信号带宽内功率）
%   2. 噪声功率计算（积分指定噪声带功率）
%   3. SNR计算与可视化
% 输入参数：
%   PSD       - 功率谱密度估计（单位：W/Hz），行向量[1×nSamples]
%   Fs        - 采样频率（Hz）或频率向量
%   Rs        - 符号率（Baud）
%   rollOff   - 升余弦滚降系数
%   fNoise    - 噪声频率区间（相对或绝对频率），矩阵[nNoiseBins×2]
%   fCenter   - 信号中心频率（Hz），默认为0
%   debugPlots- 是否生成调试图，默认为false
% 输出参数：
%   SNR_dB    - 信噪比（dB）
%   Ps        - 信号功率（W）
%   Pn        - 噪声功率（W）

%% 输入参数解析 ======================================================
nSamples = size(PSD,2); % 获取PSD样本数

% 生成频率向量
if numel(Fs) == 1
    % 当Fs为标量时，生成对称频率轴
    f = (-nSamples/2:nSamples/2-1)*(Fs/nSamples); 
else
    % 当Fs为向量时，直接使用作为频率轴
    f = Fs;
    Fs = (max(f) - min(f)); % 计算实际采样率
end

% 设置默认参数
if nargin < 6
    fCenter = 0; % 默认中心频率为0Hz
end
if nargin < 5
    % 默认噪声带：[-0.9Fs/2, -0.8Fs/2]和[0.8Fs/2, 0.9Fs/2]
    fNoise = [-0.9 -0.8; 0.8 0.9]*Fs/2 + fCenter; 
else
    % 若输入为相对频率（<=1），转换为绝对频率
    if max(fNoise) <= 1
        fNoise = fNoise*Fs/2 + fCenter;
    end
end
if nargin < 7
    debugPlots = false; % 默认不显示调试图
end

%% 参数验证 ========================================================
nNoiseBins = size(fNoise,1); % 噪声区间数量
sigBW = Rs*(1+rollOff); % 计算信号带宽（考虑滚降系数）

% 检查噪声带是否与信号带宽重叠
if any(abs(fNoise) <= sigBW/2)
    warning('噪声带与信号带宽重叠，SNR估计可能不准确！');
end

%% 计算信号+噪声总功率 ==============================================
% 确定信号带宽边界索引
[~,indSig(1)] = min(abs(f - (fCenter - sigBW/2))); % 左边界
[~,indSig(2)] = min(abs(f - (fCenter + sigBW/2))); % 右边界

% 计算实际积分带宽
Psn_BW = abs(f(indSig(2)) - f(indSig(1)));

% 梯形积分计算信号带内总功率（含噪声）
Psn_SCM = trapz(f(indSig(1):indSig(2)), PSD(indSig(1):indSig(2)))...
    * sigBW/Psn_BW; % 带宽归一化

%% 计算噪声功率 ====================================================
[Pn_BW, Pn_bins] = deal(NaN(1,nNoiseBins)); % 预分配内存
indNoise = NaN(nNoiseBins,2); % 存储噪声带索引

for n = 1:nNoiseBins
    % 查找噪声带边界索引
    [~,indNoise(n,1)] = min(abs(f - fNoise(n,1)));
    [~,indNoise(n,2)] = min(abs(f - fNoise(n,2)));
    
    % 计算当前噪声带带宽
    Pn_BW(n) = abs(f(indNoise(n,2)) - f(indNoise(n,1)));
    
    % 梯形积分计算当前噪声带功率
    Pn_bins(n) = trapz(f(indNoise(n,1):indNoise(n,2)),...
        PSD(indNoise(n,1):indNoise(n,2)));
end

% 噪声功率归一化（按符号率带宽）
Pn = sum(Pn_bins)/sum(Pn_BW)*Rs;

%% 计算SNR ========================================================
Ps = Psn_SCM - Pn; % 信号功率 = 总功率 - 噪声功率
SNR_dB = 10*log10(Ps/Pn); % 转换为分贝值

%% 调试绘图 ========================================================
if debugPlots
    % 定义颜色方案
    blue = [0 0.2400 0.4310];    % 主频谱颜色
    red = [0.7300 0.0700 0.1690]; % 噪声带颜色
    
    % 自动选择频率单位
    if Fs > 1e12
        units = 'THz'; f = f*1e-12; Fs = Fs*1e-12;
    elseif Fs > 1e9
        units = 'GHz'; f = f*1e-9; Fs = Fs*1e-9;
    elseif Fs > 1e6
        units = 'MHz'; f = f*1e-6; Fs = Fs*1e-6;
    elseif Fs > 1e3
        units = 'KHz'; f = f*1e-3; Fs = Fs*1e-3;
    else
        units = 'Hz';
    end
    
    % 归一化PSD
    PSD = PSD/max(PSD);
    
    % 创建新图窗
    figure();
    
    % 绘制完整频谱
    hPlot1 = plot(f,10*log10(PSD));
    hold on;
    
    % 绘制噪声带
    for n = 1:nNoiseBins
        ind = indNoise(n,1):indNoise(n,2);
        hPlot2(n) = plot(f(ind'),10*log10(PSD(ind')));
    end
    
    % 绘制信号带
    ind = indSig(1):indSig(2);
    hPlot3 = plot(f(ind'),10*log10(PSD(ind')));
    
    %% 图形样式设置
    % 线型颜色设置
    set(hPlot1,'Color',blue); % 主频谱蓝色
    set(hPlot2,'Color',red,'LineStyle',':','LineWidth',1.5); % 噪声带红色虚线
    set(hPlot3,'Color',[0 0 0],'LineStyle','--','LineWidth',1.5); % 信号带黑色虚线
    
    % 坐标轴标签
    xlabel(['$f$ [' units ']'],'Interpreter','latex','FontSize',11);
    ylabel('归一化PSD [dB]','Interpreter','latex','FontSize',11);
    
    % 刻度设置
    xAxis = get(gca,'xaxis');
    set(xAxis,'TickLabelInterpreter','latex','FontSize',11);
    yAxis = get(gca,'yaxis');
    set(yAxis,'TickLabelInterpreter','latex','FontSize',11);
    
    % 网格线设置
    grid on;
    set(gca,'GridLineStyle','--','XMinorTick','off','XMinorGrid','off');
    axis tight;
    
    % 图例
    hLeg = legend([hPlot1,hPlot2(1),hPlot3],...
        '完整频谱','噪声功率','信号功率');
    set(hLeg,'Interpreter','latex','fontsize',9,'Location','NorthWest');
    
    % 添加SNR标注
    hText = text(f(end)-0.02*Fs,-1,['SNR = ' num2str(SNR_dB,'%1.1f') ' dB']);
    set(hText,'HorizontalAlignment','right','Color',blue,...
        'VerticalAlignment','top','Interpreter','latex','FontSize',12,...
        'BackgroundColor','w','EdgeColor',blue,'Margin',2);
end
end
