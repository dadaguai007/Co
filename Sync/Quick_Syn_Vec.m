% [DeWaveform,P,OptSampPhase]=Quick_Syn_Vec(RecWaveform,label,SampleT_Osci,SampleT_AWG,OptSampPhase)
% Origin version was got from Jane
% log: 2014.4.2 halfLen & quadLen are added to replace original fixed value
% log: 20170619 Vectorize the correlation calculation part to have 75%
%      running time reduction

function [DeWaveform,P,OptSampPhase,MaxCorrIndex]=Quick_Syn_Vec(RecWaveform,label,SampleT_Osci,SampleT_AWG,OptSampPhase)
%RecWaveform：接收到的波形数据
%label：一个参考波形
%SampleT_Osci：示波器采样时间间隔。
%SampleT_AWG：任意波形发生器（AWG）的采样时间间隔。
%OptSampPhase：可选的最优采样相位。
if nargin < 5
    OptSampPhase = [];
end

% find the optimal sampe phase
if isempty(OptSampPhase)
    seq_len = 1e5;
    % seq_len = ceil(length(RecWaveform)/2);
    SeqForSync = RecWaveform(1:seq_len);
    
    maxPhaseError = 0.05;
    if 1/maxPhaseError - floor(1/maxPhaseError) ~= 0
        error('Wrong maxPhaseError value!');
    end
    
    tSeqBefSamp=(0:numel(SeqForSync)-1).'*SampleT_Osci;
    %搜索的相位
    SampPhaseVec=(0:maxPhaseError:(1-maxPhaseError)).';
    % 重采样后的时间点，搜索最小相位值 与 信号点数间隔时间  的乘积
    % 对SeqForSync按不同相位重采样,得到AftSampWaveform
    % 实现了对信号的重采样,通过改变时间轴,可以获得不同相位下的采样信号。
    tSeqAftSamp = 0:maxPhaseError*SampleT_AWG:tSeqBefSamp(end);
    AftSampWaveform = interp1(tSeqBefSamp,SeqForSync,tSeqAftSamp,'spline');


    %根据maxPhaseError将AftSampWaveform划分为nSymbols个符号波形,每个符号波形对应一个采样相位。
    nSymbols = floor(length(AftSampWaveform)*maxPhaseError);
    AftSampWaveform = reshape(AftSampWaveform(1:nSymbols/maxPhaseError),1/maxPhaseError,[]).';
    % use "xcorr is equivalent to conv with a vector's flip"
    % 计算AftSampWaveform与本地参考信号label的卷积CorrMat。
    % 互相关是相当于数组倒置的卷积（互相关需要减去直流）
    CorrMat = convn(AftSampWaveform,flipud(label));
    CorrResult = max(abs(CorrMat)).';
    
    [~, MaxCorrIndex]=max(CorrResult);
%     拟合求取最佳采样相位而实现的二次曲线拟合。

    if MaxCorrIndex==1
        polyForQ=polyfit([SampPhaseVec(end);SampPhaseVec(1:2)+1],[CorrResult(end);CorrResult(1:2)],2);
        % warning('Optimum resampling phase near phase Zero');
    elseif MaxCorrIndex==numel(CorrResult)
        polyForQ=polyfit([SampPhaseVec(end-1:end);SampPhaseVec(1)+1],[CorrResult(end-1:end);CorrResult(1)],2);
        % warning('Optimum resampling phase near phase One');
    else
        polyForQ = polyfit(SampPhaseVec(MaxCorrIndex-1:MaxCorrIndex+1),CorrResult(MaxCorrIndex-1:MaxCorrIndex+1), 2);
    end

    % 拟合出CorrResult的峰值点相位,作为最佳采样相位OptSampPhase。
    %  最大互相关对应的相位OptSampPhase,这就是最佳采样相位
    %   一个是多项式微分，一个是求解多项式；多项式计算
    OptSampPhase=roots(polyder(polyForQ));
    OptCorr=polyval(polyForQ,OptSampPhase);
    OptSampPhase=OptSampPhase-floor(OptSampPhase);

end
% 用OptSampPhase对整个RecWaveform信号进行重采样,得到重采样后的信号AftSampWaveform。
% 计算AftSampWaveform与label的互相关Corr,找到主峰P,这表示找到了同步位置。
% 得到同步后的信号DeWaveform
tSeqBefSamp=(0:numel(RecWaveform)-1)*SampleT_Osci;
tSeqBefSamp=tSeqBefSamp.';
tSeqAftSamp=OptSampPhase*SampleT_AWG:SampleT_AWG:tSeqBefSamp(end);
tSeqAftSamp=tSeqAftSamp.';
AftSampWaveform=interp1(tSeqBefSamp,RecWaveform,tSeqAftSamp,'spline');

% 同步
Corr=abs(xcorr(label,AftSampWaveform-mean(AftSampWaveform)));
%在一个【】中赋值两个索引，都赋值为零
Corr([1:1000,length(AftSampWaveform)-1000+1:length(AftSampWaveform)+1000]) = 0;
figure(200);
plot(Corr)
close(200);
peak = max(Corr);
P = find(Corr>0.98*peak);
P = flipud(numel(AftSampWaveform)-P);
P(P<1) = [];
DeWaveform=AftSampWaveform;

end
