% 幅度分布统计函数
function [value,percent] = amplitude_distribution(input)
    [input_qtz,~,~,~] = Quantization(input,256);  % 256级量化
    distribution_table = tabulate(input_qtz);     % 统计分布
    value = distribution_table(:,1);             % 量化电平值
    percent = distribution_table(:,3)/100;       % 百分比转换
    figure
    plot(value,percent);                         % 绘制分布图
    title('amplitude_distribution');
end

% 均匀量化函数
function [output,distortion,codebook,partition] = Quantization(InputSignal,level_num)
    maximum = max(abs(InputSignal));              % 取信号绝对值最大值
    minimum = -maximum;                           % 对称负区间
    quantization_gap = (maximum-minimum)/(level_num-1); % 计算量化间隔
    codebook = minimum : quantization_gap : maximum;    % 生成码本
    partition = codebook(1:end-1) + quantization_gap/2;  % 生成判决电平
    [~,output,distortion] = quantiz(InputSignal,partition,codebook); % 执行量化
end