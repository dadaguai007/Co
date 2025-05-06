function [Stx,txSyms] = Tx_QAM(QAM,txBits)

% Last Update: 30/09/2018


%% Generate Symbols from Bits  从比特生成符号
switch QAM.encoding
    case 'normal'  %用于QAM_config（配置） Input Parser中的encoding  类型选择  例如：encoding = 'normal' 
        txSyms = bit2sym(txBits,log2(QAM.M));%Bit-to-Symbol Mapping
        Stx = symbol2signal(txSyms,QAM.IQmap);
    case 'diff-quad'%数值微分与数值积分
        txSyms = bit2sym_DiffQuad(txBits,log2(QAM.M));
        Stx = symbol2signal(txSyms,QAM.IQmap);   
end