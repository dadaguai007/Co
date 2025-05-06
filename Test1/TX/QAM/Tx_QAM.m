function [Stx,txSyms] = Tx_QAM(QAM,txBits)

% Last Update: 30/09/2018


%% Generate Symbols from Bits  �ӱ������ɷ���
switch QAM.encoding
    case 'normal'  %����QAM_config�����ã� Input Parser�е�encoding  ����ѡ��  ���磺encoding = 'normal' 
        txSyms = bit2sym(txBits,log2(QAM.M));%Bit-to-Symbol Mapping
        Stx = symbol2signal(txSyms,QAM.IQmap);
    case 'diff-quad'%��ֵ΢������ֵ����
        txSyms = bit2sym_DiffQuad(txBits,log2(QAM.M));
        Stx = symbol2signal(txSyms,QAM.IQmap);   
end