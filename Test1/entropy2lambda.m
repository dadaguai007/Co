function lambda = entropy2lambda(H,C)
%entropy2lambda     Convert entropy values to required shaping parameter,
%                   lambda using a pre-stored look-up table (LUT)
...ʹ��Ԥ�洢�Ĳ��ұ� (LUT) ����ֵת��Ϊ��������β��� lambda
%   This function performs the conversion of query entropy points, H, into
%   the correspondent shaping parameter, lambda, required by probabilistic
%   shaping. The probability distribution function is assumed to be 
%   exp(-lambda*abs(C).^2).
...�˺���ִ�н���ѯ�ص� H ת��Ϊ���������������Ӧ���β��� lambda�� ������ʷֲ�����Ϊ exp(-lambda*abs(C).^2)��
%   INPUTS:
%       H := entropy values [1 x nEntropy]
%       M := size of the QAM constellation [scalar]  %scalar������
%
%   OUTPUTS:
%       lambda := shaping parameter [1 x nEntropy]  ���β���
%
%
%   Examples:
%       lambda = entropy2lambda(3.5:0.1:4,64);
%
%
% Authors: Fernando Guiomar
% Last Update: 23/10/2017

%% Optimize lambda to Achieve Desired Entropy  �Ż� lambda ��ʵ���������
options = optimset('MaxFunEvals',1e4,'TolX',1e-4,'TolFun',1e-4,...
    'Display','none','PlotFcns',[]);  %optimset���Ż���  TolX��TolFun���޸ľ������� ��Ϊ�Ż����������ò����ĺ���
[lambda,err] = fminsearch(@(lambda) fitMaxwellBoltzman(lambda,C,H),...
    0,options);             %�����������
...                          ��õĺ����������Ϊfit���﷨Ϊ
...                          [��Ͻ�� ��Ͼ���]��fit��X���ݣ�Y���ݣ���������͡���

end

%% Maxwell-Boltzman Fitting Function  Maxwell-Boltzmann ��Ϻ���
function err = fitMaxwellBoltzman(lambda,C,H)
    symProb = exp(-lambda*abs(C).^2);
    symProb = symProb/sum(symProb);
    entropy = -sum(symProb.*log2(symProb));%MB�ֲ���Px��symProb����    
    err = abs(H-entropy);  %errΪ��Ͼ���
end
