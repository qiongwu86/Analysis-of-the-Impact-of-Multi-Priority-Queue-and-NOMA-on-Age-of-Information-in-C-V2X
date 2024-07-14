function [RBsBeacon,gammaMin_dB,NbitsHz] = findRBsBeaconSINRmin(MCS,BeaconSizeBits)
% This function calculates RBs per beacon and minimum required SINR

% Call function to find ITBS value from MCS���Ʊ��뷽��
ITBS = findITBS(MCS);

% Call function to find the modulation format (number of bits per symbol)����������M
Nbps = findModulation(MCS);

% Call function to find the number of RBs per beacon
[RBsBeacon,Nbits] = findRBsBeaconNbits(ITBS,BeaconSizeBits);

% Compute the effective code rate����=����������ʣ�λ�����ŵ�����������ʵı�ֵ
CR = Nbits/((RBsBeacon/2)*9*12*Nbps);   %2Blog2(1+SINR)    Nbit/bps = N Hz

% Compute spectral efficiency (bits/Hz)  =  CR/B  ,��ΪRth
%ÿ����ÿ���ȵ�Ƶ�״����ڿ��Դ���ı�����
NbitsHz = (12*14*Nbps*CR)/(1e-3*180e3);      

% Compute the minimum required SINR
% (alfa is taken from 3GPP)
alfa = 0.4;
gammaMin_dB = 10*log10(2^(NbitsHz/alfa)-1);      
% NbitsHz = alfa*log2(1+SINR) �� SINR=2^(NbitsHz/alfa)-1 ����ת��ΪdB

%%%%%%%%Ŀǰ���ȼ����1��֡���Դ�������������������Ӧ��SINR�������ʴ�С��Ӱ��SINR��ֵ
end