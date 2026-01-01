function [Tu,Ts,Tg ] = effectivelength(x,fs,snr,N,K)
%**************************************************************************
%功能：估计信号的有效数据长度
%x:输入的信号
%fc:信号的频率，fs为信号的采样频率
%snr:信噪比
%N:信号每帧的符号数
%短波信号Tu的理想值3.2us,Ts的理想值4us
%**************************************************************************
 p = std(x);
 x = x/p;
 sig_1 = awgn(x,snr,'measured');
 L=length(x);
 t=1/fs;
 P=(L/N);
 xcorr_len=1;    %自相关长度,以OFDM符号为单位
 [Tu,Ts]=auto_xcorr(sig_1,P, xcorr_len, N,t,K);           %计算信号的自相关函数
 Tg=Ts-Tu;   %使用固定延时的fft估计CP长度
