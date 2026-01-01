%实现基本的OFDM通信系统仿真，包括QPSK调制、串并变换、IFFT、并串变换、
%每个符号采用1024个子载波，每个载波传输6个符号，添加的保护间隔为256个。
clc
clear all
close all
%参数设置
SubCarryN=128;                   %子载波数
SubCarry_eff=SubCarryN*1/4;      %有效子载波数
ratio=1/4;                %CP length ratio
fftLen=SubCarryN;         %FFT length
SymbN=100000;             %number of OFDM symbols
GuardLen=SubCarryN*ratio;    %CP length
SymLen=SubCarryN+GuardLen;  %length of each symbol after cp
beta=1/32; %rolloff ratio
SNR_dB=0:10;             %信噪比取值，以db为单位
% SNR_dB=5;%设置信噪比db
SNR=10*log10(SNR_dB);
%%----------------------------------------发射端-----------------------------------------------------%%
%输入比特序列长度=子载波数×符号数×每符号比特数
% for index=1:length(SNR)
% for p=1:200
SignalLen=SubCarry_eff*SymbN;
M=4;                                  %调制阶数
Signal=randi([0 M-1],1,SignalLen);
if M==4
    Signal=pskmod(Signal,M,pi/4);%
else
    Signal=qammod(Signal,M);
end
ParaSymSig=reshape(Signal,SubCarry_eff,SymbN);
% plot(ParaSymSig,'*r');
%在中间补零
%x=[ParaSymSig(1:SubCarryN*3/8,1:SymbN);zeros(SubCarryN/4,SymbN);ParaSymSig(SubCarryN*3/8+1:SubCarry_eff,1:SymbN)];%有效子载波占3/4
x=[ParaSymSig(1:SubCarryN*1/8,1:SymbN);zeros(SubCarryN*3/4,SymbN);ParaSymSig(SubCarryN*1/8+1:SubCarry_eff,1:SymbN)];%有效子载波占1/4
y=ifft(x,fftLen,1);   %通过傅里叶反变换，将频域数据转化为时域数据，1024行，500列

%在数组中加循环前缀的方法
y_CP=[y(fftLen-GuardLen+1:fftLen,:);y]; %加入循环前缀

%并串变换
TrData=reshape(y_CP,1,SymLen*SymbN);%并串变换，转换成一串，变成一行，7680列

for N=1:length(SNR_dB)
%% The channel varies among different symbols
% for n=1:SymbN   
    ReData=awgn(TrData,SNR_dB(N),'measured');%awgn信道
    r1=reshape(ReData,SymLen,SymbN);           %接收数据（可以选取一部分）
%     去循环前缀
    rev_signal=r1(GuardLen+1:SymLen,:);%去掉循环前缀并恢复格式
    rece_signal=fft(rev_signal);
% end
%尝试多个符号累加
%rece_signal=reshape(rece_signal,1024,SymbN);
F1=zeros(1,SubCarryN);%F1是识别空载波的特征量
% F1=zeros(1,SubCarryN);%F1是识别空载波的特征量
%识别空载波
F1countLen=SymbN;
for i=1:SubCarryN
    [cumulants11]=func_get_cumulants(rece_signal(i,:));
    F1(i)=cumulants11(7);%=C21, the power of each subcarrier
%     F1(i)=cumulants11(8);%=C42/((C21)^2)
end
%画出F1的分布
figure;
i=1:SubCarryN;
plot(i,F1);hold on;

%%

% % 通过识别出的空载波求N
% %矩阵，用来存F1的数值
F1_B=0.1*max(F1);%门限
[F1v,F1index]=find(F1<F1_B);
EmptySub_Num=length(F1index);  %Number of empty subcarries
NoEmptySub_Num=SubCarryN-EmptySub_Num; %Number of non-empty subcarries

EmptySub=rece_signal(F1index,:); %store received empty subcarries
NoEmptySub=rece_signal(setdiff([1:SubCarryN],F1index),:);%store received non-empty subcarries

% NoisePower=mean(mean(abs(EmptySub.^2)));%power of noise
% Sum_Power=mean(mean(abs(rece_signal.^2)));%power of sum
% Sub_Power=mean(mean(abs(NoEmptySub.^2)));%power of sub
% 
% snr1=(Sum_Power-(NoisePower))/(NoisePower);
% SNR_est=10*log10(snr1)

noise_part=sum(sum(abs(EmptySub).^2));%energy of noise
all=sum(sum(abs(rece_signal).^2));%power of sum
signal_part=sum(sum(abs(NoEmptySub).^2));%power of sub

noise = noise_part * (SubCarryN/(SubCarryN - NoEmptySub_Num));
signal = signal_part - noise_part * (NoEmptySub_Num/(SubCarryN - NoEmptySub_Num));

snr=signal/noise;
SNR_est=10*log10(snr)
end




