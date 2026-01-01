%% OFDM和单载波信号识别程序
%使用峰均比测试OFDM和单载波信号
clc;
clear;
close all;
% rng default;

%% 信号参数配置
M = [4 16 32 64];                                   % 调制阶数
M_APSK = [4 12 16];                                 % APSK星座点数
acsii = [1 2.84 5.54];                              % APSK星座半径
K = log2(M);                                        % 每符号包含bit数
NS = 600;                                          % 符号数
sps = 8;                                            % 每符号采样数
EbN0_dB = 10;                                       % Eb/N0(dB)
code_method = 'gray';                               % 编码方式，此处采用格雷码，非格雷码可将此处改为'bin'
%% OFDM信号参数设置
SubCarryN=128;                                      %子载波数总数
SubCarry_eff=SubCarryN*1/4;                         %有效子载波数
ratio=1/4;                                          %CP length ratio，循环前缀比例
fftLen=SubCarryN;                                   %FFT length，FFT长度
SymbN=1000;                                        %number of OFDM symbols
GuardLen=SubCarryN*ratio;                           %CP length
SymLen=SubCarryN+GuardLen;                          %length of each symbol after cp，符号总长度
beta=1/32;                                          %rolloff ratio
fc = 20e3;
sps = 8;                                            % 每符号采样数
f  = 1e6;                                           % 信号频率
fs = sps*f;                                         % 采样频率
fc = 1.375e9;                                       % 载波频率
ppm = 2.008e-2;
fm = f*ppm;

%% 信道生成
%% 参数设置
f_max=4e-7*fc;              %最大多普勒频偏，按移速120m/s算的
% f_max=60;              %最大多普勒频偏，按移速120m/s算的
% f_max=100000;
Flat = doppler('Flat');     %多普勒频谱结构，下同
Jakes=doppler('Jakes');
GaussI=doppler('BiGaussian','PowerGains',[5 1],'NormalizedCenterFrequencies',[-0.8 0.4],'NormalizedStandardDeviations',[0.05 0.1]);%这个不用管了，我算的参数
GaussII=doppler('BiGaussian','PowerGains',[20*sqrt(10)/3 1],'NormalizedCenterFrequencies',[0.7 -0.4],'NormalizedStandardDeviations',[0.1 0.15]);
%% 瑞利衰落信道建立

%% 单径信道+多普勒
rayleighChan = comm.RayleighChannel( ...
    'SampleRate',fs, ...
    'NormalizePathGains',true, ...
    'MaximumDopplerShift',f_max, ...
    'RandomStream','mt19937ar with seed', ...
    'Visualization','Off',...
    'Seed',22, ...
    'PathGainsOutputPort',true);


%% 接收端参数
SNR =0:30;                        % 信噪比区间
TT = 50;                               % 独立实验次数

%% 接收端仿真
for index = 1:length(SNR)
    pp=0;
    for tt = 1:TT
        %% 通信链路
        dataIn2 = randi([0 M(2)-1],NS,1);            % 生成随机符号数据
        dataIn1 = randi([0 M(1)-1],NS,1);            % 生成随机符号数据
        dataIn3 = randi([0 M(3)-1],NS,1);            % 生成随机符号数据
        dataIn4 = randi([0 M(4)-1],NS,1);            % 生成随机符号数据
        %% 调制
        dataMod_QPSK = pskmod(dataIn1,M(1),pi/M(1));                % QPSK调制
        dataMod_16QAM = qammod(dataIn2,M(2));                     % QAM调制
        dataMod_64QAM = qammod(dataIn4,64);                     % 64QAM调制
        dataMod_16APSK = apskmod(dataIn2,M_APSK,acsii);         % 16APSK调制
        dataMod_32APSK = apskmod(dataIn3,M_APSK,acsii);         % 32APSK调制
        %         scatterplot(dataMod);
        %ofdm
        SignalLen=SubCarry_eff*SymbN;
        Signal=randi([0 M(1)-1],1,SignalLen);
        if M(1)==4
            Signal=pskmod(Signal,M(1),pi/4);%
        else
            Signal=qammod(Signal,M(1));
        end
        ParaSymSig=reshape(Signal,SubCarry_eff,SymbN);
        % plot(ParaSymSig,'*r');
        
        %在中间补零
        x=[ParaSymSig(1:SubCarryN*1/8,1:SymbN);zeros(SubCarryN*3/4,SymbN);ParaSymSig(SubCarryN*1/8+1:SubCarry_eff,1:SymbN)];%有效子载波占1/4
        y=ifft(x,fftLen,1);                                 %通过傅里叶反变换，将频域数据转化为时域数据，1024行，500列
        %在数组中加循环前缀的方法
        y_CP=[y(fftLen-GuardLen+1:fftLen,:);y];             %加入循环前缀
        %rcosfilter window
        window_coff=rcoswindow(beta,SymLen);
        y_window=zeros(SymLen,SymbN);
        for n=1:SymbN
            y_window(:,n)=y_CP(:,n).*window_coff.';
        end
        %并串变换
        TrData=reshape(y_window,1,SymLen*SymbN);%并串变换，转换成一串，变成一行，7680列
        
        %% 成形滤波
        span = 8;                                           % 截断符号数
        alpha = 0.35;
        h = rcosdesign(alpha,span,sps,'sqrt');              % 平方根升余弦滤波器
        signalSend_QPSK = upfirdn(dataMod_QPSK,h,sps);
        signalSend_16QAM = upfirdn(dataMod_16QAM,h,sps);
        signalSend_64QAM = upfirdn(dataMod_64QAM,h,sps);
        signalSend_16APSK = upfirdn(dataMod_16APSK,h,sps);
        signalSend_32APSK = upfirdn(dataMod_32APSK,h,sps);
        signalMod_OFDM=TrData';
        
        signalMod_QPSK=signalSend_QPSK;
        signalMod_16QAM=signalSend_16QAM;
        signalMod_64QAM=signalSend_64QAM;
        signalMod_16APSK=signalSend_16APSK;
        signalMod_32APSK=signalSend_32APSK;
        carrier_OFDM = exp(1j*2*pi*(fc+fm)/fs.*(1:length(signalMod_OFDM)));
    signalMod_OFDM = signalMod_OFDM.*carrier_OFDM';
        %%
        signalAWGN_QPSK = awgn(signalMod_QPSK,SNR(index),'measured');
        signalRecv_QPSK = signalAWGN_QPSK;
        PAPR_QPSK(index) = func_get_PAPR(signalRecv_QPSK');
        
        signalAWGN_16QAM = awgn(signalMod_16QAM,SNR(index),'measured');
        signalRecv_16QAM = signalAWGN_16QAM;
        PAPR_16QAM(index) = func_get_PAPR(signalRecv_16QAM');
        
        signalAWGN_64QAM = awgn(signalMod_64QAM,SNR(index),'measured');
        signalRecv_64QAM = signalAWGN_64QAM;
        PAPR_64QAM(index)= func_get_PAPR(signalRecv_64QAM');
        
        signalAWGN_16APSK = awgn(signalMod_16APSK,SNR(index),'measured');
        signalRecv_16APSK = signalAWGN_16APSK;
        PAPR_16APSK(index) = func_get_PAPR(signalRecv_16APSK');
        
        signalAWGN_32APSK = awgn(signalMod_32APSK,SNR(index),'measured');
        signalRecv_32APSK = signalAWGN_32APSK;
        PAPR_32APSK(index) = func_get_PAPR(signalRecv_32APSK');
        
        [signalMod_OFDM,pathr_OFDM]= rayleighChan(signalMod_OFDM);
        signalAWGN_OFDM = awgn(signalMod_OFDM,SNR(index),'measured');
        Wc=2*30000/fs;                                          %截止频率 3Hz
        [b,a]=butter(4,Wc,'low');  % 四阶的巴特沃斯低通滤波
        signalAWGN_OFDM=filter(b,a,signalAWGN_OFDM);
        signalRecv_OFDM = signalAWGN_OFDM';
        PAPR_OFDM(index)= func_get_PAPR(signalRecv_OFDM);
%         PAPR_OFDM_short(index)= func_get_PAPR_short(signalRecv_OFDM);
        for o=1:6
            PAPR_OFDM_short(o)= func_get_PAPR(signalRecv_OFDM(100*(o-1)+1:100*o));
        end
         PAPR_OFDM_short(index)=mean(PAPR_OFDM_short);
        if PAPR_OFDM(index)>9.5 & (PAPR_OFDM_short(index))<8
            ppp(tt)=1;
        else
            ppp(tt)=0;
        end
        pp=pp+ppp(tt);
    end
    acc(index)=pp/TT;
end
figure(1)
subplot(221)
plot(SNR,acc,'r-d')
xlabel('SNR(dB)','FontName','Times New Roman');
ylabel('Acc','FontName','Times New Roman');
legend('acc-OFDM','FontName','Times New Roman');
title("估计准确率")
subplot(222)
plot(SNR,PAPR_OFDM,'r-d')
xlabel('SNR(dB)','FontName','Times New Roman');
ylabel('PAPR','FontName','Times New Roman');
legend('PAPR-OFDM','FontName','Times New Roman');
title("PAPR-截短前")
subplot(223)
plot(SNR,PAPR_OFDM_short,'r-d')
xlabel('SNR(dB)','FontName','Times New Roman');
ylabel('PAPR_short','FontName','Times New Roman');
legend('PAPR_short-OFDM','FontName','Times New Roman');
title("PAPR-截短后")
subplot(224)
plot(SNR,PAPR_OFDM-PAPR_OFDM_short,'r-d')
xlabel('SNR(dB)','FontName','Times New Roman');
ylabel('delta_PAPR','FontName','Times New Roman');
legend('delta_PAPR-OFDM','FontName','Times New Roman');
title("截短前的PAPR和截短后的PAPR做差结果")


function PAPR_dB = func_get_PAPR(signalData)
% PAPR=abs(mean(max(signalData.^2)/mean(signalData.^2)));
% PAPR=max(abs(signalData).^2)/mean(abs(signalData).^2);
Nx=length(signalData);
xI=real(signalData);
xQ=imag(signalData);
Power=xI.*xI+xQ.*xQ;
PeakP=max(Power);
PeakP_dB=10*log10(PeakP);
AvgP=sum(Power)/Nx;
AvgP_dB=10*log10(AvgP);
PAPR_dB=10*log10(PeakP/AvgP);
end

function PAPR_dB_short = func_get_PAPR_short(signalData)
% PAPR=abs(mean(max(signalData.^2)/mean(signalData.^2)));
% PAPR=max(abs(signalData).^2)/mean(abs(signalData).^2);
signalData=signalData(1:100);
Nx=length(signalData);
xI=real(signalData);
xQ=imag(signalData);
Power=xI.*xI+xQ.*xQ;
PeakP=max(Power);
PeakP_dB=10*log10(PeakP);
AvgP=sum(Power)/Nx;
AvgP_dB=10*log10(AvgP);
PAPR_dB_short=10*log10(PeakP/AvgP);
end