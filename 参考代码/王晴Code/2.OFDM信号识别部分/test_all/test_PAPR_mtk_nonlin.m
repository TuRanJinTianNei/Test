%% OFDM和单载波信号识别程序
%使用峰均比测试OFDM和单载波信号
%非线性信道下的OFDM和各个单载波的PAPR分布
%蒙特卡洛仿真100次
%% 
clc;
clear;
close all;
% rng default;

%% 信号参数配置
M = [4 16 32 64];                                             % 调制阶数\

APSK_16 = [4 12];                                 % APSK星座点数
acsii_16 = [1 2];                              % APSK星座半径

APSK_32 = [4 8 20];                                 % APSK星座点数
acsii_32 = [0.3 0.7 1.2];                              % APSK星座半径

APSK_64 = [8 12 16 28];                                 % APSK星座点数
acsii_64 = [0.5 1 1.3 2];                              % APSK星座半径
K = log2(M);                                        % 每符号包含bit数
NS = 60;                                          % 符号数
sps = 8;                                            % 每符号采样数
EbN0_dB = 10;                                       % Eb/N0(dB)
code_method = 'gray';                               % 编码方式，此处采用格雷码，非格雷码可将此处改为'bin'
%% OFDM信号参数设置
SubCarryN=128;                                      %子载波数总数
SubCarry_eff=SubCarryN*1;                         %有效子载波数
ratio=1/4;                                          %CP length ratio，循环前缀比例
fftLen=SubCarryN;                                   %FFT length，FFT长度
SymbN=100;                                        %number of OFDM symbols
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

tRate = sps/2;                                      % 步进参数
f  = 1e6;                                           % 信号频率
fs = sps*f;                                         % 采样频率
fc = 1.375e9;                                       % 载波频率
ppm = 2.008e-2;
fm = f*ppm;
%% 信道生成
%% 参数设置
f_max=4e-7*fs;              %最大多普勒频偏，按移速120m/s算的
% f_max=1000;
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
TT = 10;                               % 独立实验次数
%% 接收端仿真
for index = 1:length(SNR)
        for tt = 1:TT
    %% 通信链路
    dataIn2 = randi([0 M(2)-1],NS,1);            % 生成随机符号数据
    dataIn1 = randi([0 M(1)-1],NS,1);            % 生成随机符号数据
    dataIn3 = randi([0 M(3)-1],NS,1);            % 生成随机符号数据
    dataIn4 = randi([0 M(4)-1],NS,1);            % 生成随机符号数据
    %% 调制
    dataMod_QPSK = pskmod(dataIn1,M(1),pi/M(1));                % QPSK调制
    dataMod_16QAM = qammod(dataIn2,M(2));                     % QAM调制
    dataMod_32QAM = qammod(dataIn3,M(3));                     % 32QAM调制
    dataMod_64QAM = qammod(dataIn4,M(4));                     % 64QAM调制
    dataMod_16APSK = apskmod(dataIn2,APSK_16,acsii_16);         % 16APSK调制
    dataMod_32APSK = apskmod(dataIn3,APSK_32,acsii_32);         % 32APSK调制
    dataMod_64APSK = apskmod(dataIn4,APSK_64,acsii_64);         % 64APSK调制
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
    
    %%
%     PP_QPSK = sum(abs(dataMod_QPSK).^2)/NS;
%     dataMod_QPSK = dataMod_QPSK./sqrt(PP_QPSK);
%     
%     PP_16QAM = sum(abs(dataMod_16QAM).^2)/NS;
%     dataMod_16QAM = dataMod_16QAM./sqrt(PP_16QAM);
%     
%     PP_64QAM = sum(abs(dataMod_64QAM).^2)/NS;
%     dataMod_64QAM = dataMod_64QAM./sqrt(PP_64QAM);
%     
%     PP_16APSK = sum(abs(dataMod_16APSK).^2)/NS;
%     dataMod_16APSK = dataMod_16APSK./sqrt(PP_16APSK);
%     
%     PP_32APSK = sum(abs(dataMod_32APSK).^2)/NS;
%     dataMod_32APSK = dataMod_32APSK./sqrt(PP_32APSK);
    %% 成形滤波
    span = 8;                                           % 截断符号数
    alpha = 0.35;
    h = rcosdesign(alpha,span,sps,'sqrt');              % 平方根升余弦滤波器
signalSend_QPSK = upfirdn(dataMod_QPSK,h,sps);
    signalSend_16QAM = upfirdn(dataMod_16QAM,h,sps);
    signalSend_32QAM = upfirdn(dataMod_32QAM,h,sps);
    signalSend_64QAM = upfirdn(dataMod_64QAM,h,sps);
    signalSend_16APSK = upfirdn(dataMod_16APSK,h,sps);
    signalSend_32APSK = upfirdn(dataMod_32APSK,h,sps);
    signalSend_64APSK = upfirdn(dataMod_64APSK,h,sps);
    %     signalSend_OFDM = upfirdn(TrData,h,sps);
    signalMod_OFDM=TrData';
    
    signalMod_QPSK=signalSend_QPSK;
    signalMod_16QAM=signalSend_16QAM;
    signalMod_32QAM=signalSend_32QAM;
    signalMod_64QAM=signalSend_64QAM;
    signalMod_16APSK=signalSend_16APSK;
    signalMod_32APSK=signalSend_32APSK;
    signalMod_64APSK=signalSend_64APSK;

    signalChannel_QPSK = Channel_Nonlinear_Volterra(signalMod_QPSK);
    signalAWGN_QPSK = awgn(signalChannel_QPSK,SNR(index),'measured');
    %signalAWGN_QPSK=rayleighChan(signalAWGN_QPSK');
%     Wc=2*30000/fs;                                          %截止频率 3Hz
%     [b,a]=butter(4,Wc,'low');  % 四阶的巴特沃斯低通滤波
%     signalAWGN_QPSK=filter(b,a,signalAWGN_QPSK);
    signalRecv_QPSK = signalAWGN_QPSK;
    PAPR_QPSK(tt)= func_get_PAPR(signalRecv_QPSK);
    
    signalChannel_16QAM = Channel_Nonlinear_Volterra(signalMod_16QAM);
    signalAWGN_16QAM = awgn(signalChannel_16QAM,SNR(index),'measured');
    %signalAWGN_16QAM=rayleighChan(signalAWGN_16QAM');
%              Wc=2*30000/fs;                                          %截止频率 3Hz
%     [b,a]=butter(4,Wc,'low');  % 四阶的巴特沃斯低通滤波
%     signalAWGN_16QAM=filter(b,a,signalAWGN_16QAM);
    signalRecv_16QAM = signalAWGN_16QAM;
    PAPR_16QAM(tt)= func_get_PAPR(signalRecv_16QAM);
    
    signalChannel_32QAM = Channel_Nonlinear_Volterra(signalMod_32QAM);
    signalAWGN_32QAM = awgn(signalChannel_32QAM,SNR(index),'measured');
    signalRecv_32QAM = signalAWGN_32QAM;
    PAPR_32QAM(tt) = func_get_PAPR(signalRecv_32QAM);
    
    
    signalChannel_64QAM = Channel_Nonlinear_Volterra(signalMod_64QAM);
    signalAWGN_64QAM = awgn(signalChannel_64QAM,SNR(index),'measured');
    %signalAWGN_64QAM=rayleighChan(signalAWGN_64QAM');
%      Wc=2*30000/fs;                                          %截止频率 3Hz
%     [b,a]=butter(4,Wc,'low');  % 四阶的巴特沃斯低通滤波
%     signalAWGN_64QAM=filter(b,a,signalAWGN_64QAM);
    signalRecv_64QAM = signalAWGN_64QAM;
    PAPR_64QAM(tt) = func_get_PAPR(signalRecv_64QAM);

    
    signalChannel_16APSK = Channel_Nonlinear_Volterra(signalMod_16APSK);
    signalAWGN_16APSK = awgn(signalChannel_16APSK,SNR(index),'measured');
    %signalAWGN_16APSK=rayleighChan(signalAWGN_16APSK');
%      Wc=2*30000/fs;                                          %截止频率 3Hz
%     [b,a]=butter(4,Wc,'low');  % 四阶的巴特沃斯低通滤波
%     signalAWGN_16APSK=filter(b,a,signalAWGN_16APSK);
    signalRecv_16APSK = signalAWGN_16APSK;
    PAPR_16APSK(tt)= func_get_PAPR(signalRecv_16APSK);
    
    signalChannel_32APSK = Channel_Nonlinear_Volterra(signalMod_32APSK);
    signalAWGN_32APSK = awgn(signalChannel_32APSK,SNR(index),'measured');
    %signalAWGN_32APSK=rayleighChan(signalAWGN_32APSK');
%         Wc=2*30000/fs;                                          %截止频率 3Hz
%     [b,a]=butter(4,Wc,'low');  % 四阶的巴特沃斯低通滤波
%     signalAWGN_32APSK=filter(b,a,signalAWGN_32APSK);
    signalRecv_32APSK = signalAWGN_32APSK;
    PAPR_32APSK(tt) = func_get_PAPR(signalRecv_32APSK);
    
    signalChannel_64APSK = Channel_Nonlinear_Volterra(signalMod_64APSK);
    signalAWGN_64APSK = awgn(signalChannel_64APSK,SNR(index),'measured');
    signalRecv_64APSK = signalAWGN_64APSK;
    PAPR_64APSK(tt) = func_get_PAPR(signalRecv_64APSK);
    
    signalChannel_OFDM = Channel_Nonlinear_Volterra(signalMod_OFDM);
    signalAWGN_OFDM = awgn(signalChannel_OFDM,SNR(index),'measured');
    %signalAWGN_OFDM=rayleighChan(signalAWGN_OFDM');
%     Wc=2*30000/fs;                                          %截止频率 3Hz
%     [b,a]=butter(4,Wc,'low');  % 四阶的巴特沃斯低通滤波
%     signalAWGN_OFDM=filter(b,a,signalAWGN_OFDM);
    signalRecv_OFDM = signalAWGN_OFDM';
    PAPR_OFDM(tt)= func_get_PAPR(signalRecv_OFDM);
        end
    PAPR_QPSK_la(index)=mean(PAPR_QPSK);
    PAPR_16QAM_la(index)=mean(PAPR_16QAM);
    PAPR_32QAM_la(index)=mean(PAPR_32QAM);
    PAPR_64QAM_la(index)=mean(PAPR_64QAM);
    PAPR_16APSK_la(index)=mean(PAPR_16APSK);
    PAPR_32APSK_la(index)=mean(PAPR_32APSK);
    PAPR_64APSK_la(index)=mean(PAPR_64APSK);
    PAPR_OFDM_la(index)=mean(PAPR_OFDM);
end
plot(SNR,PAPR_QPSK_la,'r-d')
hold on
plot(SNR,PAPR_16QAM_la,'b-s')
hold on
plot(SNR,PAPR_32QAM_la,'c-x')
hold on
plot(SNR,PAPR_64QAM_la,'k-p')
hold on
plot(SNR,PAPR_16APSK_la,'g-h')
hold on
plot(SNR,PAPR_32APSK_la,'k-+','LineWidth', 2)
hold on
plot(SNR,PAPR_64APSK_la,'c-v','LineWidth', 2)
hold on
plot(SNR,PAPR_OFDM_la,'m-^')
xlabel('SNR(dB)','FontName','Times New Roman');
ylabel('PAPR','FontName','Times New Roman');
legend('PAPR-QPSK','PAPR-16QAM','PAPR-32QAM','PAPR-64QAM','PAPR-16APSK','PAPR-32APSK','PAPR-64APSK','PAPR-OFDM','FontName','Times New Roman');

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