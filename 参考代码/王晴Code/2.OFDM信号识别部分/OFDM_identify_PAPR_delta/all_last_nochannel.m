%% OFDM和单载波信号识别程序
%使用峰均比测试OFDM和单载波信号
%无信道下的OFDM和各个单载波的PAPR分布
%蒙特卡洛仿真100次
%%
clc;
clear;
close all;
% rng default;

%% 信号参数配置
M = [4 16 32 64];                                             % 调制阶数
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
        dataMod_16QAM  = qammod(dataIn2,M(2));                     % QAM调制
        dataMod_64QAM  = qammod(dataIn4,64);                     % 64QAM调制
        dataMod_16APSK  = apskmod(dataIn2,M_APSK,acsii);         % 16APSK调制
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
        %%
        signalAWGN_QPSK = awgn(signalMod_QPSK,SNR(index),'measured');
        fs=1024;%采样率
        f=linspace(-fs/2,fs/2,length(signalAWGN_QPSK));%频域横坐标，注意奈奎斯特采样定理，最大原信号最大频率不超过采样频率的一半
        f=f/(2*max(f));
        X = fftshift(fft(signalAWGN_QPSK)); %用fft得出离散傅里叶变换
        Q=20*log10(abs(X));
        X=X(1:length(X)*0.25);
        signalAWGN_QPSK=reshape(signalAWGN_QPSK,1,length(signalAWGN_QPSK));
        PAPR_QPSK(tt)= func_get_PAPR(signalAWGN_QPSK);
        X=reshape(X,1,length(X));
        X=[X,X,X,X,X];
        signalAWGN_QPSK=ifft(X);
        PAPR_QPSK_long(tt)= func_get_PAPR(signalAWGN_QPSK);

        PAPR_QPSK(tt) =PAPR_QPSK_long(tt)-PAPR_QPSK(tt);
        
        signalAWGN_16QAM = awgn(signalMod_16QAM,SNR(index),'measured');
        fs=1024;%采样率
        f=linspace(-fs/2,fs/2,length(signalAWGN_16QAM));%频域横坐标，注意奈奎斯特采样定理，最大原信号最大频率不超过采样频率的一半
        f=f/(2*max(f));
        X = fftshift(fft(signalAWGN_16QAM)); %用fft得出离散傅里叶变换
        Q=20*log10(abs(X));
        X=X(1:length(X)*0.25);
        signalAWGN_16QAM=reshape(signalAWGN_16QAM,1,length(signalAWGN_16QAM));
        PAPR_16QAM(tt)= func_get_PAPR(signalAWGN_16QAM);
        X=reshape(X,1,length(X));
        X=[X,X,X,X,X];
        signalAWGN_16QAM=ifft(X);
        PAPR_16QAM_long(tt)= func_get_PAPR(signalAWGN_16QAM);

        PAPR_16QAM(tt) = PAPR_16QAM_long(tt)-PAPR_16QAM(tt);
        
        signalAWGN_64QAM = awgn(signalMod_64QAM,SNR(index),'measured');
        fs=1024;%采样率
        f=linspace(-fs/2,fs/2,length(signalAWGN_64QAM));%频域横坐标，注意奈奎斯特采样定理，最大原信号最大频率不超过采样频率的一半
        f=f/(2*max(f));
        X = fftshift(fft(signalAWGN_64QAM)); %用fft得出离散傅里叶变换
        Q=20*log10(abs(X));
        X=X(1:length(X)*0.25);
        signalAWGN_64QAM=reshape(signalAWGN_64QAM,1,length(signalAWGN_64QAM));
        PAPR_64QAM(tt)= func_get_PAPR(signalAWGN_64QAM);
        X=reshape(X,1,length(X));
        X=[X,X,X,X,X];
        signalAWGN_64QAM=ifft(X);
        PAPR_64QAM_long(tt)= func_get_PAPR(signalAWGN_64QAM);
        PAPR_64QAM(tt)= PAPR_64QAM_long(tt)-PAPR_64QAM(tt);
        
        signalAWGN_16APSK = awgn(signalMod_16APSK,SNR(index),'measured');
        fs=1024;%采样率
        f=linspace(-fs/2,fs/2,length(signalAWGN_16APSK));%频域横坐标，注意奈奎斯特采样定理，最大原信号最大频率不超过采样频率的一半
        f=f/(2*max(f));
        X = fftshift(fft(signalAWGN_16APSK)); %用fft得出离散傅里叶变换
        Q=20*log10(abs(X));
        X=X(1:length(X)*0.25);
        signalAWGN_16APSK=reshape(signalAWGN_16APSK,1,length(signalAWGN_16APSK));
        PAPR_16APSK(tt)= func_get_PAPR(signalAWGN_16APSK);
        X=reshape(X,1,length(X));
        X=[X,X,X,X,X];
        signalAWGN_16APSK=ifft(X);
        PAPR_16APSK_long(tt)= func_get_PAPR(signalAWGN_16APSK);
        PAPR_16APSK(tt) = PAPR_16APSK_long(tt)-PAPR_16APSK(tt);
        
        signalAWGN_32APSK = awgn(signalMod_32APSK,SNR(index),'measured');
        fs=1024;%采样率
        f=linspace(-fs/2,fs/2,length(signalAWGN_32APSK));%频域横坐标，注意奈奎斯特采样定理，最大原信号最大频率不超过采样频率的一半
        f=f/(2*max(f));
        X = fftshift(fft(signalAWGN_32APSK)); %用fft得出离散傅里叶变换
        Q=20*log10(abs(X));
        X=X(1:length(X)*0.25);
        signalAWGN_32APSK=reshape(signalAWGN_32APSK,1,length(signalAWGN_32APSK));
        PAPR_32APSK(tt)= func_get_PAPR(signalAWGN_32APSK);
        X=reshape(X,1,length(X));
        X=[X,X,X,X,X];
        signalAWGN_32APSK=ifft(X);
        PAPR_32APSK_long(tt)= func_get_PAPR(signalAWGN_32APSK);
        
        PAPR_32APSK(tt) = PAPR_32APSK_long(tt)-PAPR_32APSK(tt);
        
        signalAWGN_OFDM = awgn(signalMod_OFDM,SNR(index),'measured');
        fs=1024;%采样率
        f=linspace(-fs/2,fs/2,length(signalAWGN_OFDM));%频域横坐标，注意奈奎斯特采样定理，最大原信号最大频率不超过采样频率的一半
        f=f/(2*max(f));
        X = fftshift(fft(signalAWGN_OFDM)); %用fft得出离散傅里叶变换
        Q=20*log10(abs(X));
        X=X(1:length(X)*0.25);
        signalAWGN_OFDM=reshape(signalAWGN_OFDM,1,length(signalAWGN_OFDM));
        PAPR_OFDM(tt)= func_get_PAPR(signalAWGN_OFDM);
        X=reshape(X,1,length(X));
        X=[X,X,X,X,X];
        signalAWGN_OFDM=ifft(X);
        PAPR_OFDM_long(tt)= func_get_PAPR(signalAWGN_OFDM);
        PAPR_OFDM(tt)= PAPR_OFDM_long(tt)-PAPR_OFDM(tt);
    end
    PAPR_QPSK_la(index)=mean(PAPR_QPSK);
    PAPR_16QAM_la(index)=mean(PAPR_16QAM);
    PAPR_64QAM_la(index)=mean(PAPR_64QAM);
    PAPR_16APSK_la(index)=mean(PAPR_16APSK);
    PAPR_32APSK_la(index)=mean(PAPR_32APSK);
    PAPR_OFDM_la(index)=mean(PAPR_OFDM);
end
plot(SNR,PAPR_QPSK_la,'r-d')
hold on
plot(SNR,PAPR_16QAM_la,'b-s')
hold on
plot(SNR,PAPR_64QAM_la,'k-p')
hold on
plot(SNR,PAPR_16APSK_la,'g-h')
hold on
plot(SNR,PAPR_32APSK_la,'c-+','LineWidth', 2)
hold on
plot(SNR,PAPR_OFDM_la,'m-^')
xlabel('SNR(dB)','FontName','Times New Roman');
ylabel('delta_PAPR','FontName','Times New Roman');
legend('delta_PAPR-QPSK','delta_PAPR-16QAM','delta_PAPR-64QAM','delta_PAPR-16APSK','delta_PAPR-32APSK','delta_PAPR-OFDM','FontName','Times New Roman');

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
