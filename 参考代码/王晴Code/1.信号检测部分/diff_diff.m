%两个信号重叠，计算带宽的相关程序
%

%-----------------------------------方法说明-------------------------------------------------------
%两个重叠信号的功率谱图有四个“边界”从左到右依次记为1，2，3，4号边界，本程序旨在准确求出四个边界，之后再进行讨论是完全重叠还是部分重叠
%首先利用一次拐点法，求出1号边界和4号边界
%用1号边界和4号边界截取数据，计算上凸起部分的高度，记为H
%求出最大值记为max，该最大值必定在上凸起的通带内
%设置门限1=max-H*1/3,门限2=max-H-3
%利用两个门限分别计算2和3，1和4号边界
%--------------------------------------------------------------------------------------------------
clc
clear all
close all
%% 信号1
SubCarryNN=256;      %有效子载波数%%四倍过采样
SubCarryN=1024;      %子载波数
ratio=1/4;           %循环前缀比例
fftLen=1024;         %FFT长度为1024
SymbN=20;            %一帧中OFDM符号个数、每个符号传输的子载波数、6×128=768
GuardLen=SubCarryN*ratio;    %保护时隙的长度
SNR=10;             %信噪比取值，以db为单位
%输入比特序列长度=子载波数×每载波符号数×每符号比特数10240
SignalLen=SubCarryNN*SymbN*2;           %比特序列长度，256X20X2
Signal=round(rand(1,SignalLen));       %输出带调制的二进制比特流,0,1、round是为了产生0，1
ParaBitSig=[];
%length(ParaBitSig)
%这是一个把一行数变成X行Y列的数组的方法
% for i=1:SubCarryNN
%     for j=1:SymbN*2
%         ParaBitSig(i,j)=Signal(i*j);   %串/并转换为行数SubCarryNN（有效子载波），列数为2*SymbN（每载波的比特数）
%     end
% end
ParaBitSig=reshape(Signal,SubCarryNN,SymbN*2);
%进行QPSK数据调制，将数据分为两个通道
%这是一个QPSK的调制方法
for j=1:SymbN
    ich(:,j)=ParaBitSig(:,2*j-1); %同相分量，ParaBitSig奇数列
    qch(:,j)=ParaBitSig(:,2*j);   %正交分量，ParaBitSig偶数列
end
kmod=1./sqrt(2);   %做一个根号二
ich0=ich.*2-1;     %这俩变成正负一的形式，1→1，0→-1
qch0=qch.*2-1;
ich1=ich0.*kmod;   %根号二
qch1=qch0.*kmod;
x=ich1+qch1.*sqrt(-1);      %产生复信号,768行，6列
% x=[zeros(32,6);x;zeros(32,6)];%在中间补零，补了256行
x=[x(1:128,1:20);zeros(768,20);x(129:256,1:20)];
% x=[x(1:384,1:20);zeros(256,20);x(385:768,1:20)];
y=ifft(x);   %通过傅里叶反变换，将频域数据转化为时域数据，1024行，6列
ich2=real(y);   %I信道取变换后的实部
qch2=imag(y);   %Q信道取变换后的虚部
%插入保护间隔
%在数组中加循环前缀的方法
ich3=[ich2(fftLen-GuardLen+1:fftLen,:);ich2];%行数从1024-256+1到1024,在前面加了256个，总共1024+256=1280个,1280×6个
qch3=[qch2(fftLen-GuardLen+1:fftLen,:);qch2];
%并串变换
ich4=reshape(ich3,1,(fftLen+GuardLen)*SymbN);%并串变换，转换成一串，变成一行，7680列
qch4=reshape(qch3,1,(fftLen+GuardLen)*SymbN);
TrData=ich4+qch4.*sqrt(-1);  %形成复数发射数据（加完保护时隙，保护时隙也进行相应操作）
ReData=awgn(TrData,SNR,'measured');%信道，加入高斯白噪声  %每一列是一个ofdm符号

%% 信号2
SubCarryNN2=768;      %有效子载波数%%四倍过采样
SubCarryN2=1024;      %子载波数
SymbN2=20;            %一帧中OFDM符号个数、每个符号传输的子载波数、6×128=768
%%----------------------------------------发射端-----------------------------------------------------%%
%输入比特序列长度=子载波数×每载波符号数×每符号比特数10240
SignalLen2=SubCarryNN2*SymbN2*2;           %比特序列长度，256X20X2
Signal2=round(rand(1,SignalLen2));       %输出带调制的二进制比特流,0,1、round是为了产生0，1
ParaBitSig2=[];
%length(ParaBitSig)
%这是一个把一行数变成X行Y列的数组的方法
% for i=1:SubCarryNN
%     for j=1:SymbN*2
%         ParaBitSig(i,j)=Signal(i*j);   %串/并转换为行数SubCarryNN（有效子载波），列数为2*SymbN（每载波的比特数）
%     end
% end
ParaBitSig2=reshape(Signal2,SubCarryNN2,SymbN2*2);
%进行QPSK数据调制，将数据分为两个通道
%这是一个QPSK的调制方法
for j=1:SymbN2
    ich22(:,j)=ParaBitSig2(:,2*j-1); %同相分量，ParaBitSig奇数列
    qch22(:,j)=ParaBitSig2(:,2*j);   %正交分量，ParaBitSig偶数列
end
kmod=1./sqrt(2);   %做一个根号二
ich02=ich22.*2-1;     %这俩变成正负一的形式，1→1，0→-1
qch02=qch22.*2-1;
ich12=ich02.*kmod;   %根号二
qch12=qch02.*kmod;
x2=ich12+qch12.*sqrt(-1);      %产生复信号,768行，6列

% x=[zeros(32,6);x;zeros(32,6)];%在中间补零，补了256行
% x2=[x(1:128,1:20);zeros(768,20);x(129:256,1:20)];
x2=[x2(1:384,1:20);zeros(256,20);x2(385:768,1:20)];
y2=ifft(x2);   %通过傅里叶反变换，将频域数据转化为时域数据，1024行，6列
ich22=real(y2);   %I信道取变换后的实部
qch22=imag(y2);   %Q信道取变换后的虚部

%插入保护间隔
%在数组中加循环前缀的方法
ich32=[ich22(fftLen-GuardLen+1:fftLen,:);ich22];%行数从1024-256+1到1024,在前面加了256个，总共1024+256=1280个,1280×6个
qch32=[qch22(fftLen-GuardLen+1:fftLen,:);qch22];

%并串变换
ich42=reshape(ich32,1,(fftLen+GuardLen)*SymbN2);%并串变换，转换成一串，变成一行，7680列
qch42=reshape(qch32,1,(fftLen+GuardLen)*SymbN2);
TrData42=ich42+qch42.*sqrt(-1);  %形成复数发射数据（加完保护时隙，保护时隙也进行相应操作）
ReData42=awgn(TrData42,SNR,'measured');%信道，加入高斯白噪声  %每一列是一个ofdm符号
ReData_same =ReData42.*exp(1j * (0.1) * pi * 2 .* (1:length(ReData42)));
ReData=ReData+ReData_same;

% ReData=randn(1,length(ReData));
%ReData=awgn(ReData,SNR,'measured');%信道，加入高斯白噪声  %每一列是一个ofdm符号
B_carrier=200e3;  % 子载波间隔
fs=B_carrier*SubCarryN*4;         % 采样率
fft_ofdm=fftshift(fft(real(ReData).^2));
log_ofdm=20*log10(abs(fft_ofdm));
NN=length(ReData42);           % 采样点个数
n=0:NN-1;
freq=(n/NN-0.5)*fs/10e5;  %横坐标
figure(1)
plot(freq,log_ofdm);
xlabel('f/MHz');
ylabel('幅度');

%% 两个信号叠加
%画出功率谱
fs=1024;%采样率
f=linspace(-fs/2,fs/2,25600);%频域横坐标，注意奈奎斯特采样定理，最大原信号最大频率不超过采样频率的一半
f=f/(2*max(f));
Fs=1;
nfft=1024;
window=hanning(150);
noverlap=50;%每个窗口之间重叠的长度，通常取33%~50%。窗口之间重叠得越多，图像越平滑（blurred）；反之则更粗糙（blocky）
[Pxx1,f]=pwelch(ReData,window,noverlap,length(ReData),Fs);%功率谱密度
Pxx=10*log10(fftshift(Pxx1));
figure(2)
plot(f,Pxx)
title("Welch法功率谱估计")

%% 用一次拐点法求出1号边界和4号边界
Pxx=Pxx+abs(min(Pxx));%进行平移
Pxx_int=zeros(size(Pxx,1),1);
for x=1:1:size(Pxx,1)
    if x==1
    Pxx_int(x)=(1/(size(Pxx,1)-1))*Pxx(x+1);
%     Pxx_int(x)=Pxx(x+1);
    else
    Pxx_int(x)=Pxx_int(x-1)+(1/(size(Pxx,1)-1))*Pxx(x);
%     Pxx_int(x)=Pxx_int(x-1)+Pxx(x);
    end
end
figure(3)
plot(linspace(0,1,length(Pxx_int)),Pxx_int);
axis([0 1 0 max(Pxx_int)+1]);
title("Welch法功率分布函数")
y_=Pxx_int;
x_=(linspace(0,1,length(Pxx_int)))';
loc=[x_,y_];%坐标（第一列是横坐标，第二列是纵坐标）
%曲线上选一点K，作垂线KA垂直于x轴，K,A,O构成直角三角形
m1=0;%循环变量

%% 求拐点的算法（差分）
for i=4:length(y_)-4
    n1(i)=3*(loc(i,2))-(loc(i-1,2))-(loc(i-2,2))-(loc(i-3,2));
    n2(i)=loc(i+3,2)+loc(i+2,2)+loc(i+1,2)-3*(loc(i,2));
    d(i)=n2(i)-n1(i);
end
figure(4)
p=1:length(d);
plot(p,d);

%小波变换去噪
load noisbump
x=noisbump;
% x=d;
wname='sym6';
lev=5;
[c,l]=wavedec(x,lev,wname);
sigma=wnoisest(c,l,1);
alpha=2;
thr=wbmpen(c,l,sigma,alpha);
keepapp=1;
xd=wdencmp('gbl',c,l,wname,lev,thr,'s',keepapp);
figure(5)
subplot(211),plot(x)
subplot(212),plot(xd)
% dd=mean5_3(d,3);
% figure(5)
% pp=1:length(dd);
% plot(pp,dd);

%找峰值――方法一
x = linspace(0,1,length(d));
% NPeaks表示返回值中最大的峰值个数
[pks,locs,w,p] = findpeaks(d,x,...pks为输出峰值，locs为峰值位置
    'NPeaks',1,...输出的最大峰值个数
    'SortStr','ascend',...输出峰值是否进行排序，升序或降序或不进行排序
    'Threshold',0,...峰值之间的最小差值阈值，如数组为[]则未找到符合参数的阈值
    'WidthReference','halfprom',...
    ...'MiniPeakWidth','0',...最小的峰值宽度
    'MaxPeakWidth',Inf,...最大的峰值宽度
    'Annotate','peaks',...在定义输出后无效
    'MinPeakHeight',50);%最小峰值高度
figure(6)
findpeaks(d,x,'Annotate','extents')
text(locs+.02,pks,num2str((1:numel(pks))'))



M=max(dd);
[pks1,locs1] = findpeaks(abs(d),'minpeakdistance',floor(length(d)/10))%----mpd 设定两峰值间的最小间隔数
[pks2,locs2] = findpeaks(abs(d),'minpeakheight',M/5)%设定峰值的最小高度


LocN=x_1;
LocNy=y_1;
%同上，以末尾为原点
m2=0;
for i=1:length(y_)-1
    n2(i)=(max(loc(:,2))-loc(i,2))/(1-loc(i,1));
    if n2(i)>m2
        m2=n2(i);
        x_2=loc(i,1);
        y_2=loc(i,2);
    end
end
Loc1=x_2;
Loc1y=y_2;
m3=(Loc1y-min(loc(:,2)))/Loc1;
x_3=x_2;
y_3=y_2;
if Loc1<LocN
    %     for i=2:length(Loc1y)
    %         n3(i)=(loc(i,2)-min(loc(:,2)))/loc(i,1);%y/x
    %         if n3(i)<m3
    %             m3=n3(i);
    %             x_3=loc(i,1);
    %             y_3=loc(i,2);
    %         end
    %     end
    for i=1:length(Loc1y)-1
        n3(i)=(Loc1y-loc(i,2))/(Loc1-loc(i,1));
        if n3(i)>m3
            m3=n3(i);
            x_3=loc(i,1);
            y_3=loc(i,2);
        end
    end
end
Loc1=x_3;
Loc1y=y_3;
Loc1
LocN


%% 已经求出loc1和loc4，计算两个门限，并用两个门限计算出四个边界
figure(4)
f_=f(floor(length(Pxx_int)*Loc1):floor(length(Pxx_int)*LocN),1);
Pxx_=Pxx(floor(length(Pxx_int)*Loc1):floor(length(Pxx_int)*LocN),1);
plot(f_,Pxx_)
title("Welch法功率谱估计")
% MAX=max(Pxx_);
% min_P=min(Pxx_);
% H=MAX-min_P;
% Threshold2=MAX-(H+3);%门限
% bandldx3=find(Pxx>Threshold2);%寻找门限以上的点
% bandldx4=bandldx3./length(Pxx);%归一化坐标
% bw2=max(bandldx4)-min(bandldx4);%计算带宽
% Loc4=max(bandldx4)
% Loc1=min(bandldx4)

% Threshold1=MAX-H*1/3;%门限
% bandldx3=find(Pxx>Threshold1);%寻找门限以上的点
% bandldx4=bandldx3./length(Pxx);%归一化坐标
% bw2=max(bandldx4)-min(bandldx4);%计算带宽
% Loc3=max(bandldx4)
% Loc2=min(bandldx4)

