%一维峰值检测Matlab实现
%两个信号重叠，计算带宽的相关程序
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
x2=[x2(1:SubCarryNN2/2,1:20);zeros(SubCarryN2-SubCarryNN2,20);x2(SubCarryNN2/2+1:SubCarryNN2,1:20)];
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

x=d;

% row_acc=x;
% 
% figure(5)
% plot(x)
% L=length(x);
% 
% Y = fft(x);
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);% 实信号的功率谱是对称的,只用取一半即可
% P1(2:end-1) = 2*P1(2:end-1);
% f = Fs*(0:(L/2))/L;
% 
% plot(f,P1)
% title('原始信号直接进行傅里叶分析')
% grid on
% xlabel('Frequency (Hz)')
% ylabel('Amplitude')
% 
% [c,l] = wavedec(x,2,'db2');
% approx = appcoef(c,l,'db2');
% [cd1,cd2] = detcoef(c,l,[1 2]);
% figure(6)
% subplot(2,1,1)
% plot(x)
% title('原始信号')
% subplot(2,1,2)
% plot(approx)
% title('应变信号')
% % subplot(3,1,3)
% % plot(cd3)
% % title('冲击信号')
% figure(7)
% subplot(4,1,3)
% plot(cd2)
% title('Level 2 Detail Coefficients')
% subplot(4,1,4)
% plot(cd1)
% title('Level 1 Detail Coefficients')
% 
% % Hd=LowPass_Filter;
% % output=filter(Hd,x);
% % figure
% % plot(output)
% % 
% % blo = fir1(34,0.48,chebwin(35,30));
% % outlo = filter(blo,1,x);
% 
% % 使用峰值检测
% [pks,locs]=findpeaks(x,'Npeaks',1,'MinPeakHeight',250);
% jieyue=x(:,[locs-200:L]);
% figure(8)
% plot(jieyue)
% yingbian=x(:,[1:locs-200]);
% figure
% plot(yingbian)
% ylim([-50 300])
% figure
% plot(locs,pks)
% 
%% VPD方法

%%前三点均值滤波
% row_acc = abs(x);
row_acc = x;
m = length(row_acc);
row_acc1 = linspace(0,0,m);
row_acc1(1) = row_acc(1);
row_acc1(m) = row_acc(m);
for i=2:m-1
   row_acc1(i)=(row_acc(i-1) + row_acc(i)+row_acc(i+1))/3;%三点均值滤波
end
% figure;
% plot(row_acc1);

for i=m-1:-1:2
    row_acc(i) = (row_acc1(i-1) + row_acc1(i)+row_acc1(i+1))/3;%三点均值滤波
end

%%找到局部最小值和局部最大值及其对应的位置，波峰点、波谷点满足：
peaks = linspace(0,0,m);
valleys = linspace(0,0,m);
peakindexs = linspace(0,0,m);
valleyindexs = linspace(0,0,m);
peakindex = 1;
valleyindex = 1;

for i = 2:m-1
    if row_acc(i) >row_acc(i-1) && row_acc(i)>=row_acc(i+1)%找到局部最大
        peaks(peakindex)=row_acc(i);%峰值
        peakindexs(peakindex)=i;%峰值点的位置
        peakindex = peakindex+1;%峰值点的数目
    end
    if row_acc(i) < row_acc(i-1) && row_acc(i)<row_acc(i+1)%找到局部最小
        valleys(valleyindex)=row_acc(i);%谷值
        valleyindexs(valleyindex)=i;%谷值点的位置
        valleyindex=valleyindex+1;%谷值点的数目
    end
end
% 
%%计算VPD，VPD(n)表示第n个波谷点的值与第n个波峰点的值的差，VPD用来去掉那些假的波峰点
pcount = peakindex-1;
vcount = valleyindex-1;
peakindices = linspace(0,0,pcount);
for i = 1:pcount
    peakindices(i) = peakindexs(i);%峰值点的位置
end
valleyindices = linspace(0,0,vcount);
for i = 1:vcount
    valleyindices(i) = valleyindexs(i);%谷值点的位置
end
% figure(4);
% plot(row_acc,'-o', 'MarkerIndices',peakindices,'MarkerFaceColor','red','MarkerSize',5);
% hold on
% plot(row_acc,'-s', 'MarkerIndices',valleyindices,'MarkerFaceColor','green','MarkerSize',5);
if pcount>2 && vcount>2
    if peakindexs(1) < valleyindexs(1)
        peakindex=2;
    else
        peakindex=1;
    end
    vindex=1;
end

if peakindex == 2
    for i = 1:m-1
        peaks(i)=peaks(i+1);
    end
    pcount = pcount-1;
    pindex=1;
end

vpd = linspace(0,0,m);
vpd1 = linspace(0,0,m);
for i=1:pcount
    vpd(i) = peaks(i) - valleys(i);
end

dels = linspace(0,0, pcount);
peakindexs1 = linspace(0,0,pcount);
if pcount > 2
    lastcount=pcount;
    curcount = 1;
    while lastcount ~= curcount
        lastcount = curcount;
        del_count = 0;
        for i = 2:pcount-1
            if vpd(i) <= 0.0001 * (vpd(i-1) + vpd(i)+vpd(i+1)) / 3
                dels(i)=1;
            end
        end

        count = 1;
        for i = 1:pcount
          if dels(i) ~= 1
              vpd1(count) = vpd(i); 
              peakindexs1(count) = peakindexs(i);
              count = count+1;
          else
               del_count = del_count + 1;
               dels(i) = 0;
           end
        end
        pcount = pcount - del_count;
        for i = 1:pcount
            vpd(i) = vpd1(i); 
            peakindexs(i) = peakindexs1(i);
        end
        peakindexs(pcount+1) = 0;
        vpd(pcount+1) = 0;

        indices = linspace(0,0,pcount);
        for i = 1:pcount
            indices(i) = peakindexs1(i);
        end
        figure(5)
        plot(row_acc,'-o', 'MarkerIndices',peakindices,'MarkerFaceColor','red','MarkerSize',10);
        hold on
        plot(row_acc,'-o', 'MarkerIndices',valleyindices,'MarkerFaceColor','red','MarkerSize',10);
%         figure(6)
%         plot(abs(row_acc),'-o', 'MarkerIndices',peakindices,'MarkerFaceColor','red','MarkerSize',10);
%         hold on
%         plot(abs(row_acc),'-o', 'MarkerIndices',abs(valleyindices),'MarkerFaceColor','red','MarkerSize',10);
        curcount = pcount;
    end  
end
%输出的峰值点
peaks=row_acc(indices);
% [pks,locs]=findpeaks(row_acc,'Npeaks',1,'MinPeakHeight',20);
% jieyue=row_acc(:,[locs-200:L]);%从触发峰值点左边200处分离
% figure
% plot(jieyue)
% yingbian=x(:,[1:locs-200]);
% figure
% plot(yingbian)
% ylim([-50 300])
% 
% % figure
% % plot(locs,pks)
% 
figure(6)
% plot(row_acc,'-o', 'MarkerIndices',peakindices,'MarkerFaceColor','red','MarkerSize',10);
% hold on
% plot(row_acc,'-o', 'MarkerIndices',abs(valleyindices),'MarkerFaceColor','red','MarkerSize',10);
h=histogram(peaks,100);
M=mean(x);
s=std(x);%标准差
H=M+3*s;
peaks_dou=abs(peaks);
b=[];
[D, ia] = setdiff(peaks_dou, b);
peaks_2=peaks(ia);
% a1=find(peaks_dou>H)
% b1=peaks_dou(a1);
% [D, ia] = setdiff(peaks_dou, b1);
% peaks_2=peaks(ia);



figure(7)
scatter(1:length(peaks_2),peaks_2,'red','filled')