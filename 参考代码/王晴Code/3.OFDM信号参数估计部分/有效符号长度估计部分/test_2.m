%第二个实验，循环前缀比例对OFDM符号长度估计方法的性能影响实验

%实现基本的OFDM通信系统仿真，包括QPSK调制、串并变换、IFFT、并串变换、
%插入保护间隔AWGN信道传输以及相应的逆过程。
%每个符号采用1024个子载波，每个载波传输6个符号，添加的保护间隔为256个。
%基带调制采用QPSK，接收端信噪比为10db
clc
clear all
close all
SubCarryNN=768;      %有效子载波数
SubCarryN=1024;      %子载波数
ratio=1/32;
fftLen=1024;         %FFT长度为128
SymbN=6;            %一帧中OFDM符号个数、每个符号传输的子载波数、6×128=768
GuardLen=SubCarryN*ratio;    %保护时隙的长度
SNR=-15:5:15;             %信噪比取值，以db为单位

%%----------------------------------------发射端-----------------------------------------------------%%
%输入比特序列长度=子载波数×每载波符号数×每符号比特数=768×6×2=9216
for index=1:length(SNR)
for p=1:2000
SignalLen=SubCarryNN*SymbN*2;           %比特序列长度，9216
Signal=round(rand(1,SignalLen));       %输出带调制的二进制比特流,0,1、round是为了产生0，1
ParaBitSig=[];

%length(ParaBitSig)竟然是128
%这是一个把一行数变成X行Y列的数组的方法
for i=1:SubCarryNN
    for j=1:SymbN*2
        ParaBitSig(i,j)=Signal(i*j);   %串/并转换为行数SubCarryNN（有效子载波），列数为2*SymbN（每载波的比特数）
    end
end

%进行QPSK数据调制，将数据分为两个通道
%这是一个QPSK的调制方法
for j=1:SymbN
    ich(:,j)=ParaBitSig(:,2*j-1); %同相分量，ParaBitSig奇数列
    qch(:,j)=ParaBitSig(:,2*j);   %正交分量，ParaBitSig偶数列
end
kmod=1./sqrt(2);   %做一个根号二
ich0=ich.*2-1;     %这俩变成正负一的形式，1→1，0→-1
qch0=qch.*2-1;
ich1=ich0.*kmod;   %根号二了
qch1=qch0.*kmod;
x=ich1+qch1.*sqrt(-1);      %产生复信号,768行，6列


% x=[zeros(32,6);x;zeros(32,6)];%在中间补零，补了256行
x=[x(1:384,1:6);zeros(256,6);x(385:768,1:6)];
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
ReData=awgn(TrData,SNR(index),'measured');%信道，加入高斯白噪声  %每一列是一个ofdm符号

%下面尝试估计
%test1,第一个实验，基于循环前缀的OFDM符号有效长度算法的估计原理实验
%估计出了有效符号长度
r=ReData;
L=1500;
M=256;
R=zeros(1,L);
E=zeros(1,L);
MM=zeros(1,L);
for k=1:L
    for l=1:M
        R(k)=R(k)+r(l)*conj(r(l+k));
        E(k)=E(k)+abs(r(l))^2+abs(r(l+k))^2;
    end
    E(k)=E(k)/2;
    MM(k)=abs(R(k))/E(k);
end
Max=max(MM);
for k=1:L
    if MM(k)==Max
        K(p)=k;
    end
end
% k=1:L;
% plot(k,MM(k))
% xlabel("数据位移量")
% ylabel("归一化的相关值")
err_32(p)=SubCarryN-K(p);
end
mse_32(index)=mean(abs(err_32).^2)/1024^2;
end


%% radio=1/16
ratio=1/16;
GuardLen=SubCarryN*ratio;    %保护时隙的长度
for index=1:length(SNR)
for p=1:20
SignalLen=SubCarryNN*SymbN*2;           %比特序列长度，9216
Signal=round(rand(1,SignalLen));       %输出带调制的二进制比特流,0,1、round是为了产生0，1
ParaBitSig=[];

%length(ParaBitSig)竟然是128
%这是一个把一行数变成X行Y列的数组的方法
for i=1:SubCarryNN
    for j=1:SymbN*2
        ParaBitSig(i,j)=Signal(i*j);   %串/并转换为行数SubCarryNN（有效子载波），列数为2*SymbN（每载波的比特数）
    end
end

%进行QPSK数据调制，将数据分为两个通道
%这是一个QPSK的调制方法
for j=1:SymbN
    ich(:,j)=ParaBitSig(:,2*j-1); %同相分量，ParaBitSig奇数列
    qch(:,j)=ParaBitSig(:,2*j);   %正交分量，ParaBitSig偶数列
end
kmod=1./sqrt(2);   %做一个根号二
ich0=ich.*2-1;     %这俩变成正负一的形式，1→1，0→-1
qch0=qch.*2-1;
ich1=ich0.*kmod;   %根号二了
qch1=qch0.*kmod;
x=ich1+qch1.*sqrt(-1);      %产生复信号,768行，6列


% x=[zeros(32,6);x;zeros(32,6)];%在中间补零，补了256行
x=[x(1:384,1:6);zeros(256,6);x(385:768,1:6)];
% x=[x(1:128,1:6);zeros(768,6);x(129:256,1:6)];
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
ReData=awgn(TrData,SNR(index),'measured');%信道，加入高斯白噪声  %每一列是一个ofdm符号

%下面尝试估计
%test1,第一个实验，基于循环前缀的OFDM符号有效长度算法的估计原理实验
%估计出了有效符号长度
r=ReData;
L=1500;
M=256;
R=zeros(1,L);
E=zeros(1,L);
MM=zeros(1,L);
for k=1:L
    for l=1:M
        R(k)=R(k)+r(l)*conj(r(l+k));
        E(k)=E(k)+abs(r(l))^2+abs(r(l+k))^2;
    end
    E(k)=E(k)/2;
    MM(k)=abs(R(k))/E(k);
end
Max=max(MM);
for k=1:L
    if MM(k)==Max
        K(p)=k;
    end
end
% k=1:L;
% plot(k,MM(k))
% xlabel("数据位移量")
% ylabel("归一化的相关值")
err_16(p)=SubCarryN-K(p);
end
mse_16(index)=mean(abs(err_16).^2)/1024^2;
end


%% radio=1/8
ratio=1/8;
GuardLen=SubCarryN*ratio;    %保护时隙的长度
for index=1:length(SNR)
for p=1:20
SignalLen=SubCarryNN*SymbN*2;           %比特序列长度，9216
Signal=round(rand(1,SignalLen));       %输出带调制的二进制比特流,0,1、round是为了产生0，1
ParaBitSig=[];

%length(ParaBitSig)竟然是128
%这是一个把一行数变成X行Y列的数组的方法
for i=1:SubCarryNN
    for j=1:SymbN*2
        ParaBitSig(i,j)=Signal(i*j);   %串/并转换为行数SubCarryNN（有效子载波），列数为2*SymbN（每载波的比特数）
    end
end

%进行QPSK数据调制，将数据分为两个通道
%这是一个QPSK的调制方法
for j=1:SymbN
    ich(:,j)=ParaBitSig(:,2*j-1); %同相分量，ParaBitSig奇数列
    qch(:,j)=ParaBitSig(:,2*j);   %正交分量，ParaBitSig偶数列
end
kmod=1./sqrt(2);   %做一个根号二
ich0=ich.*2-1;     %这俩变成正负一的形式，1→1，0→-1
qch0=qch.*2-1;
ich1=ich0.*kmod;   %根号二了
qch1=qch0.*kmod;
x=ich1+qch1.*sqrt(-1);      %产生复信号,768行，6列


% x=[zeros(32,6);x;zeros(32,6)];%在中间补零，补了256行
x=[x(1:384,1:6);zeros(256,6);x(385:768,1:6)];
% x=[x(1:128,1:6);zeros(768,6);x(129:256,1:6)];
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
ReData=awgn(TrData,SNR(index),'measured');%信道，加入高斯白噪声  %每一列是一个ofdm符号

%下面尝试估计
%test1,第一个实验，基于循环前缀的OFDM符号有效长度算法的估计原理实验
%估计出了有效符号长度
r=ReData;
L=1500;
M=256;
R=zeros(1,L);
E=zeros(1,L);
MM=zeros(1,L);
for k=1:L
    for l=1:M
        R(k)=R(k)+r(l)*conj(r(l+k));
        E(k)=E(k)+abs(r(l))^2+abs(r(l+k))^2;
    end
    E(k)=E(k)/2;
    MM(k)=abs(R(k))/E(k);
end
Max=max(MM);
for k=1:L
    if MM(k)==Max
        K(p)=k;
    end
end
% k=1:L;
% plot(k,MM(k))
% xlabel("数据位移量")
% ylabel("归一化的相关值")
err_8(p)=SubCarryN-K(p);
end
mse_8(index)=mean(abs(err_8).^2)/1024^2;
end


%% radio=1/4
ratio=1/4;
GuardLen=SubCarryN*ratio;    %保护时隙的长度
for index=1:length(SNR)
for p=1:20
SignalLen=SubCarryNN*SymbN*2;           %比特序列长度，9216
Signal=round(rand(1,SignalLen));       %输出带调制的二进制比特流,0,1、round是为了产生0，1
ParaBitSig=[];

%length(ParaBitSig)竟然是128
%这是一个把一行数变成X行Y列的数组的方法
for i=1:SubCarryNN
    for j=1:SymbN*2
        ParaBitSig(i,j)=Signal(i*j);   %串/并转换为行数SubCarryNN（有效子载波），列数为2*SymbN（每载波的比特数）
    end
end

%进行QPSK数据调制，将数据分为两个通道
%这是一个QPSK的调制方法
for j=1:SymbN
    ich(:,j)=ParaBitSig(:,2*j-1); %同相分量，ParaBitSig奇数列
    qch(:,j)=ParaBitSig(:,2*j);   %正交分量，ParaBitSig偶数列
end
kmod=1./sqrt(2);   %做一个根号二
ich0=ich.*2-1;     %这俩变成正负一的形式，1→1，0→-1
qch0=qch.*2-1;
ich1=ich0.*kmod;   %根号二了
qch1=qch0.*kmod;
x=ich1+qch1.*sqrt(-1);      %产生复信号,768行，6列


% x=[zeros(32,6);x;zeros(32,6)];%在中间补零，补了256行
x=[x(1:384,1:6);zeros(256,6);x(385:768,1:6)];
% x=[x(1:128,1:6);zeros(768,6);x(129:256,1:6)];
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
ReData=awgn(TrData,SNR(index),'measured');%信道，加入高斯白噪声  %每一列是一个ofdm符号

%下面尝试估计
%test1,第一个实验，基于循环前缀的OFDM符号有效长度算法的估计原理实验
%估计出了有效符号长度
r=ReData;
L=1500;
M=256;
R=zeros(1,L);
E=zeros(1,L);
MM=zeros(1,L);
for k=1:L
    for l=1:M
        R(k)=R(k)+r(l)*conj(r(l+k));
        E(k)=E(k)+abs(r(l))^2+abs(r(l+k))^2;
    end
    E(k)=E(k)/2;
    MM(k)=abs(R(k))/E(k);
end
Max=max(MM);
for k=1:L
    if MM(k)==Max
        K(p)=k;
    end
end
% k=1:L;
% plot(k,MM(k))
% xlabel("数据位移量")
% ylabel("归一化的相关值")
err_4(p)=SubCarryN-K(p);
end
mse_4(index)=mean(abs(err_4).^2)/(1024^2);
end

%% 画图
figure
semilogy(SNR,mse_32,"b-*",SNR,mse_16,"r-o",SNR,mse_8,"k-s",SNR,mse_4,"c-x" );
xlabel("SNR");
ylabel("NMSE");
axis([-15,15,0.0001,1e2])
% legend("16APSK")
legend('radio=1/32','radio=1/16','radio=1/8','radio=1/4');  