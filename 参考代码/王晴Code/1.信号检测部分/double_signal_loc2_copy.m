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

%% 求拐点的算法
for i=2:length(y_)
    n1(i)=(loc(i,2)-min(loc(:,2)))/loc(i,1);%y/x
    if n1(i)>m1
        m1=n1(i);
        x_1=loc(i,1);
        y_1=loc(i,2);
    end
end
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
% m3=(Loc1y-min(loc(:,2)))/Loc1;
% x_3=x_2;
% y_3=y_2;
% if Loc1<LocN
%     %     for i=2:length(Loc1y)
%     %         n3(i)=(loc(i,2)-min(loc(:,2)))/loc(i,1);%y/x
%     %         if n3(i)<m3
%     %             m3=n3(i);
%     %             x_3=loc(i,1);
%     %             y_3=loc(i,2);
%     %         end
%     %     end
%     for i=1:length(Loc1y)-1
%         n3(i)=(Loc1y-loc(i,2))/(Loc1-loc(i,1));
%         if n3(i)>m3
%             m3=n3(i);
%             x_3=loc(i,1);
%             y_3=loc(i,2);
%         end
%     end
% end
% Loc1=x_3;
% Loc1y=y_3;
Loc1
LocN


%% 已经求出loc1和loc4，计算两个门限，并用两个门限计算出四个边界
fprintf('\n========================================\n');
fprintf('边界检测结果\n');
fprintf('1号边界 (Loc1): %.4f\n', Loc1);
fprintf('4号边界 (LocN): %.4f\n', LocN);
fprintf('========================================\n\n');

% 使用Loc1和LocN截取数据（确保Loc1 < LocN）
if Loc1 >= LocN
    fprintf('警告：Loc1 >= LocN，交换边界位置\n');
    temp = Loc1;
    Loc1 = LocN;
    LocN = temp;
end

% 截取1号和4号边界之间的功率谱数据
idx_start = max(1, floor(length(Pxx)*Loc1));
idx_end = min(length(Pxx), floor(length(Pxx)*LocN));
f_=f(idx_start:idx_end,1);
Pxx_=Pxx(idx_start:idx_end,1);

figure(4)
plot(f_,Pxx_,'b-','LineWidth',2);
xlabel('归一化频率','FontSize',11,'FontWeight','bold');
ylabel('功率谱密度 (dB)','FontSize',11,'FontWeight','bold');
title('截取的功率谱（1号和4号边界之间）','FontSize',12,'FontWeight','bold');
grid on;
hold on;

% 计算上凸起部分的高度H
MAX=max(Pxx_);
min_P=min(Pxx_);
H=MAX-min_P;
fprintf('截取区域统计信息:\n');
fprintf('  最大值 (MAX): %.4f dB\n', MAX);
fprintf('  最小值 (min_P): %.4f dB\n', min_P);
fprintf('  高度 (H): %.4f dB\n', H);
fprintf('\n');

% 设置两个门限
Threshold1=MAX-H*1/3;  % 门限1：用于计算2号和3号边界
Threshold2=MAX-(H+3); % 门限2：用于验证1号和4号边界

fprintf('门限设置:\n');
fprintf('  门限1 (Threshold1 = MAX - H*1/3): %.4f dB\n', Threshold1);
fprintf('  门限2 (Threshold2 = MAX - H - 3): %.4f dB\n', Threshold2);
fprintf('\n');

% 在图上标注门限
yline(Threshold1,'--r','LineWidth',2,'Label','门限1','FontSize',10);
yline(Threshold2,'--g','LineWidth',2,'Label','门限2','FontSize',10);
yline(MAX,'--m','LineWidth',1.5,'Label','最大值','FontSize',9);
hold off;

%% 使用门限1计算2号和3号边界
bandldx1=find(Pxx>Threshold1);  % 寻找门限1以上的点
if ~isempty(bandldx1)
    bandldx1_norm=bandldx1./length(Pxx);  % 归一化坐标
    Loc2=min(bandldx1_norm);  % 2号边界：最小位置
    Loc3=max(bandldx1_norm);  % 3号边界：最大位置
    
    fprintf('使用门限1计算的边界:\n');
    fprintf('  2号边界 (Loc2): %.4f\n', Loc2);
    fprintf('  3号边界 (Loc3): %.4f\n', Loc3);
    fprintf('\n');
    
    % 在图上标注2号和3号边界
    figure(4);
    hold on;
    idx2 = floor(length(Pxx)*Loc2);
    idx3 = floor(length(Pxx)*Loc3);
    if idx2 >= 1 && idx2 <= length(Pxx)
        plot(f(idx2), Pxx(idx2), 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'red', ...
            'LineWidth', 2, 'DisplayName', '2号边界');
    end
    if idx3 >= 1 && idx3 <= length(Pxx)
        plot(f(idx3), Pxx(idx3), 'go', 'MarkerSize', 12, 'MarkerFaceColor', 'green', ...
            'LineWidth', 2, 'DisplayName', '3号边界');
    end
    hold off;
    legend('功率谱','门限1','门限2','最大值','2号边界','3号边界','Location','best');
else
    fprintf('警告：未找到超过门限1的点，无法计算2号和3号边界\n');
    Loc2 = NaN;
    Loc3 = NaN;
end

%% 使用门限2验证或重新计算1号和4号边界
bandldx2=find(Pxx>Threshold2);  % 寻找门限2以上的点
if ~isempty(bandldx2)
    bandldx2_norm=bandldx2./length(Pxx);  % 归一化坐标
    Loc1_verify=min(bandldx2_norm);  % 验证的1号边界
    Loc4_verify=max(bandldx2_norm);   % 验证的4号边界
    
    fprintf('使用门限2验证的边界:\n');
    fprintf('  1号边界 (Loc1_verify): %.4f (原始值: %.4f)\n', Loc1_verify, Loc1);
    fprintf('  4号边界 (Loc4_verify): %.4f (原始值: %.4f)\n', Loc4_verify, LocN);
    fprintf('\n');
    
    % 在图上标注验证的边界
    figure(4);
    hold on;
    idx1_v = floor(length(Pxx)*Loc1_verify);
    idx4_v = floor(length(Pxx)*Loc4_verify);
    if idx1_v >= 1 && idx1_v <= length(Pxx)
        plot(f(idx1_v), Pxx(idx1_v), 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'none', ...
            'LineWidth', 2, 'DisplayName', '1号边界(验证)');
    end
    if idx4_v >= 1 && idx4_v <= length(Pxx)
        plot(f(idx4_v), Pxx(idx4_v), 'gs', 'MarkerSize', 10, 'MarkerFaceColor', 'none', ...
            'LineWidth', 2, 'DisplayName', '4号边界(验证)');
    end
    hold off;
end

%% 计算带宽
if ~isnan(Loc2) && ~isnan(Loc3)
    bw_overlap = Loc3 - Loc2;  % 重叠区域带宽
    fprintf('带宽计算结果:\n');
    fprintf('  重叠区域带宽 (Loc3 - Loc2): %.4f\n', bw_overlap);
    fprintf('\n');
end

%% 绘制完整的功率谱图，标注所有四个边界
figure(6)
plot(f,Pxx,'b-','LineWidth',2);
hold on;
yline(Threshold1,'--r','LineWidth',1.5,'Label','门限1','FontSize',9);
yline(Threshold2,'--g','LineWidth',1.5,'Label','门限2','FontSize',9);

% 标注四个边界
idx1_plot = floor(length(Pxx)*Loc1);
idx2_plot = floor(length(Pxx)*Loc2);
idx3_plot = floor(length(Pxx)*Loc3);
idx4_plot = floor(length(Pxx)*LocN);

if idx1_plot >= 1 && idx1_plot <= length(Pxx)
    plot(f(idx1_plot), Pxx(idx1_plot), 'mo', 'MarkerSize', 12, 'MarkerFaceColor', 'magenta', ...
        'LineWidth', 2, 'DisplayName', sprintf('1号边界 (%.4f)', Loc1));
end
if ~isnan(Loc2) && idx2_plot >= 1 && idx2_plot <= length(Pxx)
    plot(f(idx2_plot), Pxx(idx2_plot), 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'red', ...
        'LineWidth', 2, 'DisplayName', sprintf('2号边界 (%.4f)', Loc2));
end
if ~isnan(Loc3) && idx3_plot >= 1 && idx3_plot <= length(Pxx)
    plot(f(idx3_plot), Pxx(idx3_plot), 'go', 'MarkerSize', 12, 'MarkerFaceColor', 'green', ...
        'LineWidth', 2, 'DisplayName', sprintf('3号边界 (%.4f)', Loc3));
end
if idx4_plot >= 1 && idx4_plot <= length(Pxx)
    plot(f(idx4_plot), Pxx(idx4_plot), 'co', 'MarkerSize', 12, 'MarkerFaceColor', 'cyan', ...
        'LineWidth', 2, 'DisplayName', sprintf('4号边界 (%.4f)', LocN));
end

hold off;
xlabel('归一化频率','FontSize',12,'FontWeight','bold');
ylabel('功率谱密度 (dB)','FontSize',12,'FontWeight','bold');
title('完整功率谱 - 四个边界标注','FontSize',13,'FontWeight','bold');
legend('Location','best','FontSize',10);
grid on;
grid minor;

fprintf('========================================\n');
fprintf('所有边界检测完成！\n');
fprintf('========================================\n');

