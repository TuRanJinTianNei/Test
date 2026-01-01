%两个信号重叠，计算带宽的相关程序
%-----------------------------------修改说明-------------------------------------------------------
% 1. 核心检测算法：使用差分法 + VPD (Valley-Peak Difference) 提取四个边界。
% 2. 新增功能：根据信号生成参数计算“理论边界值”。
% 3. 输出对比：输出检测值、理论值及绝对误差。
%--------------------------------------------------------------------------------------------------
clc
clear all
close all
%% 信号1参数
SubCarryNN=256;      %有效子载波数
SubCarryN=1024;      %FFT长度
ratio=1/4;           %循环前缀比例
fftLen=1024;         
SymbN=20;            
GuardLen=SubCarryN*ratio;   
SNR=10;             
% 生成信号1 (中心频率为0，基带)
SignalLen=SubCarryNN*SymbN*2;           
Signal=round(rand(1,SignalLen));       
ParaBitSig=reshape(Signal,SubCarryNN,SymbN*2);
for j=1:SymbN
    ich(:,j)=ParaBitSig(:,2*j-1); 
    qch(:,j)=ParaBitSig(:,2*j);   
end
kmod=1./sqrt(2);   
ich0=ich.*2-1;    
qch0=qch.*2-1;
ich1=ich0.*kmod;   
qch1=qch0.*kmod;
x=ich1+qch1.*sqrt(-1);      
x=[x(1:128,1:20);zeros(768,20);x(129:256,1:20)]; % 标准OFDM映射
y=ifft(x);   
ich2=real(y);   
qch2=imag(y);   
ich3=[ich2(fftLen-GuardLen+1:fftLen,:);ich2];
qch3=[qch2(fftLen-GuardLen+1:fftLen,:);qch2];
ich4=reshape(ich3,1,(fftLen+GuardLen)*SymbN);
qch4=reshape(qch3,1,(fftLen+GuardLen)*SymbN);
TrData=ich4+qch4.*sqrt(-1);  
ReData=awgn(TrData,SNR,'measured');

%% 信号2参数
SubCarryNN2=768;      %有效子载波数
SubCarryN2=1024;      
SymbN2=20;            
% 生成信号2
SignalLen2=SubCarryNN2*SymbN2*2;           
Signal2=round(rand(1,SignalLen2));       
ParaBitSig2=reshape(Signal2,SubCarryNN2,SymbN2*2);
for j=1:SymbN2
    ich22(:,j)=ParaBitSig2(:,2*j-1); 
    qch22(:,j)=ParaBitSig2(:,2*j);   
end
kmod=1./sqrt(2);   
ich02=ich22.*2-1;     
qch02=qch22.*2-1;
ich12=ich02.*kmod;   
qch12=qch02.*kmod;
x2=ich12+qch12.*sqrt(-1);      
x2=[x2(1:384,1:20);zeros(256,20);x2(385:768,1:20)];
y2=ifft(x2);   
ich22=real(y2);   
qch22=imag(y2);   
ich32=[ich22(fftLen-GuardLen+1:fftLen,:);ich22];
qch32=[qch22(fftLen-GuardLen+1:fftLen,:);qch22];
ich42=reshape(ich32,1,(fftLen+GuardLen)*SymbN2);
qch42=reshape(qch32,1,(fftLen+GuardLen)*SymbN2);
TrData42=ich42+qch42.*sqrt(-1);  
ReData42=awgn(TrData42,SNR,'measured');

% 关键步骤：信号2频移
% exp(1j * (0.1) * pi * 2 * n) 意味着归一化频移为 0.1 (相对于采样率fs)
freq_shift_ratio = 0.1; 
ReData_same = ReData42.*exp(1j * freq_shift_ratio * pi * 2 .* (1:length(ReData42)));

%% =========================================================================
%% 步骤0: 绘制重叠前的两个信号频谱
%% =========================================================================
fs=1024;
Fs=1;
nfft=1024;
window=hanning(150);
noverlap=50;

% 计算信号1的功率谱（基带，中心频率为0）
[Pxx1_sig1, ~] = pwelch(ReData, window, noverlap, length(ReData), Fs);
Pxx_sig1 = 10*log10(fftshift(Pxx1_sig1));
% 生成归一化频率轴（与原始代码保持一致）
f_sig1 = linspace(-fs/2, fs/2, length(Pxx_sig1));
f_sig1 = f_sig1 / (2*max(f_sig1));  % 归一化频率，范围[-0.5, 0.5]

% 计算信号2的功率谱（频移前，基带）
[Pxx1_sig2, ~] = pwelch(ReData42, window, noverlap, length(ReData42), Fs);
Pxx_sig2 = 10*log10(fftshift(Pxx1_sig2));
% 生成归一化频率轴
f_sig2 = linspace(-fs/2, fs/2, length(Pxx_sig2));
f_sig2 = f_sig2 / (2*max(f_sig2));  % 归一化频率，范围[-0.5, 0.5]

% 计算信号2频移后的功率谱
[Pxx1_sig2_shifted, ~] = pwelch(ReData_same, window, noverlap, length(ReData_same), Fs);
Pxx_sig2_shifted = 10*log10(fftshift(Pxx1_sig2_shifted));
% 生成归一化频率轴
f_sig2_shifted = linspace(-fs/2, fs/2, length(Pxx_sig2_shifted));
f_sig2_shifted = f_sig2_shifted / (2*max(f_sig2_shifted));  % 归一化频率，范围[-0.5, 0.5]

% 绘制重叠前的两个信号频谱对比图
figure('Name', '重叠前的两个信号频谱', 'Position', [100, 100, 1400, 600]);

% 子图1：信号1和信号2（频移前）的对比
subplot(1, 2, 1);
plot(f_sig1, Pxx_sig1, 'b-', 'LineWidth', 1.5, 'DisplayName', sprintf('信号1（基带，%d子载波）', SubCarryNN));
hold on;
plot(f_sig2, Pxx_sig2, 'r-', 'LineWidth', 1.5, 'DisplayName', sprintf('信号2（基带，%d子载波）', SubCarryNN2));
% 标记理论边界
Theory_Sig1_BW_temp = SubCarryNN / SubCarryN;
Theory_Sig1_Left_temp = 0.5 - Theory_Sig1_BW_temp/2;
Theory_Sig1_Right_temp = 0.5 + Theory_Sig1_BW_temp/2;
Theory_Sig2_BW_temp = SubCarryNN2 / SubCarryN2;
Theory_Sig2_Left_temp = 0.5 - Theory_Sig2_BW_temp/2;
Theory_Sig2_Right_temp = 0.5 + Theory_Sig2_BW_temp/2;
xline(Theory_Sig1_Left_temp, '--b', 'LineWidth', 1, 'DisplayName', '信号1左边界');
xline(Theory_Sig1_Right_temp, '--b', 'LineWidth', 1, 'DisplayName', '信号1右边界');
xline(Theory_Sig2_Left_temp, '--r', 'LineWidth', 1, 'DisplayName', '信号2左边界');
xline(Theory_Sig2_Right_temp, '--r', 'LineWidth', 1, 'DisplayName', '信号2右边界');
xlabel('归一化频率');
ylabel('功率谱密度 (dB)');
title('重叠前的两个信号频谱（频移前）');
legend('Location', 'best');
grid on;
hold off;

% 子图2：信号1和信号2（频移后）的对比
subplot(1, 2, 2);
plot(f_sig1, Pxx_sig1, 'b-', 'LineWidth', 1.5, 'DisplayName', sprintf('信号1（基带，%d子载波）', SubCarryNN));
hold on;
plot(f_sig2_shifted, Pxx_sig2_shifted, 'r-', 'LineWidth', 1.5, 'DisplayName', sprintf('信号2（频移后，%d子载波）', SubCarryNN2));
% 标记理论边界
Theory_Sig1_BW_temp = SubCarryNN / SubCarryN;
Theory_Sig1_Left_temp = 0.5 - Theory_Sig1_BW_temp/2;
Theory_Sig1_Right_temp = 0.5 + Theory_Sig1_BW_temp/2;
Theory_Sig2_Center_temp = 0.5 + freq_shift_ratio;
Theory_Sig2_BW_temp = SubCarryNN2 / SubCarryN2;
Theory_Sig2_Left_shifted = Theory_Sig2_Center_temp - Theory_Sig2_BW_temp/2;
Theory_Sig2_Right_shifted = Theory_Sig2_Center_temp + Theory_Sig2_BW_temp/2;
xline(Theory_Sig1_Left_temp, '--b', 'LineWidth', 1, 'DisplayName', '信号1左边界');
xline(Theory_Sig1_Right_temp, '--b', 'LineWidth', 1, 'DisplayName', '信号1右边界');
xline(Theory_Sig2_Left_shifted, '--r', 'LineWidth', 1, 'DisplayName', '信号2左边界');
xline(Theory_Sig2_Right_shifted, '--r', 'LineWidth', 1, 'DisplayName', '信号2右边界');
xlabel('归一化频率');
ylabel('功率谱密度 (dB)');
title(sprintf('重叠前的两个信号频谱（信号2频移后，频移比=%.2f）', freq_shift_ratio));
legend('Location', 'best');
grid on;
hold off;

fprintf('\n================== 重叠前的信号参数 ==================\n');
fprintf('信号1：%d个子载波，基带信号（中心频率=0）\n', SubCarryNN);
fprintf('信号2：%d个子载波，频移比=%.2f\n', SubCarryNN2, freq_shift_ratio);
fprintf('======================================================\n\n');

ReData=ReData+ReData_same; % 信号叠加

%% =========================================================================
%% 步骤1: Welch功率谱估计
%% =========================================================================
fs=1024;
f=linspace(-fs/2,fs/2,25600);
f=f/(2*max(f));
Fs=1;
nfft=1024;
window=hanning(150);
noverlap=50;
[Pxx1,f]=pwelch(ReData,window,noverlap,length(ReData),Fs);
Pxx=10*log10(fftshift(Pxx1));

figure(1)
plot(f,Pxx, 'b-', 'LineWidth', 1.5)
xlabel('归一化频率')
ylabel('功率谱密度 (dB)')
title('步骤1: Welch法功率谱估计')
grid on

%% =========================================================================
%% 步骤2: 计算功率分布函数（累积分布函数）
%% =========================================================================
Pxx=Pxx+abs(min(Pxx));
Pxx_int=zeros(size(Pxx,1),1);
for x=1:1:size(Pxx,1)
    if x==1
        Pxx_int(x)=(1/(size(Pxx,1)-1))*Pxx(x+1);
    else
        Pxx_int(x)=Pxx_int(x-1)+(1/(size(Pxx,1)-1))*Pxx(x);
    end
end

y_=Pxx_int;
x_=(linspace(0,1,length(Pxx_int)))';
loc=[x_,y_];

% 绘制功率分布函数
figure(2)
plot(x_, y_, 'b-', 'LineWidth', 1.5)
xlabel('归一化频率')
ylabel('累积功率分布')
title('步骤2: 功率分布函数（累积分布函数）')
grid on

%% =========================================================================
%% 算法部分：差分法 + VPD 边界检测
%% =========================================================================

% 1. 全局差分
d_all = zeros(1, length(y_));
for i=4:length(y_)-4
    n1(i)=3*(loc(i,2))-(loc(i-1,2))-(loc(i-2,2))-(loc(i-3,2));
    n2(i)=loc(i+3,2)+loc(i+2,2)+loc(i+1,2)-3*(loc(i,2));
    d_all(i)=n2(i)-n1(i);
end

% 绘制全局差分结果
figure(3)
plot(x_, d_all, 'b-', 'LineWidth', 1.5)
xlabel('归一化频率')
ylabel('差分值')
title('步骤3: 全局差分结果')
grid on
hold on
plot(x_, zeros(size(x_)), 'k--', 'LineWidth', 0.5)
hold off

% 2. VPD 方法求峰值
x_vpd = abs(d_all); 
row_acc = x_vpd;
m = length(row_acc);

% 绘制VPD处理前的结果
figure(4)
subplot(2,1,1)
plot(x_, row_acc, 'r-', 'LineWidth', 1.5)
xlabel('归一化频率')
ylabel('|差分值|')
title('步骤4a: VPD处理前（取绝对值后）')
grid on

% 均值滤波
row_acc1 = row_acc;
for i=2:m-1
    row_acc1(i)=(row_acc(i-1) + row_acc(i)+row_acc(i+1))/3;
end
for i=m-1:-1:2
    row_acc(i) = (row_acc1(i-1) + row_acc1(i)+row_acc1(i+1))/3;
end

% 绘制均值滤波后的结果
subplot(2,1,2)
plot(x_, row_acc, 'g-', 'LineWidth', 1.5)
xlabel('归一化频率')
ylabel('|差分值| (滤波后)')
title('步骤4b: 均值滤波后的结果')
grid on

% 局部极值
peaks = zeros(1,m); valleys = zeros(1,m);
peakindexs = zeros(1,m); valleyindexs = zeros(1,m);
peakindex = 1; valleyindex = 1;

for i = 2:m-1
    if row_acc(i) >row_acc(i-1) && row_acc(i)>=row_acc(i+1)
        peaks(peakindex)=row_acc(i);
        peakindexs(peakindex)=i;
        peakindex = peakindex+1;
    end
    if row_acc(i) < row_acc(i-1) && row_acc(i)<row_acc(i+1)
        valleys(valleyindex)=row_acc(i);
        valleyindexs(valleyindex)=i;
        valleyindex=valleyindex+1;
    end
end

% 绘制峰值和谷值检测结果
figure(5)
plot(x_, row_acc, 'b-', 'LineWidth', 1.5)
hold on
valid_peak_idx = peakindexs(1:peakindex-1);
valid_peak_idx = valid_peak_idx(valid_peak_idx > 0);
valid_valley_idx = valleyindexs(1:valleyindex-1);
valid_valley_idx = valid_valley_idx(valid_valley_idx > 0);
if ~isempty(valid_peak_idx)
    plot(x_(valid_peak_idx), row_acc(valid_peak_idx), 'r^', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', '峰值')
end
if ~isempty(valid_valley_idx)
    plot(x_(valid_valley_idx), row_acc(valid_valley_idx), 'gv', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'DisplayName', '谷值')
end
xlabel('归一化频率')
ylabel('|差分值| (滤波后)')
title('步骤5: 峰值和谷值检测结果')
legend('Location', 'best')
grid on
hold off

% VPD迭代筛选
pcount = peakindex-1;
vpd = peaks(1:pcount) - valleys(1:pcount);
indices = peakindexs(1:pcount);
indices = indices(indices>0);

% 3. 选取最强4个峰值作为检测边界
if length(indices) >= 4
    peak_values = row_acc(indices);
    [~, sort_idx] = sort(peak_values, 'descend');
    top4_indices = indices(sort_idx(1:4));
    boundary_indices = sort(top4_indices); % 按频率位置排序
    
    Detected_Locs = [loc(boundary_indices(1), 1), ...
                     loc(boundary_indices(2), 1), ...
                     loc(boundary_indices(3), 1), ...
                     loc(boundary_indices(4), 1)];
else
    fprintf('警告：检测到的边界点不足4个，无法进行有效对比。\n');
    Detected_Locs = zeros(1,4); 
end

% 绘制边界检测结果（在VPD结果上标注）
figure(6)
plot(x_, row_acc, 'b-', 'LineWidth', 1.5)
hold on
if length(indices) >= 4
    plot(x_(boundary_indices), row_acc(boundary_indices), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'LineWidth', 2, 'DisplayName', '检测边界')
    for i = 1:4
        xline(x_(boundary_indices(i)), '--r', 'LineWidth', 1.5);
    end
end
xlabel('归一化频率')
ylabel('|差分值| (滤波后)')
title('步骤6: 边界检测结果（前4个最强峰值）')
legend('Location', 'best')
grid on
hold off

%% =========================================================================
%% 理论边界计算与对比
%% =========================================================================

% --- 理论计算说明 ---
% 归一化频率范围 [0, 1] 对应物理频率 [-fs/2, fs/2] (经fftshift后)
% 中心点 0.5 对应 DC (0 Hz)

% 信号1 (基带):
% 带宽占比 = 256 / 1024 = 0.25
% 范围: [0.5 - 0.25/2,  0.5 + 0.25/2] = [0.375, 0.625]
Theory_Sig1_BW = SubCarryNN / SubCarryN;
Theory_Sig1_Left  = 0.5 - Theory_Sig1_BW/2;
Theory_Sig1_Right = 0.5 + Theory_Sig1_BW/2;

% 信号2 (频移):
% 带宽占比 = 768 / 1024 = 0.75
% 频移占比 = 0.1 (向右)
% 中心位置 = 0.5 + 0.1 = 0.6
% 范围: [0.6 - 0.75/2, 0.6 + 0.75/2] = [0.225, 0.975]
Theory_Sig2_BW = SubCarryNN2 / SubCarryN2;
Theory_Sig2_Center = 0.5 + freq_shift_ratio;
Theory_Sig2_Left  = Theory_Sig2_Center - Theory_Sig2_BW/2;
Theory_Sig2_Right = Theory_Sig2_Center + Theory_Sig2_BW/2;

% 汇总所有理论边界并排序 (从小到大)
Theory_Locs_Raw = [Theory_Sig1_Left, Theory_Sig1_Right, Theory_Sig2_Left, Theory_Sig2_Right];
Theory_Locs = sort(Theory_Locs_Raw);

%% =========================================================================
%% 步骤7: 最终结果对比图
%% =========================================================================
figure(7)
plot(f, Pxx, 'b-', 'LineWidth', 1.5)
hold on
% 画出检测到的边界（红色虚线）
for i=1:4
    xline(Detected_Locs(i), '--r', 'LineWidth', 2, 'DisplayName', sprintf('检测边界%d', i));
end
% 画出理论边界（绿色实线）
for i=1:4
    xline(Theory_Locs(i), '-g', 'LineWidth', 2, 'DisplayName', sprintf('理论边界%d', i));
end
xlabel('归一化频率')
ylabel('功率谱密度 (dB)')
title('步骤7: 最终结果对比（功率谱 + 检测边界 vs 理论边界）')
legend('功率谱', '检测边界', '理论边界', 'Location', 'best')
grid on
hold off

%% 输出详细对比表格
fprintf('\n================== 边界检测精度分析 ==================\n');
fprintf('%-10s %-12s %-12s %-12s\n', '边界序号', '检测值(Loc)', '理论值(Theory)', '绝对误差');
fprintf('------------------------------------------------------\n');

Labels = {'Sig2 Left', 'Sig1 Left', 'Sig1 Right', 'Sig2 Right'}; % 基于理论排序的推断标签

for i = 1:4
    error_val = abs(Detected_Locs(i) - Theory_Locs(i));
    fprintf('%d (%s)   %-12.4f %-12.4f %-12.4f\n', ...
        i, Labels{i}, Detected_Locs(i), Theory_Locs(i), error_val);
end
fprintf('======================================================\n');

% 计算平均误差
mean_error = mean(abs(Detected_Locs - Theory_Locs));
fprintf('平均检测误差: %.4f (归一化频率)\n\n', mean_error);