% 两个信号部分重叠，高精度检测版本 (集成能量门控修正版)
%-----------------------------------修改报告--------------------------------------------------
% 1. 核心算法升级：引入“能量门控(Energy Gating)”。
%    原理：计算 Pxx 的能量分布，生成掩膜 (Mask)。只有当某点的功率 > (最大值 - 25dB) 时，
%    才允许该点参与边界竞选。这彻底解决了“噪声区斜率变化大导致误判”的问题。
% 2. 基础参数：保持 NFFT=4096，Window=512 以维持高频域分辨率。
%--------------------------------------------------------------------------------------------------
clc
clear all
close all

%% =========================================================================
%% 信号生成部分 (保持不变)
%% =========================================================================
%% 信号1参数 (基准信号)
SubCarryNN=256;      
SubCarryN=1024;      
ratio=1/4;           
fftLen=1024;         
SymbN=60;            % 高统计量
GuardLen=SubCarryN*ratio;   
SNR=15;              
% 生成信号1
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
x=[x(1:128,1:SymbN);zeros(768,SymbN);x(129:256,1:SymbN)]; 
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
SubCarryNN2=256;      
SubCarryN2=1024;      
SymbN2=SymbN;         
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
x2=[x2(1:128,1:SymbN);zeros(768,SymbN);x2(129:256,1:SymbN)]; 
y2=ifft(x2);   
ich22=real(y2);   
qch22=imag(y2);   
ich32=[ich22(fftLen-GuardLen+1:fftLen,:);ich22];
qch32=[qch22(fftLen-GuardLen+1:fftLen,:);qch22];
ich42=reshape(ich32,1,(fftLen+GuardLen)*SymbN2);
qch42=reshape(qch32,1,(fftLen+GuardLen)*SymbN2);
TrData42=ich42+qch42.*sqrt(-1);  
ReData42=awgn(TrData42,SNR,'measured');

% 幅度缩放与频移
Amp_Scale_Sig2 = 0.6; 
ReData42 = ReData42 * Amp_Scale_Sig2;
freq_shift_ratio = 0.125; 
ReData_same = ReData42.*exp(1j * freq_shift_ratio * pi * 2 .* (1:length(ReData42)));
ReData1 = ReData; % 保存信号1（基准信号）
ReData=ReData+ReData_same; % 混合信号 

%% =========================================================================
%% 步骤1: Welch功率谱估计
%% =========================================================================
fs=1024;
Fs=1;

window_len = 512; 
window = hanning(window_len);
noverlap = window_len / 2; 
nfft_val = 4096; % 高分辨率 FFT

[Pxx1,f]=pwelch(ReData,window,noverlap,nfft_val,Fs);

% 频率轴对齐
Pxx1_shifted = fftshift(Pxx1);
Pxx1_shifted(Pxx1_shifted <= 0) = eps; 
Pxx=10*log10(Pxx1_shifted);
Pxx(isnan(Pxx) | isinf(Pxx)) = min(Pxx(isfinite(Pxx)));

% 生成对应的频率轴 [0, 1]
f_axis = linspace(0, 1, length(Pxx)); 

figure(1)
plot(f_axis, Pxx, 'b-', 'LineWidth', 1.5)
xlabel('归一化频率 (0 ~ 1)')
ylabel('功率谱密度 (dB)')
title('步骤1: 高精度功率谱估计')
grid on
xlim([0, 1]) 

%% =========================================================================
%% 步骤1.5: 分别计算两个信号的功率谱
%% =========================================================================
% 计算信号1的功率谱
[Pxx1_sig1, ~] = pwelch(ReData1, window, noverlap, nfft_val, Fs);
Pxx1_sig1_shifted = fftshift(Pxx1_sig1);
Pxx1_sig1_shifted(Pxx1_sig1_shifted <= 0) = eps;
Pxx_sig1 = 10*log10(Pxx1_sig1_shifted);
Pxx_sig1(isnan(Pxx_sig1) | isinf(Pxx_sig1)) = min(Pxx_sig1(isfinite(Pxx_sig1)));

% 计算信号2的功率谱
[Pxx1_sig2, ~] = pwelch(ReData_same, window, noverlap, nfft_val, Fs);
Pxx1_sig2_shifted = fftshift(Pxx1_sig2);
Pxx1_sig2_shifted(Pxx1_sig2_shifted <= 0) = eps;
Pxx_sig2 = 10*log10(Pxx1_sig2_shifted);
Pxx_sig2(isnan(Pxx_sig2) | isinf(Pxx_sig2)) = min(Pxx_sig2(isfinite(Pxx_sig2)));

% 绘制两个信号的频谱对比
figure(3)
plot(f_axis, Pxx_sig1, 'b-', 'LineWidth', 1.5, 'DisplayName', '信号1 (基准信号)')
hold on
plot(f_axis, Pxx_sig2, 'r-', 'LineWidth', 1.5, 'DisplayName', '信号2 (干扰信号)')
xlabel('归一化频率 (0 ~ 1)')
ylabel('功率谱密度 (dB)')
title('两个信号的频谱对比')
legend('show', 'Location', 'best')
grid on
xlim([0, 1])
hold off

%% =========================================================================
%% 步骤2: 计算功率分布函数
%% =========================================================================
Pxx_base = Pxx - min(Pxx); 
Pxx_base_sum = sum(Pxx_base);
if Pxx_base_sum == 0 || isnan(Pxx_base_sum) || isinf(Pxx_base_sum)
    error('Pxx_base求和为0或无效，请检查输入信号');
end
Pxx_int=cumsum(Pxx_base) ./ Pxx_base_sum;

y_=Pxx_int;
x_=f_axis'; 
loc=[x_,y_];

figure(2)
plot(x_, y_, 'b-', 'LineWidth', 1.5)
title('步骤2: 功率分布函数')
xlabel('归一化频率')
grid on
xlim([0, 1])

%% =========================================================================
%% 算法部分：能量门控 + 差分法 + VPD 边界检测
%% =========================================================================
% 1. 全局差分 
d_all = zeros(1, length(y_));
n1 = zeros(1, length(y_));
n2 = zeros(1, length(y_));

for i=4:length(y_)-4
    n1(i)=3*(loc(i,2))-(loc(i-1,2))-(loc(i-2,2))-(loc(i-3,2));
    n2(i)=loc(i+3,2)+loc(i+2,2)+loc(i+1,2)-3*(loc(i,2));
    d_all(i)=n2(i)-n1(i);
end

% 2. 原始特征提取
row_acc = abs(d_all);
m = length(row_acc);

%% ================= 核心修正：能量门控机制 =================
% 1. 设定门限 (动态门限)
%    从图中看，信号最高约-28dB，噪声-50dB。
%    我们取"最大值向下 20dB"作为安全区，这样-48dB以下的噪声全部被屏蔽。
Pxx_Max = max(Pxx); 
Gate_Threshold_dB = 20; % 门限深度，推荐 15~20dB
Power_Limit = Pxx_Max - Gate_Threshold_dB;

% 2. 生成二值掩膜 (Mask)
%    找出所有功率大于门限的索引位置
%    注意：Pxx的长度可能与差分数组 row_acc 长度略有不同(取决于差分算法)，
%    通常需要保证维度一致。假设 x_ 和 row_acc 是一一对应的。
Mask = zeros(size(row_acc));
for i = 1:length(row_acc)
    % 找到 row_acc(i) 对应的频率点在 Pxx 中的功率值
    % 如果您的 row_acc 是基于 y_ (即 Pxx_int) 算出来的，它们索引是一一对应的
    if Pxx(i) > Power_Limit
        Mask(i) = 1;
    else
        Mask(i) = 0; % 噪声区强制关门
    end
end

% 3. 应用掩膜：清洗特征值
%    噪声区的剧烈跳变现在全部变为了 0
row_acc_Cleaned = row_acc .* Mask; 

% =========================================================

% 均值滤波平滑 (在清洗后的数据上进行平滑)
row_acc_smooth = row_acc_Cleaned;
for k=1:2 % 增加一次平滑循环，使特征更稳健
    temp = row_acc_smooth;
    for i=2:m-1
        row_acc_smooth(i)=(temp(i-1) + temp(i)+temp(i+1))/3;
    end
end

% 最终特征值（已清洗和平滑）
row_acc_final = row_acc_smooth;

% 4. 局部极值检测 (在清洗后的 row_acc_final 上进行)
peaks = zeros(1,m); 
peakindexs = zeros(1,m); 
peakindex = 1;

% 动态阈值：仅针对有效区域
row_acc_max = max(row_acc_final);
if row_acc_max == 0
    threshold_detect = 0;
else
    % 因为噪声已被 Mask 剔除，这里的判定门限可以设得很低，以捕获微弱的重叠峰
    threshold_detect = row_acc_max * 0.05; 
end 

for i = 2:m-1
    if row_acc_final(i) > row_acc_final(i-1) && ...
       row_acc_final(i) >= row_acc_final(i+1) && ...
       row_acc_final(i) > threshold_detect
   
        peaks(peakindex)=row_acc_final(i);
        peakindexs(peakindex)=i;
        peakindex = peakindex+1;
    end
end

% 5. 选取最强4个峰值
indices = peakindexs(1:peakindex-1);
indices = indices(indices>0);

if length(indices) >= 4
    peak_values = row_acc_final(indices);
    [~, sort_idx] = sort(peak_values, 'descend');
    top4_indices = indices(sort_idx(1:4));
    boundary_indices = sort(top4_indices); 
    Detected_Locs = loc(boundary_indices(1:4), 1)';
else
    % 容错逻辑
    Detected_Locs = zeros(1,4); 
    if ~isempty(indices)
         Detected_Locs(1:length(indices)) = loc(indices, 1)';
    end
    warning('检测到的有效峰值不足4个，可能重叠过大或门限过高。');
end

figure(6)
% 归一化显示，方便观察
norm_factor = max(row_acc_final);
if norm_factor == 0, norm_factor=1; end
plot(x_, row_acc_final ./ norm_factor, 'b-', 'LineWidth', 1.5)
hold on
% 画出功率谱轮廓(灰色虚线)作为参考，直观显示为什么某些区域被 Mask 了
Pxx_norm_display = (Pxx - min(Pxx)) / (max(Pxx) - min(Pxx));
plot(x_, Pxx_norm_display, 'k:', 'LineWidth', 0.5);

if length(indices) >= 4
    plot(x_(boundary_indices), row_acc_final(boundary_indices)./norm_factor, ...
        'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r')
    for i=1:4, xline(Detected_Locs(i), '--r'); end
end
title('步骤6: 能量门控后的特征检测 (蓝色=特征, 黑色虚线=功率谱)')
legend('加权特征值', '功率谱轮廓', '检测点')
xlabel('归一化频率')
grid on
xlim([0, 1])
hold off

%% =========================================================================
%% 理论对比
%% =========================================================================
Theory_Sig1_Left  = 0.375; Theory_Sig1_Right = 0.625;
Theory_Sig2_Left  = 0.500; Theory_Sig2_Right = 0.750;
Theory_Locs = sort([Theory_Sig1_Left, Theory_Sig1_Right, Theory_Sig2_Left, Theory_Sig2_Right]);

figure(7)
plot(f_axis, Pxx, 'b-', 'LineWidth', 1.5)
hold on
% 智能调整 Y 轴显示范围
Pxx_valid = Pxx(isfinite(Pxx));
ylim_low = max(-50, min(Pxx_valid) - 5); % 下限稍微放宽，看清噪底
ylim_high = max(Pxx_valid) + 5;
ylim([ylim_low, ylim_high]); 

for i=1:4, xline(Detected_Locs(i), '--r', 'LineWidth', 2); end
for i=1:4, xline(Theory_Locs(i), '-g', 'LineWidth', 2); end
title('步骤7: 最终检测结果 (能量门控修正版)')
legend('功率谱', '检测边界', '理论边界')
xlabel('归一化频率')
grid on
xlim([0, 1]) 

fprintf('\n================== 能量门控检测结果分析 ==================\n');
fprintf('策略：能量门控机制 - 功率 > Peak-20dB 的区域保留，噪声区特征值已强制置零\n');
fprintf('------------------------------------------------------\n');
fprintf('%-10s %-12s %-12s %-12s\n', 'Idx', 'Detected', 'Theory', 'Error');
for i = 1:4
    error_val = abs(Detected_Locs(i) - Theory_Locs(i));
    fprintf('%d          %-12.4f %-12.4f %-12.4f\n', ...
        i, Detected_Locs(i), Theory_Locs(i), error_val);
end
fprintf('------------------------------------------------------\n');
mean_error = mean(abs(Detected_Locs - Theory_Locs));
fprintf('平均误差: %.5f \n', mean_error);
fprintf('======================================================\n');