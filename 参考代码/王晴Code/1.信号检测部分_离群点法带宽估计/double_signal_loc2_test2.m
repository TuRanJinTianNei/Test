% 两个信号部分重叠，高精度检测版本 (适配 .mat 文件数据版)
%-----------------------------------修改报告--------------------------------------------------
% 1. 数据源：已修改为加载 'OFDM_Rx_signal_180subcarriers_SNR15.0dB_20251225_111938.mat'。
% 2. 理论计算：根据文件中的 carrier_count(180) 和 IFFT_len(256) 动态计算理论带宽。
%    注意：新信号较宽 (70% 占比)，重叠区域会非常大。
% 3. 算法核心：保留了“能量门控”机制，确保在宽带信号下依然能准确锁定边缘。
%--------------------------------------------------------------------------------------------------
clc
clear all
close all

%% =========================================================================
%% 步骤0: 加载并构建重叠信号
%% =========================================================================
% 加载数据文件
filename = 'RX_signal.mat';
if exist(filename, 'file')
    load(filename);
    fprintf('成功加载文件: %s\n', filename);
else
    error('文件不存在，请检查路径');
end

% 提取文件中的参数
% 注意：HDF5加载后可能是变量直接存在 workspace
if ~exist('Rx_data', 'var')
    error('文件中未找到 Rx_data 变量');
end

% 确保数据为列向量或行向量统一处理 (转为行向量)
Rx_data = reshape(Rx_data, 1, []); 

% --- 构建信号 1 (基准信号) ---
% 直接使用文件中的接收数据。假设文件数据已包含噪声 (SNR=15dB)
ReData1 = Rx_data; 

% --- 构建信号 2 (干扰信号) ---
% 要求：功率幅度不同，且部分重叠
% 1. 幅度缩放 (模拟功率差异)
% 功率差异5dB: 功率比 = 10^(-5/10) ≈ 0.316, 幅度比 = sqrt(0.316) ≈ 0.562
Amp_Scale_Sig2 = 0.562; % 幅度变为 0.562倍 (功率约为 0.316倍, -5dB)
ReData2 = ReData1 * Amp_Scale_Sig2;

% 2. 频率搬移 (模拟频谱重叠)
% 由于信号很宽 (180/256), 偏移量不能太大，否则会发生混叠
% 设定偏移量为 0.125 (归一化频率)，中心从 0.5 移到 0.625
freq_shift_ratio = 0.125; 
phase_shift = exp(1j * freq_shift_ratio * pi * 2 .* (1:length(ReData2)));
ReData2_Shifted = ReData2 .* phase_shift;

% --- 合成总信号 ---
% 将两个信号叠加
ReData = ReData1 + ReData2_Shifted;

%% =========================================================================
%% 步骤1: Welch功率谱估计
%% =========================================================================
% 参数设置
fs_phy = fs; % 物理采样率 (从文件读取)
Fs_norm = 1; % 归一化采样率 (用于绘图和计算)

% Welch 参数
window_len = 512; 
window = hanning(window_len);
noverlap = window_len / 2; 
nfft_val = 4096; % 保持高分辨率

% 计算功率谱
[Pxx1, f] = pwelch(ReData, window, noverlap, nfft_val, Fs_norm);

% 频率轴对齐与归一化
Pxx1_shifted = fftshift(Pxx1);
Pxx1_shifted(Pxx1_shifted <= 0) = eps; 
Pxx = 10*log10(Pxx1_shifted);
Pxx(isnan(Pxx) | isinf(Pxx)) = min(Pxx(isfinite(Pxx)));

% 生成对应的频率轴 [0, 1]
f_axis = linspace(0, 1, length(Pxx)); 

figure(1)
plot(f_axis, Pxx, 'b-', 'LineWidth', 1.5)
xlabel('归一化频率 (0 ~ 1)')
ylabel('功率谱密度 (dB)')
title(['步骤1: 混合信号功率谱 (File: ' strrep(filename,'_','\_') ')'])
grid on
xlim([0, 1]) 

%% =========================================================================
%% 步骤1.5: 分别计算两个信号的功率谱
%% =========================================================================
% 计算信号1的功率谱
[Pxx1_sig1, ~] = pwelch(ReData1, window, noverlap, nfft_val, Fs_norm);
Pxx1_sig1_shifted = fftshift(Pxx1_sig1);
Pxx1_sig1_shifted(Pxx1_sig1_shifted <= 0) = eps;
Pxx_sig1 = 10*log10(Pxx1_sig1_shifted);
Pxx_sig1(isnan(Pxx_sig1) | isinf(Pxx_sig1)) = min(Pxx_sig1(isfinite(Pxx_sig1)));

% 计算信号2的功率谱
[Pxx1_sig2, ~] = pwelch(ReData2_Shifted, window, noverlap, nfft_val, Fs_norm);
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
if Pxx_base_sum == 0 
    warning('Pxx_base求和为0，可能是纯静默信号');
    Pxx_int = zeros(size(Pxx_base));
else
    Pxx_int = cumsum(Pxx_base) ./ Pxx_base_sum;
end

y_ = Pxx_int;
x_ = f_axis'; 
loc = [x_, y_];

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

% 3. 【能量门控 (Energy Gating)】
% 动态门限：最高峰值向下 25dB
Pxx_max = max(Pxx);
Gate_Threshold = 25; 
Pxx_Threshold = Pxx_max - Gate_Threshold; 

% 生成掩膜 (Mask)
Energy_Mask = zeros(1, m);
Energy_Mask(Pxx > Pxx_Threshold) = 1;

% 均值滤波平滑
row_acc_smooth = row_acc;
for k=1:2 
    temp = row_acc_smooth;
    for i=2:m-1
        row_acc_smooth(i)=(temp(i-1) + temp(i)+temp(i+1))/3;
    end
end

% 应用掩膜：滤除噪声区伪峰
row_acc_final = row_acc_smooth .* Energy_Mask;

% 4. 局部极值检测
peaks = zeros(1,m); 
peakindexs = zeros(1,m); 
peakindex = 1;

row_acc_max = max(row_acc_final);
if row_acc_max == 0
    threshold_detect = 0;
else
    % 门限设为最大特征值的 5%
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
    Detected_Locs = zeros(1,4); 
    if ~isempty(indices)
         Detected_Locs(1:length(indices)) = loc(indices, 1)';
    end
    warning('检测到的有效边界不足4个，可能是重叠区域没有形成明显的台阶。');
end

figure(6)
norm_factor = max(row_acc_final);
if norm_factor == 0, norm_factor=1; end
plot(x_, row_acc_final ./ norm_factor, 'b-', 'LineWidth', 1.5)
hold on
% 绘制功率谱轮廓供参考
Pxx_norm_display = (Pxx - min(Pxx)) / (max(Pxx) - min(Pxx));
plot(x_, Pxx_norm_display, 'k:', 'LineWidth', 0.5);

if length(indices) >= 4
    plot(x_(boundary_indices), row_acc_final(boundary_indices)./norm_factor, ...
        'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r')
    for i=1:4, xline(Detected_Locs(i), '--r'); end
end
title('步骤6: 能量门控后的特征检测')
legend('检测特征值', '功率谱轮廓', '检测点')
grid on
xlim([0, 1])
hold off

%% =========================================================================
%% 理论对比 (基于文件参数动态计算)
%% =========================================================================
% 从文件参数计算理论带宽
% 归一化带宽 = 有效子载波 / IFFT长度
Bandwidth_Norm = carrier_count / IFFT_bin_length; % 180/256 ≈ 0.703
Half_BW = Bandwidth_Norm / 2;

% 信号1 (基准) 中心在 0.5 (DC)
Theory_Sig1_Center = 0.5;
Theory_Sig1_Left  = Theory_Sig1_Center - Half_BW;
Theory_Sig1_Right = Theory_Sig1_Center + Half_BW;

% 信号2 (干扰) 中心在 0.5 + shift
Theory_Sig2_Center = 0.5 + freq_shift_ratio;
Theory_Sig2_Left  = Theory_Sig2_Center - Half_BW;
Theory_Sig2_Right = Theory_Sig2_Center + Half_BW;

Theory_Locs = sort([Theory_Sig1_Left, Theory_Sig1_Right, Theory_Sig2_Left, Theory_Sig2_Right]);

figure(7)
plot(f_axis, Pxx, 'b-', 'LineWidth', 1.5)
hold on
% 智能调整 Y 轴
Pxx_valid = Pxx(isfinite(Pxx));
if isempty(Pxx_valid)
    ylim([-50 0]);
else
    ylim_low = min(Pxx_valid) - 5; % 动态下限，防止信号功率过低时显示不全
    ylim_high = max(Pxx_valid) + 5;
    ylim([ylim_low, ylim_high]); 
end

for i=1:4, xline(Detected_Locs(i), '--r', 'LineWidth', 2); end
for i=1:4, xline(Theory_Locs(i), '-g', 'LineWidth', 2); end
title('步骤7: 最终检测结果 (真实数据验证)')
legend('功率谱', '检测边界', '理论边界')
xlabel('归一化频率')
grid on
xlim([0, 1]) 

fprintf('\n================== 真实信号检测报告 ==================\n');
fprintf('文件参数：Carrier=%d, IFFT=%d, Occupancy=%.1f%%\n', ...
    carrier_count, IFFT_bin_length, Bandwidth_Norm*100);
fprintf('------------------------------------------------------\n');
fprintf('%-10s %-12s %-12s %-12s\n', 'Idx', 'Detected', 'Theory', 'Error');
for i = 1:4
    if i <= length(Detected_Locs)
        error_val = abs(Detected_Locs(i) - Theory_Locs(i));
        fprintf('%d          %-12.4f %-12.4f %-12.4f\n', ...
            i, Detected_Locs(i), Theory_Locs(i), error_val);
    end
end
fprintf('------------------------------------------------------\n');