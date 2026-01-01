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
% 从两个独立的mat文件中加载信号
filename1 = '1.mat';
filename2 = '2.mat';

% 加载信号1
if exist(filename1, 'file')
    load(filename1);
    fprintf('成功加载文件: %s\n', filename1);
else
    error('文件不存在，请检查路径: %s', filename1);
end

% 提取信号1数据
if exist('Rx_data', 'var')
    Rx_data1 = Rx_data;
else
    % 尝试其他可能的变量名
    vars = who('-file', filename1);
    if isempty(vars)
        error('文件 %s 中未找到数据变量', filename1);
    else
        % 使用第一个变量
        temp = load(filename1);
        Rx_data1 = temp.(vars{1});
        fprintf('使用变量 %s 作为信号1数据\n', vars{1});
    end
end

% 加载信号2
if exist(filename2, 'file')
    load(filename2);
    fprintf('成功加载文件: %s\n', filename2);
else
    error('文件不存在，请检查路径: %s', filename2);
end

% 提取信号2数据
if exist('Rx_data', 'var')
    Rx_data2 = Rx_data;
else
    % 尝试其他可能的变量名
    vars = who('-file', filename2);
    if isempty(vars)
        error('文件 %s 中未找到数据变量', filename2);
    else
        % 使用第一个变量
        temp = load(filename2);
        Rx_data2 = temp.(vars{1});
        fprintf('使用变量 %s 作为信号2数据\n', vars{1});
    end
end

% 确保数据为行向量统一处理
ReData1 = reshape(Rx_data1, 1, []); 
ReData2 = reshape(Rx_data2, 1, []); 

% 检查两个信号长度是否一致，如果不一致则截断到较短的长度
min_len = min(length(ReData1), length(ReData2));
if length(ReData1) ~= length(ReData2)
    warning('两个信号长度不一致 (信号1: %d, 信号2: %d)，将截断到较短长度: %d', ...
        length(ReData1), length(ReData2), min_len);
    ReData1 = ReData1(1:min_len);
    ReData2 = ReData2(1:min_len);
end

% --- 合成总信号 ---
% 将两个信号直接叠加
ReData = ReData1 + ReData2;

%% =========================================================================
%% 步骤1: Welch功率谱估计
%% =========================================================================
% 参数设置
% 尝试从文件中获取采样率，如果不存在则使用默认值
try
    file1_vars = load(filename1);
    if isfield(file1_vars, 'fs')
        fs_phy = file1_vars.fs;
    else
        fs_phy = 1024; % 默认采样率
        warning('未找到采样率参数fs，使用默认值1024');
    end
catch
    fs_phy = 1024; % 默认采样率
    warning('无法读取文件参数，使用默认采样率1024');
end
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
title(['步骤1: 混合信号功率谱 (From: ' strrep(filename1,'_','\_') ' + ' strrep(filename2,'_','\_') ')'])
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
[Pxx1_sig2, ~] = pwelch(ReData2, window, noverlap, nfft_val, Fs_norm);
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
%% 步骤1.6: 循环谱分析
%% =========================================================================
fprintf('计算循环谱...\n');

% 循环谱参数设置
N_fft = 512;              % FFT点数（频率分辨率）
N_alpha = 128;            % 循环频率分辨率
N_data = length(ReData);

% 计算循环频率范围（归一化）
% 循环频率通常与符号速率相关，对于OFDM信号，循环频率在符号速率附近
alpha_max = 0.05;  % 最大循环频率（归一化）
alpha_vec = linspace(-alpha_max, alpha_max, N_alpha);

% 频率范围（归一化）
f_vec = linspace(-0.5, 0.5, N_fft);

% 初始化循环谱矩阵
S_alpha_f = zeros(N_alpha, N_fft);

% 使用FFT累积量方法（FAM）计算循环谱
% 对于每个循环频率α，计算循环自相关函数，然后做FFT
fprintf('  正在计算循环谱，进度: ');
for alpha_idx = 1:N_alpha
    if mod(alpha_idx, 20) == 0 || alpha_idx == 1
        fprintf('%d/%d ', alpha_idx, N_alpha);
    end
    
    alpha = alpha_vec(alpha_idx);
    
    % 计算循环自相关函数 R_x^α(τ)
    % R_x^α(τ) = E[x(t+τ/2) x*(t-τ/2) e^{-j2παt}]
    tau_max = min(256, floor(N_data/4));
    tau_vec = -tau_max:tau_max;
    R_alpha = zeros(size(tau_vec));
    
    % 对每个时延τ计算循环自相关
    for tau_idx = 1:length(tau_vec)
        tau = tau_vec(tau_idx);
        
        % 计算有效数据范围
        t_start = max(1, abs(tau) + 1);
        t_end = min(N_data, N_data - abs(tau));
        
        if t_end > t_start
            % 提取数据段
            if tau >= 0
                x1 = ReData(t_start:t_end);
                x2 = ReData(t_start + tau:t_end + tau);
            else
                x1 = ReData(t_start - tau:t_end - tau);
                x2 = ReData(t_start:t_end);
            end
            
            % 计算循环自相关（使用共轭）
            t_indices = (t_start:t_end)';
            phase_term = exp(-1j * 2 * pi * alpha * t_indices);
            R_alpha(tau_idx) = mean(x1(:) .* conj(x2(:)) .* phase_term(:));
        end
    end
    
    % 对循环自相关函数做FFT得到循环谱
    % 使用零填充以提高频率分辨率
    if length(R_alpha) < N_fft
        R_alpha_padded = [zeros(1, floor((N_fft - length(R_alpha))/2)), ...
                          R_alpha, ...
                          zeros(1, ceil((N_fft - length(R_alpha))/2))];
    else
        R_alpha_padded = R_alpha(1:N_fft);
    end
    
    % FFT并移位
    S_alpha_f(alpha_idx, :) = fftshift(fft(R_alpha_padded, N_fft));
end
fprintf('\n');

% 转换为dB
S_alpha_f_abs = abs(S_alpha_f);
S_alpha_f_max = max(S_alpha_f_abs(:));
S_alpha_f_dB = 10*log10(S_alpha_f_abs + eps * S_alpha_f_max);

% 绘制循环谱（3D表面图）
figure(8)
surf(f_vec, alpha_vec, S_alpha_f_dB, 'EdgeColor', 'none');
view(2);  % 俯视图
colorbar;
xlabel('归一化频率 f', 'FontSize', 12);
ylabel('循环频率 α', 'FontSize', 12);
title('循环谱 (Cyclic Spectrum)', 'FontSize', 14);
colormap('jet');
shading interp;
axis tight;
grid on;

% 绘制循环谱（2D等高线图）
figure(9)
contourf(f_vec, alpha_vec, S_alpha_f_dB, 30);
colorbar;
xlabel('归一化频率 f', 'FontSize', 12);
ylabel('循环频率 α', 'FontSize', 12);
title('循环谱等高线图 (Cyclic Spectrum Contour)', 'FontSize', 14);
colormap('jet');
grid on;

% 绘制循环频率切片（在α=0处的切片，即传统功率谱）
figure(10)
alpha_zero_idx = find(abs(alpha_vec) == min(abs(alpha_vec)), 1);
plot(f_vec, S_alpha_f_dB(alpha_zero_idx, :), 'b-', 'LineWidth', 1.5);
hold on;
% 对比传统功率谱（需要插值到相同的频率点）
Pxx_interp = interp1(f_axis, Pxx, f_vec + 0.5, 'linear', 'extrap');
Pxx_interp_norm = Pxx_interp - max(Pxx_interp);
plot(f_vec, Pxx_interp_norm, 'r--', 'LineWidth', 1.5);
xlabel('归一化频率', 'FontSize', 12);
ylabel('功率 (dB)', 'FontSize', 12);
title('循环谱在α=0处的切片 vs 传统功率谱', 'FontSize', 14);
legend('循环谱 (α=0)', '传统功率谱', 'Location', 'best');
grid on;
hold off;

% 绘制循环频率剖面（在f=0处的剖面）
figure(11)
f_zero_idx = find(abs(f_vec) == min(abs(f_vec)), 1);
plot(alpha_vec, S_alpha_f_dB(:, f_zero_idx), 'b-', 'LineWidth', 1.5);
xlabel('循环频率 α', 'FontSize', 12);
ylabel('功率 (dB)', 'FontSize', 12);
title('循环谱在f=0处的剖面', 'FontSize', 14);
grid on;

% 绘制循环谱的峰值检测（寻找循环频率特征）
figure(12)
% 对每个频率f，找到最大循环谱值
[max_val_per_f, max_alpha_idx_per_f] = max(S_alpha_f_dB, [], 1);
plot(f_vec, max_val_per_f, 'b-', 'LineWidth', 1.5);
hold on;
% 标记峰值位置
[peaks, peak_locs] = findpeaks(max_val_per_f, 'MinPeakHeight', max(max_val_per_f) - 10);
plot(f_vec(peak_locs), peaks, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
xlabel('归一化频率 f', 'FontSize', 12);
ylabel('最大循环谱功率 (dB)', 'FontSize', 12);
title('循环谱峰值检测（各频率处的最大循环谱值）', 'FontSize', 14);
legend('最大循环谱值', '峰值位置', 'Location', 'best');
grid on;
hold off;

fprintf('循环谱计算完成！\n');
fprintf('  循环频率范围: [%.4f, %.4f]\n', min(alpha_vec), max(alpha_vec));
fprintf('  频率分辨率: %.4f\n', f_vec(2) - f_vec(1));
fprintf('  循环频率分辨率: %.6f\n', alpha_vec(2) - alpha_vec(1));
fprintf('\n');

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
% 尝试从文件参数计算理论带宽
% 如果文件中有相关参数则计算，否则跳过理论对比
calc_theory = false;
Theory_Locs = [];

% 重新加载文件以获取参数（如果存在）
try
    % 加载文件1的参数
    file1_vars = load(filename1);
    % 加载文件2的参数
    file2_vars = load(filename2);
    
    % 检查是否存在必要的参数
    if isfield(file1_vars, 'carrier_count') && isfield(file1_vars, 'IFFT_bin_length')
        carrier_count1 = file1_vars.carrier_count;
        IFFT_bin_length1 = file1_vars.IFFT_bin_length;
        Bandwidth_Norm1 = carrier_count1 / IFFT_bin_length1;
        Half_BW1 = Bandwidth_Norm1 / 2;
        
        % 检查信号1的中心频率（如果文件中有）
        if isfield(file1_vars, 'center_freq_norm')
            Theory_Sig1_Center = file1_vars.center_freq_norm;
        else
            Theory_Sig1_Center = 0.5; % 默认中心频率
        end
        Theory_Sig1_Left  = Theory_Sig1_Center - Half_BW1;
        Theory_Sig1_Right = Theory_Sig1_Center + Half_BW1;
        
        % 检查信号2的参数
        if isfield(file2_vars, 'carrier_count') && isfield(file2_vars, 'IFFT_bin_length')
            carrier_count2 = file2_vars.carrier_count;
            IFFT_bin_length2 = file2_vars.IFFT_bin_length;
            Bandwidth_Norm2 = carrier_count2 / IFFT_bin_length2;
            Half_BW2 = Bandwidth_Norm2 / 2;
            
            % 检查信号2的中心频率
            if isfield(file2_vars, 'center_freq_norm')
                Theory_Sig2_Center = file2_vars.center_freq_norm;
            else
                Theory_Sig2_Center = 0.5; % 默认中心频率
            end
            Theory_Sig2_Left  = Theory_Sig2_Center - Half_BW2;
            Theory_Sig2_Right = Theory_Sig2_Center + Half_BW2;
            
            Theory_Locs = sort([Theory_Sig1_Left, Theory_Sig1_Right, Theory_Sig2_Left, Theory_Sig2_Right]);
            calc_theory = true;
        end
    end
catch
    warning('无法从文件中读取理论参数，将跳过理论对比');
end

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
if calc_theory && ~isempty(Theory_Locs)
    for i=1:4, xline(Theory_Locs(i), '-g', 'LineWidth', 2); end
    title('步骤7: 最终检测结果 (真实数据验证)')
    legend('功率谱', '检测边界', '理论边界')
else
    title('步骤7: 最终检测结果 (真实数据验证)')
    legend('功率谱', '检测边界')
end
xlabel('归一化频率')
grid on
xlim([0, 1]) 

fprintf('\n================== 真实信号检测报告 ==================\n');
fprintf('数据来源：%s + %s\n', filename1, filename2);
fprintf('信号1长度：%d, 信号2长度：%d, 混合信号长度：%d\n', ...
    length(ReData1), length(ReData2), length(ReData));
if calc_theory
    fprintf('文件参数：信号1 (Carrier=%d, IFFT=%d, Occupancy=%.1f%%), ', ...
        carrier_count1, IFFT_bin_length1, Bandwidth_Norm1*100);
    fprintf('信号2 (Carrier=%d, IFFT=%d, Occupancy=%.1f%%)\n', ...
        carrier_count2, IFFT_bin_length2, Bandwidth_Norm2*100);
end
fprintf('------------------------------------------------------\n');
fprintf('%-10s %-12s', 'Idx', 'Detected');
if calc_theory
    fprintf(' %-12s %-12s\n', 'Theory', 'Error');
else
    fprintf('\n');
end
for i = 1:4
    if i <= length(Detected_Locs)
        if calc_theory && i <= length(Theory_Locs)
            error_val = abs(Detected_Locs(i) - Theory_Locs(i));
            fprintf('%d          %-12.4f %-12.4f %-12.4f\n', ...
                i, Detected_Locs(i), Theory_Locs(i), error_val);
        else
            fprintf('%d          %-12.4f\n', i, Detected_Locs(i));
        end
    end
end
fprintf('------------------------------------------------------\n');