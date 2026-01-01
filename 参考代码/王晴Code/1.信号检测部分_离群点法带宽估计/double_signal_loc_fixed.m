% double_signal_loc_fixed.m
% 功能：双信号重叠盲检测（修复衰落信道误判版）
% 核心改进：引入形态学闭运算填平频率选择性衰落的深坑
% =========================================================================

clc;
clear all;
close all;

%% =========================================================================
%% 步骤0: 加载数据与构建测试场景
%% =========================================================================
filename = 'RX_signal.mat';
if exist(filename, 'file')
    load(filename);
    fprintf('成功加载文件: %s\n', filename);
else
    error('文件 RX_signal.mat 不存在，请先运行生成器生成数据。');
end

% 确保数据维度正确
Rx_data = reshape(Rx_data, 1, []);

% --- 构建部分重叠的双信号场景 ---
% 信号1: 原始接收信号 (基准)
Sig1 = Rx_data;

% 信号2: 模拟第二个用户 (功率较低，频率偏移)
Amp_Scale = 0.6;          % 功率约为信号1的 36%
Freq_Shift_Ratio = 0.125; % 频率偏移量 (归一化频率 0~1)

% 构造频移相位因子
t_vec = 0:(length(Sig1)-1);
phase_shift = exp(1j * 2 * pi * Freq_Shift_Ratio * t_vec);
Sig2 = (Sig1 .* phase_shift) * Amp_Scale;

% 合成重叠信号
ReData = Sig1 + Sig2;

%% =========================================================================
%% 步骤1: 功率谱估计 (PSD)
%% =========================================================================
% 使用 Welch 法估计功率谱
NFFT = 4096;              % 高分辨率 FFT
Window_Len = 512;
Window = hanning(Window_Len);
N_Overlap = Window_Len / 2;

[Pxx_raw, f] = pwelch(ReData, Window, N_Overlap, NFFT, 'centered');

% 归一化频率轴 [0, 1] 用于显示和计算
f_norm = linspace(0, 1, length(Pxx_raw));

% 转为对数谱 (dB)
Pxx_dB = 10 * log10(Pxx_raw + eps); 
% 简单底噪钳位，防止负无穷干扰
min_dB = max(Pxx_dB) - 60; 
Pxx_dB(Pxx_dB < min_dB) = min_dB;

%% =========================================================================
%% 步骤1.5: 【关键修复】形态学滤波处理
%% =========================================================================
% 问题：瑞利衰落导致频谱出现深坑，二阶差分算法会将其误判为边缘。
% 解决：使用形态学闭运算 (Closing) = 先膨胀 (Dilation) 后腐蚀 (Erosion)。
% 作用：填平深坑，但保持外包络边缘位置不变。

% 1. 定义结构元素大小 (Structure Element Size)
% 盲估计参数：取 NFFT 的 1.5% 左右，足以覆盖衰落坑，但小于信号带宽。
SE_Size = round(NFFT * 0.015); % 4096 * 0.015 ≈ 61 个点

% 2. 膨胀 (Dilation) - 使用移动最大值填充低谷
Pxx_Dilated = movmax(Pxx_dB, SE_Size);

% 3. 腐蚀 (Erosion) - 使用移动最小值恢复轮廓
Pxx_Closed = movmin(Pxx_Dilated, SE_Size);

% 4. 转换回线性域用于计算 CDF (因为 CDF 基于能量累积)
Pxx_Processed_Lin = 10.^(Pxx_Closed / 10);

% --- 绘图验证处理效果 ---
figure('Name', '形态学滤波效果对比', 'Position', [100, 100, 1000, 400]);
plot(f_norm, Pxx_dB, 'Color', [0.7 0.7 0.7], 'LineWidth', 1); hold on;
plot(f_norm, Pxx_Closed, 'r-', 'LineWidth', 1.5);
legend('原始衰落频谱 (含深坑)', '处理后频谱 (闭运算包络)');
title('预处理：填平衰落深坑以防止误判');
grid on;
xlim([0, 1]);

%% =========================================================================
%% 步骤2: 计算功率分布函数 (CDF) 与 特征提取
%% =========================================================================
% 使用处理后的线性谱计算 CDF
Pxx_Base = Pxx_Processed_Lin - min(Pxx_Processed_Lin);
Pxx_Sum = sum(Pxx_Base);
if Pxx_Sum == 0, Pxx_Sum = eps; end

% 计算累积分布 (CDF)
Pxx_CDF = cumsum(Pxx_Base) / Pxx_Sum;
y = Pxx_CDF; 
m = length(y);

% 计算二阶差分特征 (模拟曲率)
% 算法：d(i) = [Sum(右3点) - 3*Center] - [3*Center - Sum(左3点)]
d_all = zeros(1, m);
for i = 4 : m-4
    n1 = 3*y(i) - y(i-1) - y(i-2) - y(i-3);
    n2 = y(i+1) + y(i+2) + y(i+3) - 3*y(i);
    d_all(i) = n2 - n1;
end
feature_curve = abs(d_all);

%% =========================================================================
%% 步骤3: 能量门控与峰值检测
%% =========================================================================
% 1. 能量门控 (Energy Gating)
% 目的：只在信号能量存在的区域检测边缘，去除纯噪声区的抖动
Gate_Threshold_dB = 25; % 相对于峰值下降 25dB
Mask_Threshold = max(Pxx_Closed) - Gate_Threshold_dB;
Energy_Mask = (Pxx_Closed > Mask_Threshold)'; % 生成掩膜 (转置为行向量匹配 feature_curve)

% 应用掩膜
feature_masked = feature_curve .* Energy_Mask;

% 2. 峰值搜索
% 门限：最大特征值的 10% (提高鲁棒性)
Peak_Threshold = max(feature_masked) * 0.10; 
[pks, locs] = findpeaks(feature_masked, 'MinPeakHeight', Peak_Threshold, 'MinPeakDistance', SE_Size);

% 3. 提取前4个最显著的峰 (假设双信号重叠最多有4个边界)
% 排序：按峰值高度降序
[~, sort_idx] = sort(pks, 'descend');
num_peaks_to_keep = min(4, length(sort_idx));
top_indices = locs(sort_idx(1:num_peaks_to_keep));

% 按频率位置重新排序 (从左到右)
Detected_Indices = sort(top_indices);
Detected_Freqs = f_norm(Detected_Indices);

%% =========================================================================
%% 步骤4: 结果验证与可视化
%% =========================================================================
% 计算理论值 (基于 .mat 文件中的参数)
% 归一化带宽 = 有效子载波 / IFFT长度
BW_Norm = carrier_count / IFFT_bin_length; 
Half_BW = BW_Norm / 2;

% 信号1 (基准) 中心在 0.5
Theory_Locs_Sig1 = [0.5 - Half_BW, 0.5 + Half_BW];
% 信号2 (频移) 中心在 0.5 + Offset
Theory_Locs_Sig2 = [0.5 + Freq_Shift_Ratio - Half_BW, 0.5 + Freq_Shift_Ratio + Half_BW];

Theory_Locs = sort([Theory_Locs_Sig1, Theory_Locs_Sig2]);

% --- 绘制最终检测结果 ---
figure('Name', '最终检测结果', 'Position', [100, 550, 1000, 500]);
yyaxis left
plot(f_norm, Pxx_dB, 'Color', [0.8 0.8 1.0], 'LineWidth', 1); hold on;
plot(f_norm, Pxx_Closed, 'b-', 'LineWidth', 1.5);
ylabel('功率谱密度 (dB)');
ylim([min(Pxx_Closed)-10, max(Pxx_Closed)+5]);

yyaxis right
% 归一化特征曲线以便显示
feat_norm = feature_masked / max(feature_masked);
plot(f_norm, feat_norm, 'm-', 'LineWidth', 1.5);
ylabel('归一化特征值');
ylim([0, 1.2]);

% 标记检测点
for i = 1:length(Detected_Freqs)
    xline(Detected_Freqs(i), '--r', 'LineWidth', 1.5);
    text(Detected_Freqs(i), 1.05, sprintf('P%d', i), 'Color', 'r', 'HorizontalAlignment', 'center');
end

title('双信号重叠检测 (盲估计 + 形态学修复)');
xlabel('归一化频率');
legend('原始谱', '修复谱', '检测特征曲线', '检测边界');
grid on;
xlim([0, 1]);

% --- 命令行输出报告 ---
fprintf('\n================== 真实信号检测报告 ==================\n');
fprintf('文件参数：Carrier=%d, IFFT=%d, Occupancy=%.1f%%\n', ...
    carrier_count, IFFT_bin_length, BW_Norm*100);
fprintf('------------------------------------------------------\n');
fprintf('%-10s %-12s %-12s %-12s\n', 'Idx', 'Detected', 'Theory', 'Error');

for i = 1:4
    d_val = 0; t_val = 0; err_val = NaN;
    
    if i <= length(Detected_Freqs)
        d_val = Detected_Freqs(i);
    else
        d_val = NaN;
    end
    
    if i <= length(Theory_Locs)
        t_val = Theory_Locs(i);
    end
    
    if ~isnan(d_val)
        err_val = abs(d_val - t_val);
        fprintf('%d          %-12.4f %-12.4f %-12.4f\n', i, d_val, t_val, err_val);
    else
        fprintf('%d          %-12s %-12.4f %-12s\n', i, 'Missed', t_val, '---');
    end
end
fprintf('------------------------------------------------------\n');
fprintf('注: 误差 < 0.01 (约1%%) 通常被认为检测精确。\n');