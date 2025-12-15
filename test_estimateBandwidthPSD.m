%===============================================================================
% test_estimateBandwidthPSD.m - 测试estimateBandwidthPSD.m函数
% 
% 功能说明:
%   1. 运行signalGenerate.m生成OFDM信号（发送和接收信号）
%   2. 使用estimateBandwidthPSD.m对接收信号进行功率谱密度估计
%   3. 基于PSD估计结果计算带宽（使用阈值法）
%   4. 显示估计结果和性能对比
%
% 依赖关系:
%   - signalGenerate.m：生成OFDM信号
%   - estimateBandwidthPSD.m：功率谱密度估计
%
% 使用示例:
%   test_estimateBandwidthPSD
%
% 创建日期: 2025.12.10
%===============================================================================

clc;
clear all;
close all;

fprintf('========================================\n');
fprintf('测试 estimateBandwidthPSD.m 函数\n');
fprintf('========================================\n\n');

%===============================================================================
% 步骤1: 生成OFDM信号（发送和接收）
%===============================================================================
fprintf('步骤1: 生成OFDM信号...\n');
fprintf('----------------------------------------\n');

% 运行signalGenerate.m生成信号
% 注意：需要确保signalGenerate.m已运行，或者直接调用它
try
    % 如果workspace中已有信号，直接使用
    if ~exist('Rx_data', 'var') || ~exist('Tx_data', 'var') || ~exist('fs', 'var')
        fprintf('未找到信号变量，运行 signalGenerate.m...\n');
        signalGenerate;
    else
        fprintf('使用workspace中已有的信号变量\n');
    end
    
    % 获取信号参数
    if ~exist('fs', 'var')
        fs = 15.36e6;  % 默认采样频率（与signalGenerate.m一致）
    end
    if ~exist('targetSNRdB', 'var')
        targetSNRdB = 15;  % 默认SNR
    end
    if ~exist('carrier_count', 'var')
        carrier_count = 300;  % 默认子载波数
    end
    if ~exist('subcarrier_spacing', 'var')
        subcarrier_spacing = 15e3;  % 默认子载波间隔
    end
    
    % 计算理论带宽
    B_ideal = carrier_count * subcarrier_spacing;  % 理论带宽（Hz）
    
    fprintf('信号参数：\n');
    fprintf('  - 采样频率: %.2f MHz\n', fs/1e6);
    fprintf('  - 接收信号长度: %d 样本\n', length(Rx_data));
    fprintf('  - SNR: %.1f dB\n', targetSNRdB);
    fprintf('  - 理论带宽: %.3f MHz (%.0f kHz)\n', B_ideal/1e6, B_ideal/1e3);
    fprintf('  完成！\n\n');
    
catch ME
    error('错误：无法生成或获取信号。请确保signalGenerate.m已正确运行。\n错误信息：%s', ME.message);
end

%===============================================================================
% 步骤2: 使用estimateBandwidthPSD.m进行功率谱密度估计
%===============================================================================
fprintf('步骤2: 使用estimateBandwidthPSD.m进行功率谱密度估计...\n');
fprintf('----------------------------------------\n');

% 调用功率谱密度估计函数
[Pxx_welch, f_welch, Pxx_ar, f_ar] = estimateBandwidthPSD(...
    Rx_data, fs, targetSNRdB, 'plot', true);

fprintf('功率谱密度估计完成！\n');
fprintf('  - Welch算法：频率点数 = %d\n', length(f_welch));
fprintf('  - AR模型法：频率点数 = %d\n', length(f_ar));
fprintf('  完成！\n\n');

%===============================================================================
% 步骤3: 基于PSD估计结果计算带宽（使用阈值法）
%===============================================================================
fprintf('步骤3: 基于PSD估计结果计算带宽...\n');
fprintf('----------------------------------------\n');

% Welch算法带宽计算：使用固定-3dB阈值
threshold_welch = -3;
fprintf('Welch算法：使用固定阈值 %.1f dB\n', threshold_welch);

% 找到峰值位置
[~, peak_idx_welch] = max(Pxx_welch);
peak_freq_welch = f_welch(peak_idx_welch);

% 在峰值左侧（低频侧）找阈值点
left_side_welch = Pxx_welch(1:peak_idx_welch);
left_indices_welch = find(left_side_welch <= threshold_welch);
if ~isempty(left_indices_welch)
    left_idx_welch = left_indices_welch(end);
    band_lower_welch = f_welch(left_idx_welch);
else
    band_lower_welch = f_welch(1);
    fprintf('  警告：左侧未找到阈值点，使用最低频率\n');
end

% 在峰值右侧（高频侧）找阈值点
right_side_welch = Pxx_welch(peak_idx_welch:end);
right_indices_welch = find(right_side_welch <= threshold_welch);
if ~isempty(right_indices_welch)
    right_idx_welch = peak_idx_welch + right_indices_welch(1) - 1;
    band_upper_welch = f_welch(right_idx_welch);
else
    band_upper_welch = f_welch(end);
    fprintf('  警告：右侧未找到阈值点，使用最高频率\n');
end

B_welch = abs(band_upper_welch - band_lower_welch);

% AR模型带宽计算：根据SNR自适应选择阈值
if targetSNRdB > 4
    threshold_ar = -6;
    fprintf('AR模型法：使用自适应阈值 %.1f dB (SNR > 4 dB)\n', threshold_ar);
elseif targetSNRdB > 0 && targetSNRdB <= 4
    threshold_ar = -5;
    fprintf('AR模型法：使用自适应阈值 %.1f dB (0 < SNR <= 4 dB)\n', threshold_ar);
else
    threshold_ar = -3;
    fprintf('AR模型法：使用自适应阈值 %.1f dB (SNR <= 0 dB)\n', threshold_ar);
end

% 找到峰值位置
[~, peak_idx_ar] = max(Pxx_ar);
peak_freq_ar = f_ar(peak_idx_ar);

% 在峰值左侧（低频侧）找阈值点
left_side_ar = Pxx_ar(1:peak_idx_ar);
left_indices_ar = find(left_side_ar <= threshold_ar);
if ~isempty(left_indices_ar)
    left_idx_ar = left_indices_ar(end);
    band_lower_ar = f_ar(left_idx_ar);
else
    band_lower_ar = f_ar(1);
    fprintf('  警告：左侧未找到阈值点，使用最低频率\n');
end

% 在峰值右侧（高频侧）找阈值点
right_side_ar = Pxx_ar(peak_idx_ar:end);
right_indices_ar = find(right_side_ar <= threshold_ar);
if ~isempty(right_indices_ar)
    right_idx_ar = peak_idx_ar + right_indices_ar(1) - 1;
    band_upper_ar = f_ar(right_idx_ar);
else
    band_upper_ar = f_ar(end);
    fprintf('  警告：右侧未找到阈值点，使用最高频率\n');
end

B_ar = abs(band_upper_ar - band_lower_ar);

fprintf('  完成！\n\n');

%===============================================================================
% 步骤4: 显示结果和性能分析
%===============================================================================
fprintf('========================================\n');
fprintf('带宽估计结果\n');
fprintf('========================================\n');
fprintf('理论带宽: %.3f MHz (%.0f kHz)\n', B_ideal/1e6, B_ideal/1e3);
fprintf('\n');
fprintf('Welch算法估计结果：\n');
fprintf('  - 估计带宽: %.3f MHz (%.0f kHz)\n', B_welch/1e6, B_welch/1e3);
fprintf('  - 下边界频率: %.3f MHz\n', band_lower_welch/1e6);
fprintf('  - 上边界频率: %.3f MHz\n', band_upper_welch/1e6);
fprintf('  - 峰值频率: %.3f MHz\n', peak_freq_welch/1e6);
fprintf('  - 绝对误差: %.3f MHz (%.0f kHz)\n', abs(B_welch - B_ideal)/1e6, abs(B_welch - B_ideal)/1e3);
fprintf('  - 相对误差: %.2f%%\n', abs(B_welch - B_ideal) / B_ideal * 100);
fprintf('  - 检测率: %.2f%%\n', (1 - abs(B_welch - B_ideal) / B_ideal) * 100);
fprintf('\n');
fprintf('AR模型法估计结果：\n');
fprintf('  - 估计带宽: %.3f MHz (%.0f kHz)\n', B_ar/1e6, B_ar/1e3);
fprintf('  - 下边界频率: %.3f MHz\n', band_lower_ar/1e6);
fprintf('  - 上边界频率: %.3f MHz\n', band_upper_ar/1e6);
fprintf('  - 峰值频率: %.3f MHz\n', peak_freq_ar/1e6);
fprintf('  - 绝对误差: %.3f MHz (%.0f kHz)\n', abs(B_ar - B_ideal)/1e6, abs(B_ar - B_ideal)/1e3);
fprintf('  - 相对误差: %.2f%%\n', abs(B_ar - B_ideal) / B_ideal * 100);
fprintf('  - 检测率: %.2f%%\n', (1 - abs(B_ar - B_ideal) / B_ideal) * 100);
fprintf('========================================\n\n');

%===============================================================================
% 步骤5: 绘制带宽估计结果对比图
%===============================================================================
fprintf('正在绘制带宽估计结果对比图...\n');

figure('Name', '带宽估计结果对比', 'Position', [100, 100, 1400, 800]);

% 子图1: Welch算法PSD和带宽标记
subplot(2, 2, 1);
plot(f_welch/1e6, Pxx_welch, 'b-', 'LineWidth', 1.5);
hold on;
plot([band_lower_welch/1e6, band_lower_welch/1e6], ylim, 'r--', 'LineWidth', 2, 'DisplayName', '下边界');
plot([band_upper_welch/1e6, band_upper_welch/1e6], ylim, 'r--', 'LineWidth', 2, 'DisplayName', '上边界');
plot([peak_freq_welch/1e6, peak_freq_welch/1e6], ylim, 'g--', 'LineWidth', 1.5, 'DisplayName', '峰值');
plot([band_lower_welch/1e6, band_upper_welch/1e6], [threshold_welch, threshold_welch], 'k:', 'LineWidth', 1, 'DisplayName', sprintf('阈值 (%.1f dB)', threshold_welch));
grid on;
xlabel('频率 (MHz)');
ylabel('PSD (dB)');
title(sprintf('Welch算法带宽估计\n估计带宽: %.3f MHz (理论: %.3f MHz)', B_welch/1e6, B_ideal/1e6));
legend('Location', 'best');
hold off;

% 子图2: AR模型PSD和带宽标记
subplot(2, 2, 2);
plot(f_ar/1e6, Pxx_ar, 'r-', 'LineWidth', 1.5);
hold on;
plot([band_lower_ar/1e6, band_lower_ar/1e6], ylim, 'r--', 'LineWidth', 2, 'DisplayName', '下边界');
plot([band_upper_ar/1e6, band_upper_ar/1e6], ylim, 'r--', 'LineWidth', 2, 'DisplayName', '上边界');
plot([peak_freq_ar/1e6, peak_freq_ar/1e6], ylim, 'g--', 'LineWidth', 1.5, 'DisplayName', '峰值');
plot([band_lower_ar/1e6, band_upper_ar/1e6], [threshold_ar, threshold_ar], 'k:', 'LineWidth', 1, 'DisplayName', sprintf('阈值 (%.1f dB)', threshold_ar));
grid on;
xlabel('频率 (MHz)');
ylabel('PSD (dB)');
title(sprintf('AR模型法带宽估计\n估计带宽: %.3f MHz (理论: %.3f MHz)', B_ar/1e6, B_ideal/1e6));
legend('Location', 'best');
hold off;

% 子图3: 两种方法PSD对比
subplot(2, 2, 3);
plot(f_welch/1e6, Pxx_welch, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Welch算法');
hold on;
% 需要插值AR模型的PSD到Welch的频率点进行对比，或者分别绘制
% 这里简单绘制，注意频率范围可能不同
f_ar_limited = f_ar(f_ar <= max(f_welch));
Pxx_ar_limited = Pxx_ar(1:length(f_ar_limited));
if ~isempty(f_ar_limited)
    plot(f_ar_limited/1e6, Pxx_ar_limited, 'r-', 'LineWidth', 1.5, 'DisplayName', 'AR模型法');
end
grid on;
xlabel('频率 (MHz)');
ylabel('PSD (dB)');
title('功率谱密度估计对比');
legend('Location', 'best');
hold off;

% 子图4: 带宽估计结果柱状图
subplot(2, 2, 4);
methods = {'Welch算法', 'AR模型法', '理论值'};
bandwidths = [B_welch/1e6, B_ar/1e6, B_ideal/1e6];
colors = [0.2 0.6 0.8; 0.8 0.2 0.2; 0.2 0.8 0.2];
b = bar(bandwidths, 'FaceColor', 'flat');
b.CData = colors;
set(gca, 'XTickLabel', methods);
ylabel('带宽 (MHz)');
title('带宽估计结果对比');
grid on;
hold on;
% 添加数值标签
for i = 1:length(bandwidths)
    text(i, bandwidths(i) + 0.05, sprintf('%.3f MHz', bandwidths(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', 10);
end
hold off;

fprintf('图形绘制完成！\n\n');

fprintf('========================================\n');
fprintf('测试完成！\n');
fprintf('========================================\n');
fprintf('所有结果已显示在上方的图形窗口中\n');
fprintf('变量说明：\n');
fprintf('  - B_welch: Welch算法估计的带宽 (Hz)\n');
fprintf('  - B_ar: AR模型法估计的带宽 (Hz)\n');
fprintf('  - B_ideal: 理论带宽 (Hz)\n');
fprintf('  - Pxx_welch, f_welch: Welch算法PSD估计结果\n');
fprintf('  - Pxx_ar, f_ar: AR模型法PSD估计结果\n');
fprintf('========================================\n');

