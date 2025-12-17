%===============================================================================
% example_bandwidthEstimation.m - 带宽估计方法示例
% 
% 功能说明:
%   演示如何使用 estimateBandwidthFromPSD.m 从功率谱密度估计带宽
%   展示多种估计方法的使用和对比
%
% 依赖关系:
%   - PSD_SNR15.0dB_300subcarriers_20251216_111043.mat：保存的PSD数据文件（优先使用）
%   - signalGenerate.m：生成OFDM信号（备用）
%   - estimateBandwidthPSD.m：功率谱密度估计（备用）
%   - estimateBandwidthFromPSD.m：带宽估计（新函数）
%
% 使用示例:
%   example_bandwidthEstimation
%
% 创建日期: 2025.12.10
% 修改日期: 2025.12.16 - 添加从.mat文件加载PSD数据的功能
%===============================================================================

clc;
clear all;
close all;

fprintf('========================================\n');
fprintf('带宽估计方法示例\n');
fprintf('========================================\n\n');

%===============================================================================
% 步骤1: 初始化参数（备用，如果从.mat文件加载PSD则参数会从文件中获取）
%===============================================================================
% 功能：设置默认参数，如果从.mat文件加载PSD，这些参数会被文件中的参数覆盖
fprintf('步骤1: 初始化参数...\n');
fprintf('----------------------------------------\n');

% 设置默认参数（如果变量不存在）
if ~exist('fs', 'var')
    fs = 15.36e6;  % 默认采样频率：15.36 MHz
end
if ~exist('targetSNRdB', 'var')
    targetSNRdB = 15;  % 默认SNR：15 dB
end
if ~exist('carrier_count', 'var')
    carrier_count = 300;  % 默认子载波数：300
end
if ~exist('subcarrier_spacing', 'var')
    subcarrier_spacing = 15e3;  % 默认子载波间隔：15 kHz
end

% 计算理论带宽（如果从.mat文件加载，会被文件中的B_ideal覆盖）
B_ideal = carrier_count * subcarrier_spacing;

fprintf('默认参数：\n');
fprintf('  - 采样频率: %.2f MHz\n', fs/1e6);
fprintf('  - SNR: %.1f dB\n', targetSNRdB);
fprintf('  - 子载波数: %d\n', carrier_count);
fprintf('  - 子载波间隔: %.1f kHz\n', subcarrier_spacing/1e3);
fprintf('  - 理论带宽: %.3f MHz (%.0f kHz)\n', B_ideal/1e6, B_ideal/1e3);
fprintf('  注意：如果从.mat文件加载PSD，参数将从文件中读取\n');
fprintf('  完成！\n\n');

%===============================================================================
% 步骤2: 加载或计算功率谱密度
%===============================================================================
% 功能：优先从.mat文件加载已保存的PSD数据，如果文件不存在则重新计算
fprintf('步骤2: 加载或计算功率谱密度...\n');
fprintf('----------------------------------------\n');

psd_filename = 'PSD_SNR15.0dB_300subcarriers_20251216_111043.mat';

try
    % 优先从.mat文件加载PSD数据
    if exist(psd_filename, 'file')
        fprintf('从文件加载PSD数据: %s\n', psd_filename);
        loaded_psd = load(psd_filename);
        
        % 检查文件中是否包含必要的PSD变量
        if isfield(loaded_psd, 'Pxx_welch') && isfield(loaded_psd, 'f_welch') && ...
           isfield(loaded_psd, 'Pxx_ar') && isfield(loaded_psd, 'f_ar')
            Pxx_welch = loaded_psd.Pxx_welch;
            f_welch = loaded_psd.f_welch;
            Pxx_ar = loaded_psd.Pxx_ar;
            f_ar = loaded_psd.f_ar;
            
            fprintf('  成功加载 Pxx_welch, f_welch (Welch算法PSD)\n');
            fprintf('  成功加载 Pxx_ar, f_ar (AR模型法PSD)\n');
            
            % 如果文件中有参数，使用文件中的参数
            if isfield(loaded_psd, 'fs')
                fs = loaded_psd.fs;
            end
            if isfield(loaded_psd, 'targetSNRdB')
                targetSNRdB = loaded_psd.targetSNRdB;
            end
            if isfield(loaded_psd, 'carrier_count')
                carrier_count = loaded_psd.carrier_count;
            end
            if isfield(loaded_psd, 'subcarrier_spacing')
                subcarrier_spacing = loaded_psd.subcarrier_spacing;
            end
            if isfield(loaded_psd, 'B_ideal')
                B_ideal = loaded_psd.B_ideal;
            else
                B_ideal = carrier_count * subcarrier_spacing;
            end
            
            fprintf('  使用文件中的参数：SNR=%.1f dB, 子载波数=%d\n', ...
                targetSNRdB, carrier_count);
            
        else
            error('文件中缺少必要的PSD变量（Pxx_welch, f_welch, Pxx_ar, f_ar）');
        end
        
    % 如果文件不存在，重新计算PSD
    else
        fprintf('未找到PSD数据文件，重新计算功率谱密度...\n');
        fprintf('  文件路径: %s\n', psd_filename);
        
        % 确保有接收信号数据
        if ~exist('Rx_data', 'var')
            fprintf('  未找到Rx_data，运行signalGenerate.m生成信号...\n');
            signalGenerate;
        end
        
        % 调用PSD估计函数
        [Pxx_welch, f_welch, Pxx_ar, f_ar] = estimateBandwidthPSD(...
            Rx_data, fs, targetSNRdB, 'plot', false);
        
        fprintf('  功率谱密度计算完成！\n');
    end
    
    fprintf('PSD数据准备完成！\n');
    fprintf('  - Welch算法：频率点数 = %d\n', length(f_welch));
    fprintf('  - AR模型法：频率点数 = %d\n', length(f_ar));
    fprintf('  完成！\n\n');
    
catch ME
    error('错误：无法加载或计算PSD数据。\n错误信息：%s\n\n请确保：\n  1. .mat文件存在且包含PSD变量\n  2. 或者可以运行signalGenerate.m和estimateBandwidthPSD.m', ME.message);
end

%===============================================================================
% 步骤3: 使用多种方法估计带宽
%===============================================================================
fprintf('步骤3: 使用多种方法估计带宽...\n');
fprintf('----------------------------------------\n');

% 方法1：阈值法（固定阈值-3dB）
fprintf('\n方法1：阈值法（固定阈值-3dB）\n');
[B_threshold_welch, results_threshold_welch] = estimateBandwidthFromPSD(...
    Pxx_welch, f_welch, 'method', 'threshold', 'threshold', -3, 'plot', false);
fprintf('  Welch算法：估计带宽 = %.3f MHz\n', B_threshold_welch/1e6);
fprintf('    下边界: %.3f MHz, 上边界: %.3f MHz\n', ...
    results_threshold_welch.lower_freq/1e6, results_threshold_welch.upper_freq/1e6);

[B_threshold_ar, results_threshold_ar] = estimateBandwidthFromPSD(...
    Pxx_ar, f_ar, 'method', 'threshold', 'threshold', -3, 'plot', false);
fprintf('  AR模型法：估计带宽 = %.3f MHz\n', B_threshold_ar/1e6);
fprintf('    下边界: %.3f MHz, 上边界: %.3f MHz\n', ...
    results_threshold_ar.lower_freq/1e6, results_threshold_ar.upper_freq/1e6);

% 方法2：能量百分比法（90%能量带宽）
fprintf('\n方法2：能量百分比法（90%%能量带宽）\n');
[B_energy_welch, results_energy_welch] = estimateBandwidthFromPSD(...
    Pxx_welch, f_welch, 'method', 'energy', 'energy_percent', 0.9, 'plot', false);
fprintf('  Welch算法：估计带宽 = %.3f MHz\n', B_energy_welch/1e6);
fprintf('    实际能量占比: %.2f%%\n', results_energy_welch.actual_energy_percent*100);

[B_energy_ar, results_energy_ar] = estimateBandwidthFromPSD(...
    Pxx_ar, f_ar, 'method', 'energy', 'energy_percent', 0.9, 'plot', false);
fprintf('  AR模型法：估计带宽 = %.3f MHz\n', B_energy_ar/1e6);
fprintf('    实际能量占比: %.2f%%\n', results_energy_ar.actual_energy_percent*100);

% 方法3：峰值下降法（从峰值下降3dB）
fprintf('\n方法3：峰值下降法（从峰值下降3dB）\n');
[B_peak_drop_welch, results_peak_drop_welch] = estimateBandwidthFromPSD(...
    Pxx_welch, f_welch, 'method', 'peak_drop', 'peak_drop_db', 3, 'plot', false);
fprintf('  Welch算法：估计带宽 = %.3f MHz\n', B_peak_drop_welch/1e6);
fprintf('    峰值功率: %.2f dB, 阈值: %.2f dB\n', ...
    results_peak_drop_welch.peak_power, results_peak_drop_welch.threshold);

[B_peak_drop_ar, results_peak_drop_ar] = estimateBandwidthFromPSD(...
    Pxx_ar, f_ar, 'method', 'peak_drop', 'peak_drop_db', 3, 'plot', false);
fprintf('  AR模型法：估计带宽 = %.3f MHz\n', B_peak_drop_ar/1e6);
fprintf('    峰值功率: %.2f dB, 阈值: %.2f dB\n', ...
    results_peak_drop_ar.peak_power, results_peak_drop_ar.threshold);

% 方法4：RMS带宽法
fprintf('\n方法4：RMS带宽法（二阶矩法）\n');
[B_rms_welch, results_rms_welch] = estimateBandwidthFromPSD(...
    Pxx_welch, f_welch, 'method', 'rms', 'plot', false);
fprintf('  Welch算法：估计带宽 = %.3f MHz\n', B_rms_welch/1e6);
fprintf('    中心频率: %.3f MHz\n', results_rms_welch.center_freq/1e6);

[B_rms_ar, results_rms_ar] = estimateBandwidthFromPSD(...
    Pxx_ar, f_ar, 'method', 'rms', 'plot', false);
fprintf('  AR模型法：估计带宽 = %.3f MHz\n', B_rms_ar/1e6);
fprintf('    中心频率: %.3f MHz\n', results_rms_ar.center_freq/1e6);

% 方法5：积分法（累积功率法）
fprintf('\n方法5：积分法（累积功率法）\n');
[B_integral_welch, results_integral_welch] = estimateBandwidthFromPSD(...
    Pxx_welch, f_welch, 'method', 'integral', 'plot', false);
fprintf('  Welch算法：估计带宽 = %.3f MHz\n', B_integral_welch/1e6);

[B_integral_ar, results_integral_ar] = estimateBandwidthFromPSD(...
    Pxx_ar, f_ar, 'method', 'integral', 'plot', false);
fprintf('  AR模型法：估计带宽 = %.3f MHz\n', B_integral_ar/1e6);

% 方法6：使用所有方法（自动对比）
fprintf('\n方法6：使用所有方法（自动对比）\n');
[B_all_welch, results_all_welch] = estimateBandwidthFromPSD(...
    Pxx_welch, f_welch, 'method', 'all', 'plot', true);
fprintf('  Welch算法 - 各方法估计结果：\n');
fprintf('    阈值法: %.3f MHz\n', B_all_welch.threshold/1e6);
fprintf('    能量法: %.3f MHz\n', B_all_welch.energy/1e6);
fprintf('    峰值下降法: %.3f MHz\n', B_all_welch.peak_drop/1e6);
fprintf('    RMS法: %.3f MHz\n', B_all_welch.rms/1e6);
fprintf('    积分法: %.3f MHz\n', B_all_welch.integral/1e6);
fprintf('    平均带宽: %.3f MHz (标准差: %.3f MHz)\n', ...
    B_all_welch.mean/1e6, B_all_welch.std/1e6);

[B_all_ar, results_all_ar] = estimateBandwidthFromPSD(...
    Pxx_ar, f_ar, 'method', 'all', 'plot', true);
fprintf('  AR模型法 - 各方法估计结果：\n');
fprintf('    阈值法: %.3f MHz\n', B_all_ar.threshold/1e6);
fprintf('    能量法: %.3f MHz\n', B_all_ar.energy/1e6);
fprintf('    峰值下降法: %.3f MHz\n', B_all_ar.peak_drop/1e6);
fprintf('    RMS法: %.3f MHz\n', B_all_ar.rms/1e6);
fprintf('    积分法: %.3f MHz\n', B_all_ar.integral/1e6);
fprintf('    平均带宽: %.3f MHz (标准差: %.3f MHz)\n', ...
    B_all_ar.mean/1e6, B_all_ar.std/1e6);

fprintf('\n  完成！\n\n');

%===============================================================================
% 步骤4: 结果汇总和性能分析
%===============================================================================
fprintf('========================================\n');
fprintf('带宽估计结果汇总\n');
fprintf('========================================\n');
fprintf('理论带宽: %.3f MHz (%.0f kHz)\n', B_ideal/1e6, B_ideal/1e3);
fprintf('\n');

fprintf('Welch算法估计结果：\n');
methods_welch = {'阈值法', '能量法', '峰值下降法', 'RMS法', '积分法', '平均'};
B_welch_values = [B_threshold_welch, B_energy_welch, B_peak_drop_welch, ...
    B_rms_welch, B_integral_welch, B_all_welch.mean];
for i = 1:length(methods_welch)
    error_abs = abs(B_welch_values(i) - B_ideal);
    error_rel = error_abs / B_ideal * 100;
    fprintf('  %s: %.3f MHz, 误差: %.3f MHz (%.2f%%)\n', ...
        methods_welch{i}, B_welch_values(i)/1e6, error_abs/1e6, error_rel);
end

fprintf('\nAR模型法估计结果：\n');
B_ar_values = [B_threshold_ar, B_energy_ar, B_peak_drop_ar, ...
    B_rms_ar, B_integral_ar, B_all_ar.mean];
for i = 1:length(methods_welch)
    error_abs = abs(B_ar_values(i) - B_ideal);
    error_rel = error_abs / B_ideal * 100;
    fprintf('  %s: %.3f MHz, 误差: %.3f MHz (%.2f%%)\n', ...
        methods_welch{i}, B_ar_values(i)/1e6, error_abs/1e6, error_rel);
end

fprintf('========================================\n\n');

%===============================================================================
% 步骤5: 绘制综合对比图
%===============================================================================
fprintf('正在绘制综合对比图...\n');

figure('Name', '带宽估计方法综合对比', 'Position', [50, 50, 1600, 900]);

% 子图1：Welch算法PSD和所有方法的带宽标记
subplot(2, 3, 1);
plot(f_welch/1e6, Pxx_welch, 'b-', 'LineWidth', 1.5);
hold on;
colors = lines(5);
methods_short = {'threshold', 'energy', 'peak_drop', 'rms', 'integral'};
for i = 1:length(methods_short)
    r = results_all_welch.(methods_short{i});
    plot([r.lower_freq/1e6, r.lower_freq/1e6], ylim, '--', ...
        'Color', colors(i,:), 'LineWidth', 1.5);
    plot([r.upper_freq/1e6, r.upper_freq/1e6], '--', ...
        'Color', colors(i,:), 'LineWidth', 1.5);
end
grid on;
xlabel('频率 (MHz)');
ylabel('PSD (dB)');
title('Welch算法 - 所有方法带宽估计');
legend(['PSD', methods_short], 'Location', 'best');
hold off;

% 子图2：AR模型PSD和所有方法的带宽标记
subplot(2, 3, 2);
plot(f_ar/1e6, Pxx_ar, 'r-', 'LineWidth', 1.5);
hold on;
for i = 1:length(methods_short)
    r = results_all_ar.(methods_short{i});
    plot([r.lower_freq/1e6, r.lower_freq/1e6], ylim, '--', ...
        'Color', colors(i,:), 'LineWidth', 1.5);
    plot([r.upper_freq/1e6, r.upper_freq/1e6], '--', ...
        'Color', colors(i,:), 'LineWidth', 1.5);
end
grid on;
xlabel('频率 (MHz)');
ylabel('PSD (dB)');
title('AR模型法 - 所有方法带宽估计');
legend(['PSD', methods_short], 'Location', 'best');
hold off;

% 子图3：Welch算法各方法估计值对比
subplot(2, 3, 3);
bar(B_welch_values(1:5)/1e6);
set(gca, 'XTickLabel', methods_welch(1:5));
ylabel('带宽 (MHz)');
title('Welch算法 - 各方法估计值');
grid on;
hold on;
plot([0.5, 5.5], [B_ideal/1e6, B_ideal/1e6], 'r--', 'LineWidth', 2, ...
    'DisplayName', '理论值');
plot([0.5, 5.5], [B_welch_values(6)/1e6, B_welch_values(6)/1e6], 'g--', ...
    'LineWidth', 2, 'DisplayName', '平均值');
legend('Location', 'best');
hold off;

% 子图4：AR模型各方法估计值对比
subplot(2, 3, 4);
bar(B_ar_values(1:5)/1e6);
set(gca, 'XTickLabel', methods_welch(1:5));
ylabel('带宽 (MHz)');
title('AR模型法 - 各方法估计值');
grid on;
hold on;
plot([0.5, 5.5], [B_ideal/1e6, B_ideal/1e6], 'r--', 'LineWidth', 2, ...
    'DisplayName', '理论值');
plot([0.5, 5.5], [B_ar_values(6)/1e6, B_ar_values(6)/1e6], 'g--', ...
    'LineWidth', 2, 'DisplayName', '平均值');
legend('Location', 'best');
hold off;

% 子图5：误差对比（Welch算法）
subplot(2, 3, 5);
errors_welch = abs(B_welch_values(1:5) - B_ideal) / B_ideal * 100;
bar(errors_welch);
set(gca, 'XTickLabel', methods_welch(1:5));
ylabel('相对误差 (%)');
title('Welch算法 - 各方法估计误差');
grid on;

% 子图6：误差对比（AR模型法）
subplot(2, 3, 6);
errors_ar = abs(B_ar_values(1:5) - B_ideal) / B_ideal * 100;
bar(errors_ar);
set(gca, 'XTickLabel', methods_welch(1:5));
ylabel('相对误差 (%)');
title('AR模型法 - 各方法估计误差');
grid on;

fprintf('图形绘制完成！\n\n');

fprintf('========================================\n');
fprintf('示例完成！\n');
fprintf('========================================\n');
fprintf('提示：\n');
fprintf('  1. 阈值法：简单快速，适合信噪比较好的情况\n');
fprintf('  2. 能量法：物理意义明确，适合能量集中的信号\n');
fprintf('  3. 峰值下降法：与阈值法类似，但相对于峰值\n');
fprintf('  4. RMS法：基于统计特性，适合高斯型频谱\n');
fprintf('  5. 积分法：基于累积功率，鲁棒性较好\n');
fprintf('  6. 建议：使用多种方法，取平均值或选择最稳定的方法\n');
fprintf('========================================\n');

