%===============================================================================
% test_bandwidth_methods_comparison.m - 对比多种带宽估计方法
% 
% 功能说明:
%   对比多种带宽估计方法的效果，帮助选择最适合的方法。
%   优先从.mat文件加载已保存的功率谱密度数据，如果文件不存在则重新计算。
%
% 依赖关系:
%   - PSD_SNR15.0dB_300subcarriers_20251216_111043.mat：保存的PSD数据文件（优先使用）
%   - signalGenerate.m：生成OFDM信号（备用）
%   - estimateBandwidthPSD.m：功率谱密度估计（备用）
%   - estimateBandwidthFromPSD.m：多种带宽估计方法实现
%
% 使用示例:
%   test_bandwidth_methods_comparison
%
% 创建日期: 2025.12.10
% 修改日期: 2025.12.15 - 添加从.mat文件加载PSD数据的功能
%===============================================================================

clc;
clear all;
close all;

fprintf('========================================\n');
fprintf('带宽估计方法对比测试\n');
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
% 指定要加载的PSD数据文件名
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

% 定义要测试的方法列表（对应estimateBandwidthFromPSD.m支持的方法）
methods = {'threshold', 'energy', 'peak_drop', 'rms', 'integral'};
method_names = {'阈值法(-3dB)', '能量百分比法(90%)', ...
                '峰值下降法(-3dB)', 'RMS带宽法', '积分法'};

results_welch = struct();
results_ar = struct();

for i = 1:length(methods)
    method = methods{i};
    fprintf('\n方法 %d/%d: %s\n', i, length(methods), method_names{i});
    
    try
        % 对Welch算法估计的PSD进行带宽估计
        switch method
            case 'threshold'
                [B_w, r_w] = estimateBandwidthFromPSD(...
                    Pxx_welch, f_welch, 'method', 'threshold', ...
                    'threshold', -3, 'plot', false);
            case 'energy'
                [B_w, r_w] = estimateBandwidthFromPSD(...
                    Pxx_welch, f_welch, 'method', 'energy', ...
                    'energy_percent', 0.9, 'plot', false);
            case 'peak_drop'
                [B_w, r_w] = estimateBandwidthFromPSD(...
                    Pxx_welch, f_welch, 'method', 'peak_drop', ...
                    'peak_drop_db', 3, 'plot', false);
            case 'rms'
                [B_w, r_w] = estimateBandwidthFromPSD(...
                    Pxx_welch, f_welch, 'method', 'rms', 'plot', false);
            case 'integral'
                [B_w, r_w] = estimateBandwidthFromPSD(...
                    Pxx_welch, f_welch, 'method', 'integral', 'plot', false);
        end
        
        % 对AR模型法估计的PSD进行带宽估计
        switch method
            case 'threshold'
                [B_a, r_a] = estimateBandwidthFromPSD(...
                    Pxx_ar, f_ar, 'method', 'threshold', ...
                    'threshold', -3, 'plot', false);
            case 'energy'
                [B_a, r_a] = estimateBandwidthFromPSD(...
                    Pxx_ar, f_ar, 'method', 'energy', ...
                    'energy_percent', 0.9, 'plot', false);
            case 'peak_drop'
                [B_a, r_a] = estimateBandwidthFromPSD(...
                    Pxx_ar, f_ar, 'method', 'peak_drop', ...
                    'peak_drop_db', 3, 'plot', false);
            case 'rms'
                [B_a, r_a] = estimateBandwidthFromPSD(...
                    Pxx_ar, f_ar, 'method', 'rms', 'plot', false);
            case 'integral'
                [B_a, r_a] = estimateBandwidthFromPSD(...
                    Pxx_ar, f_ar, 'method', 'integral', 'plot', false);
        end
        
        % 保存Welch算法结果
        results_welch.(method).band_lower = r_w.lower_freq;
        results_welch.(method).band_upper = r_w.upper_freq;
        results_welch.(method).bandwidth = B_w;
        results_welch.(method).error = abs(B_w - B_ideal) / B_ideal * 100;
        results_welch.(method).peak_freq = r_w.peak_freq;
        results_welch.(method).peak_power = r_w.peak_power;
        
        % 保存AR模型法结果
        results_ar.(method).band_lower = r_a.lower_freq;
        results_ar.(method).band_upper = r_a.upper_freq;
        results_ar.(method).bandwidth = B_a;
        results_ar.(method).error = abs(B_a - B_ideal) / B_ideal * 100;
        results_ar.(method).peak_freq = r_a.peak_freq;
        results_ar.(method).peak_power = r_a.peak_power;
        
        % 显示估计结果
        fprintf('  Welch: %.3f MHz (误差: %.2f%%)\n', B_w/1e6, results_welch.(method).error);
        fprintf('  AR:    %.3f MHz (误差: %.2f%%)\n', B_a/1e6, results_ar.(method).error);
        
    catch ME
        % 如果方法执行失败，记录错误信息
        fprintf('  错误: %s\n', ME.message);
        results_welch.(method).bandwidth = NaN;
        results_welch.(method).error = NaN;
        results_ar.(method).bandwidth = NaN;
        results_ar.(method).error = NaN;
    end
end

fprintf('\n  完成！\n\n');

%===============================================================================
% 步骤4: 显示对比结果
%===============================================================================
fprintf('========================================\n');
fprintf('带宽估计方法对比结果\n');
fprintf('========================================\n');
fprintf('理论带宽: %.3f MHz (%.0f kHz)\n\n', B_ideal/1e6, B_ideal/1e3);

fprintf('%-25s %12s %12s %12s %12s\n', '方法', 'Welch带宽(MHz)', 'Welch误差(%)', ...
        'AR带宽(MHz)', 'AR误差(%)');
fprintf('%s\n', repmat('-', 1, 75));

for i = 1:length(methods)
    method = methods{i};
    if ~isnan(results_welch.(method).bandwidth)
        fprintf('%-25s %12.3f %12.2f %12.3f %12.2f\n', ...
            method_names{i}, ...
            results_welch.(method).bandwidth/1e6, ...
            results_welch.(method).error, ...
            results_ar.(method).bandwidth/1e6, ...
            results_ar.(method).error);
    else
        fprintf('%-25s %12s %12s %12s %12s\n', method_names{i}, 'N/A', 'N/A', 'N/A', 'N/A');
    end
end

fprintf('========================================\n\n');

%===============================================================================
% 步骤5: 绘制对比图
%===============================================================================
fprintf('正在绘制对比图...\n');

figure('Name', '带宽估计方法对比', 'Position', [100, 100, 1600, 1000]);

% 子图1: Welch算法 - 所有方法对比
subplot(2, 3, 1);
plot(f_welch/1e6, Pxx_welch, 'k-', 'LineWidth', 1);
hold on;
colors = lines(length(methods));
for i = 1:length(methods)
    method = methods{i};
    if ~isnan(results_welch.(method).bandwidth)
        bl = results_welch.(method).band_lower / 1e6;
        bu = results_welch.(method).band_upper / 1e6;
        plot([bl, bl], ylim, '--', 'Color', colors(i,:), 'LineWidth', 1.5);
        plot([bu, bu], '--', 'Color', colors(i,:), 'LineWidth', 1.5);
    end
end
grid on;
xlabel('频率 (MHz)');
ylabel('PSD (dB)');
title('Welch算法 - 各方法带宽边界对比');
% 构建图例：PSD + 有效的方法名称
legend_entries = {'PSD'};
for i = 1:length(methods)
    if ~isnan(results_welch.(methods{i}).bandwidth)
        legend_entries{end+1} = method_names{i};
    end
end
legend(legend_entries, 'Location', 'best', 'FontSize', 8);
hold off;

% 子图2: AR模型 - 所有方法对比
subplot(2, 3, 2);
plot(f_ar/1e6, Pxx_ar, 'k-', 'LineWidth', 1);
hold on;
for i = 1:length(methods)
    method = methods{i};
    if ~isnan(results_ar.(method).bandwidth)
        bl = results_ar.(method).band_lower / 1e6;
        bu = results_ar.(method).band_upper / 1e6;
        plot([bl, bl], ylim, '--', 'Color', colors(i,:), 'LineWidth', 1.5);
        plot([bu, bu], '--', 'Color', colors(i,:), 'LineWidth', 1.5);
    end
end
grid on;
xlabel('频率 (MHz)');
ylabel('PSD (dB)');
title('AR模型法 - 各方法带宽边界对比');
% 构建图例：PSD + 有效的方法名称
legend_entries_ar = {'PSD'};
for i = 1:length(methods)
    if ~isnan(results_ar.(methods{i}).bandwidth)
        legend_entries_ar{end+1} = method_names{i};
    end
end
legend(legend_entries_ar, 'Location', 'best', 'FontSize', 8);
hold off;

% 子图3: 带宽估计值对比（Welch）
subplot(2, 3, 3);
bw_welch = [];
method_names_valid = {};
for i = 1:length(methods)
    method = methods{i};
    if ~isnan(results_welch.(method).bandwidth)
        bw_welch(end+1) = results_welch.(method).bandwidth / 1e6;
        method_names_valid{end+1} = method_names{i};
    end
end
if ~isempty(bw_welch)
    bar(bw_welch);
    hold on;
    plot([1, length(bw_welch)], [B_ideal/1e6, B_ideal/1e6], 'r--', 'LineWidth', 2);
    set(gca, 'XTickLabel', method_names_valid, 'XTickLabelRotation', 45);
    ylabel('带宽 (MHz)');
    title('Welch算法 - 带宽估计值对比');
    legend('估计值', '理论值', 'Location', 'best');
    grid on;
    hold off;
end

% 子图4: 带宽估计值对比（AR）
subplot(2, 3, 4);
bw_ar = [];
method_names_valid_ar = {};
for i = 1:length(methods)
    method = methods{i};
    if ~isnan(results_ar.(method).bandwidth)
        bw_ar(end+1) = results_ar.(method).bandwidth / 1e6;
        method_names_valid_ar{end+1} = method_names{i};
    end
end
if ~isempty(bw_ar)
    bar(bw_ar);
    hold on;
    plot([1, length(bw_ar)], [B_ideal/1e6, B_ideal/1e6], 'r--', 'LineWidth', 2);
    set(gca, 'XTickLabel', method_names_valid_ar, 'XTickLabelRotation', 45);
    ylabel('带宽 (MHz)');
    title('AR模型法 - 带宽估计值对比');
    legend('估计值', '理论值', 'Location', 'best');
    grid on;
    hold off;
end

% 子图5: 相对误差对比（Welch）
subplot(2, 3, 5);
err_welch = [];
for i = 1:length(methods)
    method = methods{i};
    if ~isnan(results_welch.(method).bandwidth)
        err_welch(end+1) = results_welch.(method).error;
    end
end
if ~isempty(err_welch)
    bar(err_welch);
    set(gca, 'XTickLabel', method_names_valid, 'XTickLabelRotation', 45);
    ylabel('相对误差 (%)');
    title('Welch算法 - 相对误差对比');
    grid on;
end

% 子图6: 相对误差对比（AR）
subplot(2, 3, 6);
err_ar = [];
for i = 1:length(methods)
    method = methods{i};
    if ~isnan(results_ar.(method).bandwidth)
        err_ar(end+1) = results_ar.(method).error;
    end
end
if ~isempty(err_ar)
    bar(err_ar);
    set(gca, 'XTickLabel', method_names_valid_ar, 'XTickLabelRotation', 45);
    ylabel('相对误差 (%)');
    title('AR模型法 - 相对误差对比');
    grid on;
end

fprintf('图形绘制完成！\n\n');

%===============================================================================
% 步骤6: 单独绘制功率谱密度
%===============================================================================
fprintf('正在绘制功率谱密度图...\n');

figure('Name', '功率谱密度估计结果', 'Position', [200, 200, 1400, 600]);

% 子图1: Welch算法估计的PSD
subplot(1, 2, 1);
plot(f_welch/1e6, Pxx_welch, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('频率 (MHz)', 'FontSize', 12);
ylabel('功率谱密度 (dB)', 'FontSize', 12);
title(sprintf('Welch算法估计的功率谱密度\n(SNR=%.1f dB, 子载波数=%d, 理论带宽=%.3f MHz)', ...
    targetSNRdB, carrier_count, B_ideal/1e6), 'FontSize', 13);
legend('Welch PSD', 'Location', 'best', 'FontSize', 10);

% 子图2: AR模型法估计的PSD
subplot(1, 2, 2);
plot(f_ar/1e6, Pxx_ar, 'g-', 'LineWidth', 1.5);
grid on;
xlabel('频率 (MHz)', 'FontSize', 12);
ylabel('功率谱密度 (dB)', 'FontSize', 12);
title(sprintf('AR模型法估计的功率谱密度\n(SNR=%.1f dB, 子载波数=%d, 理论带宽=%.3f MHz)', ...
    targetSNRdB, carrier_count, B_ideal/1e6), 'FontSize', 13);
legend('AR PSD', 'Location', 'best', 'FontSize', 10);

fprintf('功率谱密度图绘制完成！\n\n');

fprintf('========================================\n');
fprintf('测试完成！\n');
fprintf('========================================\n');
fprintf('结果说明：\n');
fprintf('  1. 命令窗口显示了各方法的估计结果和误差对比表\n');
fprintf('  2. 图形窗口1显示了6个子图的详细可视化对比\n');
fprintf('  3. 图形窗口2单独显示了功率谱密度估计结果\n');
fprintf('  4. 建议根据误差对比结果选择最适合您信号的方法\n');
fprintf('\n');
fprintf('方法选择建议：\n');
fprintf('  - 高SNR（>10dB）：推荐使用阈值法或峰值下降法\n');
fprintf('  - 低SNR（<5dB）：推荐使用能量法或积分法\n');
fprintf('  - 需要统计特性：推荐使用RMS带宽法\n');
fprintf('  - 最佳实践：使用多种方法，取平均值或选择最稳定的方法\n');
fprintf('========================================\n');
