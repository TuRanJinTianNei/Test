%===============================================================================
% test_bandwidth_methods_comparison.m - 对比多种带宽估计方法
% 
% 功能说明:
%   对比多种带宽估计方法的效果，帮助选择最适合的方法
%
% 依赖关系:
%   - signalGenerate.m：生成OFDM信号
%   - estimateBandwidthPSD.m：功率谱密度估计
%   - estimateBandwidthMethods.m：多种带宽估计方法
%
% 使用示例:
%   test_bandwidth_methods_comparison
%
% 创建日期: 2025.12.10
%===============================================================================

clc;
clear all;
close all;

fprintf('========================================\n');
fprintf('带宽估计方法对比测试\n');
fprintf('========================================\n\n');

%===============================================================================
% 步骤1: 生成OFDM信号
%===============================================================================
fprintf('步骤1: 生成OFDM信号...\n');
fprintf('----------------------------------------\n');

try
    if ~exist('Rx_data', 'var') || ~exist('Tx_data', 'var') || ~exist('fs', 'var')
        fprintf('未找到信号变量，运行 signalGenerate.m...\n');
        signalGenerate;
    else
        fprintf('使用workspace中已有的信号变量\n');
    end
    
    if ~exist('fs', 'var')
        fs = 15.36e6;
    end
    if ~exist('targetSNRdB', 'var')
        targetSNRdB = 15;
    end
    if ~exist('carrier_count', 'var')
        carrier_count = 300;
    end
    if ~exist('subcarrier_spacing', 'var')
        subcarrier_spacing = 15e3;
    end
    
    B_ideal = carrier_count * subcarrier_spacing;
    
    fprintf('信号参数：\n');
    fprintf('  - 采样频率: %.2f MHz\n', fs/1e6);
    fprintf('  - SNR: %.1f dB\n', targetSNRdB);
    fprintf('  - 理论带宽: %.3f MHz (%.0f kHz)\n', B_ideal/1e6, B_ideal/1e3);
    fprintf('  完成！\n\n');
    
catch ME
    error('错误：无法生成或获取信号。\n错误信息：%s', ME.message);
end

%===============================================================================
% 步骤2: 功率谱密度估计
%===============================================================================
fprintf('步骤2: 进行功率谱密度估计...\n');
fprintf('----------------------------------------\n');

[Pxx_welch, f_welch, Pxx_ar, f_ar] = estimateBandwidthPSD(...
    Rx_data, fs, targetSNRdB, 'plot', false);

fprintf('功率谱密度估计完成！\n\n');

%===============================================================================
% 步骤3: 使用多种方法估计带宽
%===============================================================================
fprintf('步骤3: 使用多种方法估计带宽...\n');
fprintf('----------------------------------------\n');

methods = {'threshold', 'energy_percent', 'peak_drop', ...
           'rms', 'adaptive', 'gradient'};
method_names = {'阈值法(-3dB)', '能量百分比法(90%)', ...
                '峰值下降法(-3dB)', 'RMS带宽法', '自适应阈值法', '梯度边缘检测法'};

results_welch = struct();
results_ar = struct();

for i = 1:length(methods)
    method = methods{i};
    fprintf('\n方法 %d/%d: %s\n', i, length(methods), method_names{i});
    
    try
        % Welch算法结果
        switch method
            case 'threshold'
                [bl_w, bu_w, B_w, d_w] = estimateBandwidthMethods(...
                    Pxx_welch, f_welch, method, 'threshold_db', -3);
            case 'energy_percent'
                [bl_w, bu_w, B_w, d_w] = estimateBandwidthMethods(...
                    Pxx_welch, f_welch, method, 'energy_percent', 90);
            case 'peak_drop'
                [bl_w, bu_w, B_w, d_w] = estimateBandwidthMethods(...
                    Pxx_welch, f_welch, method, 'peak_drop_db', -3);
            case 'adaptive'
                [bl_w, bu_w, B_w, d_w] = estimateBandwidthMethods(...
                    Pxx_welch, f_welch, method, 'adaptive_factor', 0.1);
            otherwise
                [bl_w, bu_w, B_w, d_w] = estimateBandwidthMethods(...
                    Pxx_welch, f_welch, method);
        end
        
        % AR模型结果
        switch method
            case 'threshold'
                [bl_a, bu_a, B_a, d_a] = estimateBandwidthMethods(...
                    Pxx_ar, f_ar, method, 'threshold_db', -3);
            case 'energy_percent'
                [bl_a, bu_a, B_a, d_a] = estimateBandwidthMethods(...
                    Pxx_ar, f_ar, method, 'energy_percent', 90);
            case 'peak_drop'
                [bl_a, bu_a, B_a, d_a] = estimateBandwidthMethods(...
                    Pxx_ar, f_ar, method, 'peak_drop_db', -3);
            case 'adaptive'
                [bl_a, bu_a, B_a, d_a] = estimateBandwidthMethods(...
                    Pxx_ar, f_ar, method, 'adaptive_factor', 0.1);
            otherwise
                [bl_a, bu_a, B_a, d_a] = estimateBandwidthMethods(...
                    Pxx_ar, f_ar, method);
        end
        
        results_welch.(method).band_lower = bl_w;
        results_welch.(method).band_upper = bu_w;
        results_welch.(method).bandwidth = B_w;
        results_welch.(method).error = abs(B_w - B_ideal) / B_ideal * 100;
        results_welch.(method).details = d_w;
        
        results_ar.(method).band_lower = bl_a;
        results_ar.(method).band_upper = bu_a;
        results_ar.(method).bandwidth = B_a;
        results_ar.(method).error = abs(B_a - B_ideal) / B_ideal * 100;
        results_ar.(method).details = d_a;
        
        fprintf('  Welch: %.3f MHz (误差: %.2f%%)\n', B_w/1e6, results_welch.(method).error);
        fprintf('  AR:    %.3f MHz (误差: %.2f%%)\n', B_a/1e6, results_ar.(method).error);
    catch ME
        
        fprintf('  错误: %s\n', ME.message);
        results_welch.(method).bandwidth = NaN;
        results_ar.(method).bandwidth = NaN;
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
legend([{'PSD'}, method_names(~isnan([results_welch.(methods{1}).bandwidth]))], ...
       'Location', 'best', 'FontSize', 8);
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
legend([{'PSD'}, method_names(~isnan([results_ar.(methods{1}).bandwidth]))], ...
       'Location', 'best', 'FontSize', 8);
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

fprintf('========================================\n');
fprintf('测试完成！\n');
fprintf('========================================\n');
fprintf('建议：根据误差对比结果选择最适合您信号的方法\n');
fprintf('========================================\n');
