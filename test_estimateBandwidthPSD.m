%===============================================================================
% test_estimateBandwidthPSD.m - 测试estimateBandwidthPSD.m函数
% 
% 功能说明:
%   1. 运行signalGenerate.m生成OFDM信号（发送和接收信号）
%   2. 使用estimateBandwidthPSD.m对接收信号进行功率谱密度估计
%   3. 显示PSD估计结果
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
fprintf('  注意：对于实信号，只估计正频率部分（0到fs/2），频谱共轭对称\n');
fprintf('  单边谱转换规则：正频率部分（除DC和Nyquist）幅度×2\n');
fprintf('  完成！\n\n');

%===============================================================================
% 步骤3: 显示PSD估计结果
%===============================================================================
fprintf('步骤3: PSD估计结果分析...\n');
fprintf('----------------------------------------\n');

% 找到峰值位置
[~, peak_idx_welch] = max(Pxx_welch);
peak_freq_welch = f_welch(peak_idx_welch);
peak_value_welch = Pxx_welch(peak_idx_welch);

[~, peak_idx_ar] = max(Pxx_ar);
peak_freq_ar = f_ar(peak_idx_ar);
peak_value_ar = Pxx_ar(peak_idx_ar);

fprintf('Welch算法PSD特征：\n');
fprintf('  - 峰值频率: %.3f MHz\n', peak_freq_welch/1e6);
fprintf('  - 峰值功率: %.2f dB\n', peak_value_welch);
fprintf('  - 频率范围: %.3f - %.3f MHz (正频率部分，0到fs/2)\n', f_welch(1)/1e6, f_welch(end)/1e6);
fprintf('  - 奈奎斯特频率: %.3f MHz\n', fs/2/1e6);

fprintf('\nAR模型法PSD特征：\n');
fprintf('  - 峰值频率: %.3f MHz\n', peak_freq_ar/1e6);
fprintf('  - 峰值功率: %.2f dB\n', peak_value_ar);
fprintf('  - 频率范围: %.3f - %.3f MHz (正频率部分，0到fs/2)\n', f_ar(1)/1e6, f_ar(end)/1e6);
fprintf('  - 奈奎斯特频率: %.3f MHz\n', fs/2/1e6);

fprintf('  完成！\n\n');

%===============================================================================
% 步骤4: 显示PSD估计结果总结
%===============================================================================
fprintf('========================================\n');
fprintf('PSD估计结果总结\n');
fprintf('========================================\n');
fprintf('理论带宽: %.3f MHz (%.0f kHz)\n', B_ideal/1e6, B_ideal/1e3);
fprintf('\n');
fprintf('Welch算法PSD估计：\n');
fprintf('  - 峰值频率: %.3f MHz\n', peak_freq_welch/1e6);
fprintf('  - 峰值功率: %.2f dB\n', peak_value_welch);
fprintf('  - 频率分辨率: %.3f kHz\n', (f_welch(2) - f_welch(1))/1e3);
fprintf('\n');
fprintf('AR模型法PSD估计：\n');
fprintf('  - 峰值频率: %.3f MHz\n', peak_freq_ar/1e6);
fprintf('  - 峰值功率: %.2f dB\n', peak_value_ar);
fprintf('  - 频率分辨率: %.3f kHz\n', (f_ar(2) - f_ar(1))/1e3);
fprintf('========================================\n\n');

%===============================================================================
% 步骤5: 绘制PSD估计结果对比图
%===============================================================================
fprintf('正在绘制PSD估计结果对比图...\n');

figure('Name', 'PSD估计结果对比', 'Position', [100, 100, 1400, 600]);

% 子图1: Welch算法PSD
subplot(1, 3, 1);
plot(f_welch/1e6, Pxx_welch, 'b-', 'LineWidth', 1.5);
hold on;
plot([peak_freq_welch/1e6, peak_freq_welch/1e6], ylim, 'g--', 'LineWidth', 1.5, 'DisplayName', '峰值');
grid on;
xlabel('频率 (MHz)');
ylabel('PSD (dB)');
title(sprintf('Welch算法PSD估计\n峰值频率: %.3f MHz', peak_freq_welch/1e6));
legend('Location', 'best');
hold off;

% 子图2: AR模型PSD
subplot(1, 3, 2);
plot(f_ar/1e6, Pxx_ar, 'r-', 'LineWidth', 1.5);
hold on;
plot([peak_freq_ar/1e6, peak_freq_ar/1e6], ylim, 'g--', 'LineWidth', 1.5, 'DisplayName', '峰值');
grid on;
xlabel('频率 (MHz)');
ylabel('PSD (dB)');
title(sprintf('AR模型法PSD估计\n峰值频率: %.3f MHz', peak_freq_ar/1e6));
legend('Location', 'best');
hold off;

% 子图3: 两种方法PSD对比
subplot(1, 3, 3);
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

fprintf('图形绘制完成！\n\n');

fprintf('========================================\n');
fprintf('测试完成！\n');
fprintf('========================================\n');
fprintf('所有结果已显示在上方的图形窗口中\n');
fprintf('变量说明：\n');
fprintf('  - B_ideal: 理论带宽 (Hz)\n');
fprintf('  - Pxx_welch, f_welch: Welch算法PSD估计结果\n');
fprintf('  - Pxx_ar, f_ar: AR模型法PSD估计结果\n');
fprintf('  - peak_freq_welch: Welch算法峰值频率 (Hz)\n');
fprintf('  - peak_freq_ar: AR模型法峰值频率 (Hz)\n');
fprintf('========================================\n');

