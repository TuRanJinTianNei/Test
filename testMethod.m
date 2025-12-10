%===============================================================================
% testMethod.m - OFDM信号带宽估计测试脚本
% 
% 功能说明:
%   1. 生成OFDM信号（基带复数信号）
%   2. 使用method.m的两种算法（Welch和AR模型法）估计带宽
%   3. 对比估计结果与理论带宽
%
% 创建日期: 2025.12.10
%===============================================================================

clc;
clear all;
close all;

fprintf('========================================\n');
fprintf('OFDM信号带宽估计测试（method.m）\n');
fprintf('========================================\n\n');

%===============================================================================
% 【系统参数配置】
%===============================================================================
% 基本时间/频率参数
subcarrier_spacing = 15e3;         % 子载波间隔：15 kHz
fs = 15.36e6;                      % 采样频率：15.36 MHz

carrier_count = 300;               % 有效数据子载波数
total_symbols = 10;                % OFDM符号数（简化版本，使用较少符号）
bits_per_symbol = 4;               % 每个子载波承载的比特数（4=16QAM）

% IFFT/FFT参数
IFFT_bin_length = 1024;            % IFFT 点数

% 保护间隔参数
GI = 72;                           % 循环前缀长度（样本数）

% 信道参数
targetSNRdB = 15;                  % 目标信噪比（dB）

%===============================================================================
% 【OFDM信号生成】
%===============================================================================
fprintf('步骤1: 生成OFDM信号...\n');

% 步骤1: 随机比特流生成
baseband_out_length = carrier_count * total_symbols * bits_per_symbol;
rng('shuffle');  % 基于当前时间设置随机种子
baseband_out = randi([0 1], 1, baseband_out_length);

% 步骤2: 16QAM调制
% 注意：这里假设有qam16函数，如果没有可以使用MATLAB内置的qammod
try
    complex_carrier_matrix = qam16(baseband_out);
catch
    % 如果没有qam16函数，使用MATLAB内置的qammod
    complex_carrier_matrix = qammod(baseband_out, 16, 'bin', 'InputType', 'bit');
end
complex_carrier_matrix = reshape(complex_carrier_matrix', carrier_count, total_symbols)';

% 步骤3: 频域子载波映射（构造埃尔米特共轭对称，使IFFT输出为实信号）
carriers = (1:carrier_count) + (floor(IFFT_bin_length/4) - floor(carrier_count/2));
conjugate_carriers = IFFT_bin_length - carriers + 2;

IFFT_modulation = zeros(total_symbols, IFFT_bin_length);
IFFT_modulation(:, carriers) = complex_carrier_matrix;
IFFT_modulation(:, conjugate_carriers) = conj(complex_carrier_matrix);

% 步骤4: IFFT变换（频域 → 时域）
signal_after_IFFT = ifft(IFFT_modulation, IFFT_bin_length, 2);

% 步骤5: 添加循环前缀(CP)
XX = zeros(total_symbols, IFFT_bin_length + GI);
for k = 1:total_symbols
    % 符号主体部分
    XX(k, GI+1:IFFT_bin_length+GI) = signal_after_IFFT(k, :);
    % 循环前缀：将符号尾部复制到开头
    XX(k, 1:GI) = signal_after_IFFT(k, IFFT_bin_length-GI+1:IFFT_bin_length);
end

% 步骤6: 串并转换，生成发送信号
Tx_data = reshape(XX', 1, total_symbols*(IFFT_bin_length+GI));

fprintf('  - 子载波数: %d\n', carrier_count);
fprintf('  - OFDM符号数: %d\n', total_symbols);
fprintf('  - IFFT点数: %d\n', IFFT_bin_length);
fprintf('  - CP长度: %d\n', GI);
fprintf('  - 信号长度: %d 样本\n', length(Tx_data));
fprintf('  完成！\n\n');

%===============================================================================
% 【信道仿真（可选）】
%===============================================================================
fprintf('步骤2: 信道仿真...\n');
use_channel = true;  % 是否使用信道仿真

if use_channel
    % 添加AWGN噪声
    Tx_signal_power = var(Tx_data);
    linear_SNR = 10^(targetSNRdB/10);
    noise_sigma = Tx_signal_power / linear_SNR;
    noise_scale_factor = sqrt(noise_sigma);
    noise = randn(1, length(Tx_data)) * noise_scale_factor;
    Rx_data = Tx_data + noise;
    
    fprintf('  - 添加AWGN噪声（SNR = %.1f dB）\n', targetSNRdB);
else
    Rx_data = Tx_data;
    fprintf('  - 跳过信道仿真，使用原始信号\n');
end
fprintf('  完成！\n\n');

%===============================================================================
% 【带宽估计（使用method.m）】
%===============================================================================
fprintf('步骤3: 使用method.m进行带宽估计...\n');

% 计算理论带宽
B_ideal = carrier_count * subcarrier_spacing;  % 300 * 15e3 = 4.5 MHz
fprintf('\n理论带宽 B_ideal: %.6f Hz (%.3f MHz)\n', B_ideal, B_ideal/1e6);

% 使用method.m进行带宽估计
try
    % 使用实际SNR进行估计（如果使用信道，使用targetSNRdB；否则使用一个较大的值）
    if use_channel
        snr_for_method = targetSNRdB;
    else
        snr_for_method = 20;  % 无噪声时使用较大的SNR值
    end
    
    [B_welch, B_ar] = method(Rx_data, fs, snr_for_method);
    
    % 计算误差
    error_welch = abs(B_welch - B_ideal);
    error_welch_percent = (error_welch / B_ideal) * 100;
    
    error_ar = abs(B_ar - B_ideal);
    error_ar_percent = (error_ar / B_ideal) * 100;
    
    % 输出结果
    fprintf('\n【Welch算法估计结果】\n');
    fprintf('  估计带宽 B_welch: %.6f Hz (%.3f MHz)\n', B_welch, B_welch/1e6);
    fprintf('  绝对误差        : %.6f Hz (%.3f MHz)\n', error_welch, error_welch/1e6);
    fprintf('  相对误差        : %.2f%%\n', error_welch_percent);
    
    fprintf('\n【AR模型法估计结果】\n');
    fprintf('  估计带宽 B_ar   : %.6f Hz (%.3f MHz)\n', B_ar, B_ar/1e6);
    fprintf('  绝对误差        : %.6f Hz (%.3f MHz)\n', error_ar, error_ar/1e6);
    fprintf('  相对误差        : %.2f%%\n', error_ar_percent);
    
    fprintf('\n【方法对比】\n');
    if error_welch < error_ar
        fprintf('  Welch算法估计更准确（误差更小）\n');
        fprintf('  误差差异        : %.6f Hz (%.2f%%)\n', ...
            error_ar - error_welch, error_ar_percent - error_welch_percent);
    elseif error_ar < error_welch
        fprintf('  AR模型法估计更准确（误差更小）\n');
        fprintf('  误差差异        : %.6f Hz (%.2f%%)\n', ...
            error_welch - error_ar, error_welch_percent - error_ar_percent);
    else
        fprintf('  两种方法估计精度相同\n');
    end
    
catch ME
    fprintf('错误：带宽估计失败！\n');
    fprintf('错误信息：%s\n', ME.message);
    fprintf('可能原因：\n');
    fprintf('  1. method.m文件不在MATLAB路径中\n');
    fprintf('  2. 缺少必要的工具箱（Signal Processing Toolbox等）\n');
    fprintf('  3. 缺少qam16函数（如果使用自定义函数）\n');
    rethrow(ME);
end

fprintf('\n========================================\n');
fprintf('测试完成！\n');
fprintf('========================================\n');

%===============================================================================
% 【可视化（可选）】
%===============================================================================
plot_signal = true;  % 是否绘制信号

if plot_signal
    fprintf('\n步骤4: 绘制信号波形...\n');
    
    % 绘制时域信号
    figure('Name', 'OFDM信号时域波形');
    subplot(2, 1, 1);
    plot(1:min(1000, length(Tx_data)), Tx_data(1:min(1000, length(Tx_data))));
    grid on;
    xlabel('样本索引');
    ylabel('幅度');
    title('发送信号时域波形（前1000个样本）');
    
    subplot(2, 1, 2);
    plot(1:min(1000, length(Rx_data)), Rx_data(1:min(1000, length(Rx_data))));
    grid on;
    xlabel('样本索引');
    ylabel('幅度');
    if use_channel
        title(sprintf('接收信号时域波形（前1000个样本，SNR=%.1f dB）', targetSNRdB));
    else
        title('接收信号时域波形（前1000个样本，无噪声）');
    end
    
    % 绘制频谱
    figure('Name', 'OFDM信号频谱');
    Nfft = 2^nextpow2(length(Rx_data));
    Rx_Fz = fftshift(fft(Rx_data, Nfft));
    f_axis = (-Nfft/2:(Nfft/2-1)) * fs / Nfft;
    
    subplot(2, 1, 1);
    plot(f_axis/1e6, 20*log10(abs(Rx_Fz) / max(abs(Rx_Fz)) + eps));
    grid on;
    xlabel('频率 (MHz)');
    ylabel('幅度 (dB)');
    title('接收信号频谱（归一化）');
    xlim([-fs/2/1e6, fs/2/1e6]);
    
    subplot(2, 1, 2);
    plot(f_axis/1e6, 20*log10(abs(Rx_Fz) / max(abs(Rx_Fz)) + eps));
    grid on;
    xlabel('频率 (MHz)');
    ylabel('幅度 (dB)');
    title('接收信号频谱（局部放大，显示信号带宽）');
    xlim([-2, 2]);  % 放大到±2MHz范围
    hold on;
    % 标记理论带宽
    plot([-B_ideal/2/1e6, B_ideal/2/1e6], [-3, -3], 'r--', 'LineWidth', 2);
    text(0, -5, sprintf('理论带宽: %.1f MHz', B_ideal/1e6), ...
        'HorizontalAlignment', 'center', 'Color', 'r', 'FontSize', 12);
    hold off;
    
    fprintf('  完成！\n');
end

fprintf('\n所有步骤完成！\n');
