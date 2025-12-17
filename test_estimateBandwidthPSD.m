%===============================================================================
% test_estimateBandwidthPSD.m - 测试estimateBandwidthPSD.m函数
% 
% 功能说明:
%   1. 从.mat文件加载OFDM信号（发送和接收信号）
%     - 优先从指定的.mat文件加载信号
%     - 如果文件不存在，则从workspace加载或运行signalGenerate.m生成
%   2. 使用estimateBandwidthPSD.m对接收信号进行功率谱密度估计
%   3. 显示PSD估计结果
%
% 依赖关系:
%   - OFDM_signal_SNR15.0dB_300subcarriers_20251215_202859.mat：信号数据文件
%   - signalGenerate.m：生成OFDM信号（备用）
%   - estimateBandwidthPSD.m：功率谱密度估计
%
% 使用示例:
%   test_estimateBandwidthPSD
%
% 创建日期: 2025.12.10
% 修改日期: 2025.12.15 - 添加从.mat文件加载信号的功能
%===============================================================================

clc;
clear all;
close all;

fprintf('========================================\n');
fprintf('测试 estimateBandwidthPSD.m 函数\n');
fprintf('========================================\n\n');

%===============================================================================
% 步骤1: 加载OFDM信号（发送和接收）
%===============================================================================
% 功能：从.mat文件加载信号，如果文件不存在则使用workspace变量或生成新信号
fprintf('步骤1: 加载OFDM信号...\n');
fprintf('----------------------------------------\n');

% 指定要加载的.mat文件名
mat_filename = 'OFDM_signal_SNR15.0dB_300subcarriers_20251215_202859.mat';

try
    % 优先从.mat文件加载信号
    if exist(mat_filename, 'file')
        fprintf('从文件加载信号: %s\n', mat_filename);
        loaded_data = load(mat_filename);
        
        % 检查文件中是否包含必要的变量
        if isfield(loaded_data, 'Rx_data')
            Rx_data = loaded_data.Rx_data;
            fprintf('  成功加载 Rx_data\n');
        else
            error('文件中未找到 Rx_data 变量');
        end
        
        % 加载其他参数（如果存在）
        if isfield(loaded_data, 'Tx_data')
            Tx_data = loaded_data.Tx_data;
            fprintf('  成功加载 Tx_data\n');
        end
        if isfield(loaded_data, 'fs')
            fs = loaded_data.fs;
            fprintf('  成功加载 fs = %.2f MHz\n', fs/1e6);
        else
            fs = 15.36e6;  % 默认采样频率
            fprintf('  使用默认采样频率: %.2f MHz\n', fs/1e6);
        end
        if isfield(loaded_data, 'targetSNRdB')
            targetSNRdB = loaded_data.targetSNRdB;
            fprintf('  成功加载 targetSNRdB = %.1f dB\n', targetSNRdB);
        else
            targetSNRdB = 15;  % 默认SNR
            fprintf('  使用默认SNR: %.1f dB\n', targetSNRdB);
        end
        if isfield(loaded_data, 'carrier_count')
            carrier_count = loaded_data.carrier_count;
            fprintf('  成功加载 carrier_count = %d\n', carrier_count);
        else
            carrier_count = 300;  % 默认子载波数
            fprintf('  使用默认子载波数: %d\n', carrier_count);
        end
        if isfield(loaded_data, 'subcarrier_spacing')
            subcarrier_spacing = loaded_data.subcarrier_spacing;
            fprintf('  成功加载 subcarrier_spacing = %.1f kHz\n', subcarrier_spacing/1e3);
        else
            subcarrier_spacing = 15e3;  % 默认子载波间隔
            fprintf('  使用默认子载波间隔: %.1f kHz\n', subcarrier_spacing/1e3);
        end
        
    % 如果文件不存在，检查workspace中是否有信号变量
    elseif exist('Rx_data', 'var')
        fprintf('未找到.mat文件，使用workspace中已有的信号变量\n');
        fprintf('  使用workspace中的 Rx_data\n');
        
        % 获取或设置默认参数
        if ~exist('fs', 'var')
            fs = 15.36e6;
            fprintf('  使用默认采样频率: %.2f MHz\n', fs/1e6);
        end
        if ~exist('targetSNRdB', 'var')
            targetSNRdB = 15;
            fprintf('  使用默认SNR: %.1f dB\n', targetSNRdB);
        end
        if ~exist('carrier_count', 'var')
            carrier_count = 300;
            fprintf('  使用默认子载波数: %d\n', carrier_count);
        end
        if ~exist('subcarrier_spacing', 'var')
            subcarrier_spacing = 15e3;
            fprintf('  使用默认子载波间隔: %.1f kHz\n', subcarrier_spacing/1e3);
        end
        
    % 如果都没有，运行signalGenerate.m生成信号
    else
        fprintf('未找到.mat文件和workspace变量，运行 signalGenerate.m 生成信号...\n');
        signalGenerate;
        
        % 获取信号参数
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
    end
    
    % 计算理论带宽
    B_ideal = carrier_count * subcarrier_spacing;  % 理论带宽（Hz）
    
    fprintf('\n信号参数：\n');
    fprintf('  - 采样频率: %.2f MHz\n', fs/1e6);
    fprintf('  - 接收信号长度: %d 样本\n', length(Rx_data));
    fprintf('  - SNR: %.1f dB\n', targetSNRdB);
    fprintf('  - 子载波数: %d\n', carrier_count);
    fprintf('  - 子载波间隔: %.1f kHz\n', subcarrier_spacing/1e3);
    fprintf('  - 理论带宽: %.3f MHz (%.0f kHz)\n', B_ideal/1e6, B_ideal/1e3);
    fprintf('  完成！\n\n');
    
catch ME
    error('错误：无法加载或获取信号。\n错误信息：%s\n\n请确保：\n  1. .mat文件存在且包含Rx_data变量\n  2. 或者workspace中有Rx_data变量\n  3. 或者signalGenerate.m可以正常运行', ME.message);
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

% 清除原信号变量，只保留PSD估计结果
fprintf('清除原信号变量以节省内存...\n');
if exist('Rx_data', 'var')
    clear Rx_data;
    fprintf('  已清除 Rx_data\n');
end
if exist('Tx_data', 'var')
    clear Tx_data;
    fprintf('  已清除 Tx_data\n');
end
fprintf('  原信号已清除，仅保留PSD估计结果\n\n');

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
% 步骤4: 保存PSD估计结果到.mat文件
%===============================================================================
% 功能：将两种模型估计的功率谱密度及相关参数保存到.mat文件中
% 注意：不保存原信号（Rx_data, Tx_data），只保存PSD估计结果
% 保存内容：
%   - Pxx_welch, f_welch：Welch算法估计的PSD和频率向量
%   - Pxx_ar, f_ar：AR模型法估计的PSD和频率向量
%   - 相关参数：fs, targetSNRdB, carrier_count, subcarrier_spacing, B_ideal等
fprintf('步骤4: 保存PSD估计结果到.mat文件...\n');
fprintf('----------------------------------------\n');

% 生成文件名（包含SNR、子载波数、日期时间等信息）
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
psd_filename = sprintf('PSD_SNR%.1fdB_%dsubcarriers_%s.mat', ...
    targetSNRdB, carrier_count, timestamp);

try
    % 保存PSD估计结果和相关参数
    save(psd_filename, ...
        'Pxx_welch', 'f_welch', ...           % Welch算法PSD结果
        'Pxx_ar', 'f_ar', ...                  % AR模型法PSD结果
        'fs', 'targetSNRdB', ...               % 系统参数
        'carrier_count', 'subcarrier_spacing', ... % OFDM参数
        'B_ideal', ...                         % 理论带宽
        'peak_freq_welch', 'peak_value_welch', ... % Welch峰值信息
        'peak_freq_ar', 'peak_value_ar', ...   % AR峰值信息
        '-v7.3');                              % 使用v7.3格式以支持大文件
    
    fprintf('成功保存PSD估计结果到文件: %s\n', psd_filename);
    fprintf('保存的变量：\n');
    fprintf('  - Pxx_welch, f_welch: Welch算法PSD估计结果\n');
    fprintf('  - Pxx_ar, f_ar: AR模型法PSD估计结果\n');
    fprintf('  - fs, targetSNRdB: 系统参数\n');
    fprintf('  - carrier_count, subcarrier_spacing: OFDM参数\n');
    fprintf('  - B_ideal: 理论带宽\n');
    fprintf('  - peak_freq_welch, peak_value_welch: Welch峰值信息\n');
    fprintf('  - peak_freq_ar, peak_value_ar: AR峰值信息\n');
    fprintf('  完成！\n\n');
    
catch ME
    fprintf('警告：保存PSD结果失败。错误信息：%s\n', ME.message);
    fprintf('PSD数据仍在workspace中，可以手动保存\n\n');
end

%===============================================================================
% 步骤5: 显示PSD估计结果总结
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
% 步骤6: 绘制PSD估计结果对比图
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
fprintf('\n');
fprintf('变量说明：\n');
fprintf('  - B_ideal: 理论带宽 (Hz)\n');
fprintf('  - Pxx_welch, f_welch: Welch算法PSD估计结果\n');
fprintf('  - Pxx_ar, f_ar: AR模型法PSD估计结果\n');
fprintf('  - peak_freq_welch: Welch算法峰值频率 (Hz)\n');
fprintf('  - peak_freq_ar: AR模型法峰值频率 (Hz)\n');
fprintf('\n');
if exist('psd_filename', 'var')
    fprintf('PSD数据已保存到文件: %s\n', psd_filename);
    fprintf('可以使用以下命令加载：\n');
    fprintf('  load(''%s'')\n', psd_filename);
end
fprintf('========================================\n');

