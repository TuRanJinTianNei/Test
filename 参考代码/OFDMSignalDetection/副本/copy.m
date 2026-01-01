%===============================================================================
% signalGenerate.m - OFDM信号生成函数（可独立运行）
% 
% 功能说明:
%   1. 生成完整的OFDM发送信号（基带复数信号）
%   2. 支持16QAM调制、IFFT、循环前缀、升余弦加窗等完整流程
%   3. 采用LTE风格的符号重叠相加机制
%   4. 可独立运行，不依赖其他脚本
%
% 创建日期: 2025.12.10
% 基于: test3.m的信号生成部分
%===============================================================================

clc;
clear all;
close all;

fprintf('========================================\n');
fprintf('OFDM信号生成程序（signalGenerate.m）\n');
fprintf('========================================\n\n');

%===============================================================================
% 【系统参数配置】
%===============================================================================
% 基本时间/频率参数
subcarrier_spacing = 15e3;         % 子载波间隔：15 kHz
fs = 15.36e6;                      % 采样频率：15.36 MHz

carrier_count = 300;               % 有效数据子载波数（正频率数据，负频率为镜像）
symbols_per_frame = 50;            % 每帧 OFDM 符号数
total_symbols = 100;                % 总共要传输的 OFDM 符号数
symbols_per_carrier = total_symbols; % 为兼容后续矩阵尺寸，沿用原变量表示"总符号数"
bits_per_symbol = 4;               % 每个子载波承载的比特数（4=16QAM）

% IFFT/FFT参数
IFFT_bin_length = 1024;            % IFFT 点数（包含有效、空、镜像与DC）

% 保护间隔参数
PrefixRatio = 72/ IFFT_bin_length; % 循环前缀比例（使 GI 固定为 72）
GI = 72;                            % 保护前缀长度（样本数）

% 加窗参数（LTE风格）
beta = 1/16;                        % 加窗滚降系数（过渡带占比，越小过渡越短）
GIP = beta*(IFFT_bin_length+GI);    % 右端后缀长度：配合窗的尾部过渡
GIP = min(GIP, GI);                 % LTE要求：后缀长度不超过CP长度，确保重叠区域在CP内
GIP = floor(GIP);                   % 确保为整数

% 导频参数（可选，默认不使用）
use_pilot_equalization = false;     % 是否使用导频进行信道估计和均衡
pilot_spacing = 6;                  % 导频间隔（如果使用导频）
pilot_symbol = 1 + 1i*0;            % 导频符号（已知的复数值，归一化功率）

% 输出参数
fprintf('========== 系统参数 ==========\n');
fprintf('子载波间隔: %.1f kHz\n', subcarrier_spacing/1e3);
fprintf('采样频率: %.2f MHz\n', fs/1e6);
fprintf('有效数据子载波数: %d\n', carrier_count);
fprintf('OFDM符号总数: %d\n', total_symbols);
fprintf('每帧符号数: %d\n', symbols_per_frame);
fprintf('IFFT点数: %d\n', IFFT_bin_length);
fprintf('循环前缀长度: %d\n', GI);
fprintf('后缀长度: %d\n', GIP);
fprintf('加窗滚降系数: 1/%d\n', round(1/beta));
fprintf('================================\n\n');

%===============================================================================
% 【发送端处理流程】
%===============================================================================

%------------------------------------------------------------------------------
% 步骤1: 随机比特流生成
%------------------------------------------------------------------------------
% 计算实际数据子载波数量（如果使用导频，需要排除导频占用的子载波）
if use_pilot_equalization
    % 计算导频位置
    carriers_temp = (1:carrier_count) + (floor(IFFT_bin_length/4) - floor(carrier_count/2));
    pilot_positions_temp = 1:pilot_spacing:carrier_count;
    data_carrier_count_actual = carrier_count - length(pilot_positions_temp);
else
    data_carrier_count_actual = carrier_count;
end

baseband_out_length = data_carrier_count_actual * symbols_per_carrier * bits_per_symbol;
% 使用基于时间的随机种子，确保每次运行生成不同的随机数据
rng('shuffle');  % 基于当前时间设置随机种子，确保每次运行都不同
baseband_out = randi([0 1], 1, baseband_out_length);

fprintf('步骤1: 生成随机比特流...\n');
fprintf('  - 比特流长度: %d bits\n', baseband_out_length);
fprintf('  - 数据子载波数: %d\n', data_carrier_count_actual);
fprintf('  完成！\n\n');

%------------------------------------------------------------------------------
% 步骤2: 16QAM调制（仅对数据子载波）
%------------------------------------------------------------------------------
fprintf('步骤2: 16QAM调制...\n');
try
    complex_carrier_matrix = qam16(baseband_out);
catch
    % 如果没有qam16函数，使用MATLAB内置的qammod
    fprintf('  警告：未找到qam16函数，使用MATLAB内置qammod\n');
    complex_carrier_matrix = qammod(baseband_out, 16, 'bin', 'InputType', 'bit');
end
complex_carrier_matrix = reshape(complex_carrier_matrix', data_carrier_count_actual, symbols_per_carrier)';
fprintf('  - 调制符号数: %d\n', length(complex_carrier_matrix(:)));
fprintf('  完成！\n\n');

%------------------------------------------------------------------------------
% 步骤2.5: 导频插入准备（计算导频位置和数据子载波位置）
%------------------------------------------------------------------------------
% 子载波索引计算：进行子载波共轭映射，使得OFDM符号经过IFFT之后是实信号
carriers = (1:carrier_count) + (floor(IFFT_bin_length/4) - floor(carrier_count/2));
conjugate_carriers = IFFT_bin_length - carriers + 2;

% 计算导频位置和数据子载波位置
if use_pilot_equalization
    % 在carriers范围内，每隔pilot_spacing个子载波插入一个导频
    pilot_positions_in_carriers = 1:pilot_spacing:carrier_count;  % 导频在carriers中的相对位置
    pilot_carriers = carriers(pilot_positions_in_carriers);  % 导频的绝对bin位置
    pilot_conjugate_carriers = conjugate_carriers(pilot_positions_in_carriers);  % 导频的共轭位置
    
    % 数据子载波位置（排除导频位置）
    data_positions_in_carriers = setdiff(1:carrier_count, pilot_positions_in_carriers);
    data_carriers = carriers(data_positions_in_carriers);  % 数据子载波的绝对bin位置
    data_conjugate_carriers = conjugate_carriers(data_positions_in_carriers);  % 数据子载波的共轭位置
    
    % 数据矩阵直接使用complex_carrier_matrix（已经是data_carrier_count_actual大小）
    data_matrix = complex_carrier_matrix;
    pilot_matrix = repmat(pilot_symbol, symbols_per_carrier, length(pilot_positions_in_carriers));  % 导频矩阵
else
    % 不使用导频，所有子载波都是数据子载波
    data_carriers = carriers;
    data_conjugate_carriers = conjugate_carriers;
    pilot_carriers = [];
    pilot_conjugate_carriers = [];
    data_matrix = complex_carrier_matrix;
    pilot_matrix = [];
end

%------------------------------------------------------------------------------
% 步骤3: 频域子载波映射（构造埃尔米特共轭对称，使IFFT输出为实信号）
%------------------------------------------------------------------------------
fprintf('步骤3: 频域子载波映射...\n');
IFFT_modulation = zeros(symbols_per_carrier, IFFT_bin_length);

% 映射数据子载波
IFFT_modulation(:, data_carriers) = data_matrix;
IFFT_modulation(:, data_conjugate_carriers) = conj(data_matrix);

% 映射导频子载波（如果使用导频）
if use_pilot_equalization && ~isempty(pilot_matrix)
    IFFT_modulation(:, pilot_carriers) = pilot_matrix;
    IFFT_modulation(:, pilot_conjugate_carriers) = conj(pilot_matrix);
end
fprintf('  - 数据子载波位置: %d 个\n', length(data_carriers));
if use_pilot_equalization
    fprintf('  - 导频子载波位置: %d 个\n', length(pilot_carriers));
end
fprintf('  完成！\n\n');

%------------------------------------------------------------------------------
% 步骤4: IFFT变换（频域 → 时域）
%------------------------------------------------------------------------------
fprintf('步骤4: IFFT变换（频域 → 时域）...\n');
signal_after_IFFT = ifft(IFFT_modulation, IFFT_bin_length, 2);
time_wave_matrix = signal_after_IFFT;
fprintf('  - IFFT点数: %d\n', IFFT_bin_length);
fprintf('  完成！\n\n');

%------------------------------------------------------------------------------
% 步骤5: 添加循环前缀(CP)与后缀
%------------------------------------------------------------------------------
fprintf('步骤5: 添加循环前缀和后缀...\n');
XX = zeros(symbols_per_carrier, IFFT_bin_length+GI+GIP);
for k = 1:symbols_per_carrier
    % 符号主体部分（中间）
    for i = 1:IFFT_bin_length
        XX(k, i+GI) = signal_after_IFFT(k, i);
    end
    % 循环前缀：将符号尾部复制到开头
    for i = 1:GI
        XX(k, i) = signal_after_IFFT(k, i+IFFT_bin_length-GI);
    end
    % 后缀：将符号头部复制到末尾（用于窗的右侧过渡）
    for j = 1:GIP
        XX(k, IFFT_bin_length+GI+j) = signal_after_IFFT(k, j);
    end
end
time_wave_matrix_cp = XX;
fprintf('  - CP长度: %d\n', GI);
fprintf('  - 后缀长度: %d\n', GIP);
fprintf('  - 符号总长度: %d\n', IFFT_bin_length+GI+GIP);
fprintf('  完成！\n\n');

%------------------------------------------------------------------------------
% 步骤6: OFDM符号加窗处理（LTE风格）
%------------------------------------------------------------------------------
fprintf('步骤6: 升余弦加窗处理（LTE风格）...\n');
windowed_time_wave_matrix_cp = zeros(symbols_per_carrier, IFFT_bin_length+GI+GIP);

% 生成升余弦窗函数（覆盖CP+主体部分，长度为N+GI）
rcos_win_full = rcoswindow(beta, IFFT_bin_length+GI);  % 列向量
rcos_win = rcos_win_full(1:IFFT_bin_length+GI)';  % 只取前N+GI个元素，转置为行向量

% 对每个符号加窗
for i = 1:symbols_per_carrier
    % 符号结构：[CP(GI) | 主体(N) | 后缀(GIP)]
    % 窗函数应用于前 N+GI 个样本（CP+主体）
    windowed_time_wave_matrix_cp(i, 1:IFFT_bin_length+GI) = ...
        real(time_wave_matrix_cp(i, 1:IFFT_bin_length+GI)) .* rcos_win;
    
    % 后缀部分（GIP个样本）：保持原值，用于与下一个符号的CP重叠
    if GIP > 0
        windowed_time_wave_matrix_cp(i, IFFT_bin_length+GI+1:IFFT_bin_length+GI+GIP) = ...
            real(time_wave_matrix_cp(i, IFFT_bin_length+GI+1:IFFT_bin_length+GI+GIP));
    end
end
fprintf('  - 窗函数滚降系数: 1/%d\n', round(1/beta));
fprintf('  完成！\n\n');

%------------------------------------------------------------------------------
% 步骤7: 生成发送信号，并串变换（按帧组织，LTE风格重叠相加）
%------------------------------------------------------------------------------
fprintf('步骤7: 生成发送信号（LTE风格重叠相加）...\n');
% LTE原理：符号拼接时在CP范围内重叠相加，实现平滑过渡
num_frames = total_symbols / symbols_per_frame;
frame_len_CP_suffix = symbols_per_frame*(IFFT_bin_length+GI)+GIP; % 每帧串行长度（仅末尾一次后缀）

% 按帧构造：LTE风格重叠相加
Tx_data = zeros(1, num_frames*frame_len_CP_suffix);

write_offset = 0;
for f = 1:num_frames
    sym_start = (f-1)*symbols_per_frame + 1;
    sym_end   = f*symbols_per_frame;

    % 当前帧的加窗符号矩阵
    frame_windowed = windowed_time_wave_matrix_cp(sym_start:sym_end, :);

    % LTE风格重叠相加：构造帧级串行信号
    frame_serial_windowed = zeros(1, frame_len_CP_suffix);
    
    % 第一个符号：完整写入（包含CP+主体+后缀）
    frame_serial_windowed(1:IFFT_bin_length+GI+GIP) = frame_windowed(1, :);
    
    % 后续符号：重叠相加处理
    for i = 1:(symbols_per_frame-1)
        % 下一个符号在串行序列中的起始位置
        next_symbol_start = i*(IFFT_bin_length+GI) + 1;
        next_symbol_end = (i+1)*(IFFT_bin_length+GI) + GIP;
        
        % LTE重叠机制：
        % 重叠区域：当前符号的后缀(GIP)与下一个符号的CP前GIP个样本重叠
        if GIP > 0 && next_symbol_end <= length(frame_serial_windowed)
            % 重叠区域：当前符号的后缀（已写入）+ 下一个符号的CP前GIP个样本（待写入）
            overlap_start = next_symbol_start;
            overlap_end = overlap_start + GIP - 1;
            
            % 当前符号的后缀部分（在frame_serial_windowed中已写入）
            current_suffix = frame_windowed(i, IFFT_bin_length+GI+1:IFFT_bin_length+GI+GIP);
            
            % 下一个符号的CP前GIP个样本
            next_cp_prefix = frame_windowed(i+1, 1:GIP);
            
            % LTE重叠相加：在重叠区域将两个符号的幅度相加
            frame_serial_windowed(overlap_start:overlap_end) = ...
                current_suffix + next_cp_prefix;
            
            % 非重叠部分：写入下一个符号的剩余部分（CP的剩余部分+主体+后缀）
            if overlap_end < next_symbol_end
                non_overlap_start = overlap_end + 1;
                non_overlap_end = next_symbol_end;
                % 下一个符号从GIP+1开始到结尾
                frame_serial_windowed(non_overlap_start:non_overlap_end) = ...
                    frame_windowed(i+1, GIP+1:IFFT_bin_length+GI+GIP);
            end
        else
            % 如果没有后缀（GIP=0），直接写入下一个符号
            if next_symbol_end <= length(frame_serial_windowed)
                frame_serial_windowed(next_symbol_start:next_symbol_end) = frame_windowed(i+1, :);
            end
        end
    end

    % 写入到整帧串行序列
    Tx_data(write_offset+1:write_offset+frame_len_CP_suffix) = frame_serial_windowed;
    write_offset = write_offset + frame_len_CP_suffix;
end

fprintf('  - 发送信号长度: %d 样本\n', length(Tx_data));
fprintf('  - 信号时长: %.4f ms\n', length(Tx_data)/fs*1e3);
fprintf('  完成！\n\n');

%===============================================================================
% 【信号参数总结】
%===============================================================================
fprintf('========================================\n');
fprintf('信号生成完成！\n');
fprintf('========================================\n');
fprintf('发送信号参数：\n');
fprintf('  - 信号长度: %d 样本\n', length(Tx_data));
fprintf('  - 采样频率: %.2f MHz\n', fs/1e6);
fprintf('  - 信号时长: %.4f ms\n', length(Tx_data)/fs*1e3);
fprintf('  - 理论带宽: %.3f MHz\n', carrier_count * subcarrier_spacing / 1e6);
fprintf('  - OFDM符号数: %d\n', total_symbols);
fprintf('  - 每符号长度: %d 样本 (N+GI+GIP)\n', IFFT_bin_length+GI+GIP);
fprintf('========================================\n\n');

%===============================================================================
% 【可视化（可选）】
%===============================================================================
plot_signal = true;  % 是否绘制信号

if plot_signal
    fprintf('正在绘制信号波形和频谱...\n');
    
    % 绘制时域信号（前1000个样本）
    figure('Name', 'OFDM发送信号时域波形', 'Position', [100, 100, 1200, 600]);
    subplot(2, 1, 1);
    plot_samples = min(1000, length(Tx_data));
    plot(1:plot_samples, Tx_data(1:plot_samples), 'b-', 'LineWidth', 1);
    grid on;
    xlabel('样本索引');
    ylabel('幅度');
    title(sprintf('发送信号时域波形（前%d个样本）', plot_samples));
    
    % 绘制频谱
    subplot(2, 1, 2);
    Nfft = 2^nextpow2(length(Tx_data));
    Tx_Fz = fftshift(fft(Tx_data, Nfft));
    f_axis = (-Nfft/2:(Nfft/2-1)) * fs / Nfft;
    
    plot(f_axis/1e6, 20*log10(abs(Tx_Fz) / max(abs(Tx_Fz)) + eps));
    grid on;
    xlabel('频率 (MHz)');
    ylabel('幅度 (dB)');
    title('发送信号频谱（归一化）');
    xlim([-fs/2/1e6, fs/2/1e6]);
    
    fprintf('图形绘制完成！\n\n');
end

fprintf('所有步骤完成！\n');
fprintf('发送信号已保存在变量 Tx_data 中\n');

%===============================================================================
% 【带宽估计】
%===============================================================================
fprintf('\n========================================\n');
fprintf('开始带宽估计...\n');
fprintf('========================================\n\n');

% 设置载波频率（基带信号可以使用0，或者设置一个载波频率用于上变频）
fc = 0;  % 基带信号，载波频率为0（或者可以设置为其他值，如 fc = 5e6）
snr_estimate = 20;  % 用于带宽估计的SNR值（dB）

% 计算理论带宽
B_ideal = carrier_count * subcarrier_spacing;  % 理论带宽 = 300 * 15e3 = 4.5 MHz
fprintf('理论带宽: %.6f Hz (%.3f MHz)\n\n', B_ideal, B_ideal/1e6);

% 处理信号：上变频和加噪声（参考PSD_OFDM.m的处理方式）
fprintf('处理信号（上变频和加噪声）...\n');
sig_processed = real(Tx_data.*exp(1j*2*pi*fc/fs*(0:length(Tx_data)-1)));
sig_processed = awgn(sig_processed, snr_estimate, 'measured');
fprintf('  - SNR: %.1f dB\n', snr_estimate);
fprintf('  完成！\n\n');

% 使用PSD_OFDM.m的方法进行带宽估计
fprintf('使用Welch算法和AR模型法估计带宽...\n');
try
    % 调用PSD_OFDM函数进行带宽估计（k=1表示绘制PSD图）
    [B_welch, B_ar] = PSD_OFDM(Tx_data, fc, fs, snr_estimate, 1);
    
    % 计算估计误差
    error_welch_abs = abs(B_welch - B_ideal);
    error_welch_rel = (error_welch_abs / B_ideal) * 100;
    
    error_ar_abs = abs(B_ar - B_ideal);
    error_ar_rel = (error_ar_abs / B_ideal) * 100;
    
    % 输出估计结果
    fprintf('\n========================================\n');
    fprintf('带宽估计结果\n');
    fprintf('========================================\n');
    fprintf('理论带宽: %.6f Hz (%.3f MHz)\n', B_ideal, B_ideal/1e6);
    fprintf('\n【Welch算法估计结果】\n');
    fprintf('  估计带宽: %.6f Hz (%.3f MHz)\n', B_welch, B_welch/1e6);
    fprintf('  绝对误差: %.6f Hz (%.3f MHz)\n', error_welch_abs, error_welch_abs/1e6);
    fprintf('  相对误差: %.2f%%\n', error_welch_rel);
    
    fprintf('\n【AR模型法估计结果】\n');
    fprintf('  估计带宽: %.6f Hz (%.3f MHz)\n', B_ar, B_ar/1e6);
    fprintf('  绝对误差: %.6f Hz (%.3f MHz)\n', error_ar_abs, error_ar_abs/1e6);
    fprintf('  相对误差: %.2f%%\n', error_ar_rel);
    
    fprintf('\n【方法对比】\n');
    if error_welch_abs < error_ar_abs
        fprintf('  Welch算法估计更准确（误差更小）\n');
        fprintf('  误差差异: %.6f Hz (%.2f%%)\n', ...
            error_ar_abs - error_welch_abs, error_ar_rel - error_welch_rel);
    elseif error_ar_abs < error_welch_abs
        fprintf('  AR模型法估计更准确（误差更小）\n');
        fprintf('  误差差异: %.6f Hz (%.2f%%)\n', ...
            error_welch_abs - error_ar_abs, error_welch_rel - error_ar_rel);
    else
        fprintf('  两种方法估计精度相同\n');
    end
    fprintf('========================================\n\n');
    
catch ME
    fprintf('错误：带宽估计失败！\n');
    fprintf('错误信息：%s\n', ME.message);
    fprintf('可能原因：\n');
    fprintf('  1. PSD_OFDM.m文件不在MATLAB路径中\n');
    fprintf('  2. 缺少必要的函数（Burg.m, Proximate.m等）\n');
    fprintf('  3. 缺少必要的工具箱（Signal Processing Toolbox等）\n');
    fprintf('\n尝试使用内置方法进行带宽估计...\n');
    
    % 如果PSD_OFDM不可用，使用简化的方法
    try
        % 使用Welch算法
        [Pxx2, f1] = pwelch(sig_processed, hanning(100), 55, 4096*2, fs);
        Pxx22 = Pxx2;
        Pxx22 = Pxx22 / min(Pxx22);
        Pxx22 = 10 * log10(Pxx22);
        Pxx22 = Pxx22 - max(Pxx22);
        
        % 找到-3dB点
        L1 = ceil(length(Pxx22) / 2);
        P1 = Pxx22(1:L1, 1);
        P2 = Pxx22(L1:end, 1);
        
        % 找到最接近-3dB的点
        [~, idx1] = min(abs(P1 - (-3)));
        [~, idx2] = min(abs(P2 - (-3)));
        band1 = f1(idx1);
        band2 = f1(idx2 + L1 - 1);
        B_welch_simple = abs(band2 - band1);
        
        fprintf('简化Welch算法估计带宽: %.6f Hz (%.3f MHz)\n', ...
            B_welch_simple, B_welch_simple/1e6);
    catch
        fprintf('简化方法也失败，请检查MATLAB工具箱\n');
    end
end

fprintf('带宽估计完成！\n');
