%% OFDM信号生成与功率谱密度分析
% 使用 Welch 周期图法计算 OFDM 信号的功率谱密度（PSD）
% 通过分段平均对信号进行平滑，降低噪声影响

clear all;
close all;
clc;

%% 参数设置（根据表3-1 OFDM信号带宽估计仿真默认参数）
N = 512;             % 总子载波个数
N_cp = 32;           % CP长度（循环前缀长度）
M = 4;               % 子载波调制方式：QPSK（4-QAM）
num_symbols = 24;    % OFDM符号数
fs = 20e6;           % 采样率：20MHz
bandwidth = 12e6;    % 带宽：12MHz
SNR_dB = 5;          % 带内信噪比：-5:15 dB（使用中间值5dB，可根据需要调整）
freq_offset = 10e6;  % 频偏：10MHz
wavelet_decomp_level = 5;  % 小波分解层数：5
welch_nfft = 2048;   % Welch谱点数：2048
hamming_window_length = 512;  % 汉明窗长：512

%% 生成OFDM信号
% 1. 生成随机数据并进行QPSK调制（矩形脉冲成形）
data = randi([0 M-1], N, num_symbols);
modulated_data = qammod(data, M, 'gray');  % QPSK调制

% 2. IFFT变换（将频域信号转换为时域，矩形脉冲成形）
ofdm_symbols = ifft(modulated_data, N);

% 3. 添加循环前缀（CP）
ofdm_symbols_cp = [ofdm_symbols(end-N_cp+1:end, :); ofdm_symbols];

% 4. 串行化，生成完整的OFDM信号
ofdm_signal = reshape(ofdm_symbols_cp, [], 1);

% 5. 添加频偏（频率偏移）
t = (0:length(ofdm_signal)-1) / fs;
ofdm_signal = ofdm_signal .* exp(1j * 2 * pi * freq_offset * t(:));

% 6. 添加AWGN噪声（信道类型：AWGN）
ofdm_signal_noisy = awgn(ofdm_signal, SNR_dB, 'measured');

%% 使用Welch周期图法计算功率谱密度
% Welch方法通过分段平均实现平滑，降低方差
% 根据表3-1参数设置：
%   - 汉明窗长：512
%   - Welch谱点数（nfft）：2048
%   - 重叠率：50%（默认）
%   - 采样频率：fs = 20MHz

% 使用汉明窗
window_length = hamming_window_length;  % 汉明窗长：512
overlap_ratio = 0.5;  % 50%重叠
noverlap = round(window_length * overlap_ratio);
nfft = welch_nfft;  % Welch谱点数：2048

% 使用汉明窗进行Welch功率谱估计
hamming_win = hamming(window_length);
[psd_smooth, freq_smooth] = pwelch(ofdm_signal_noisy, hamming_win, noverlap, nfft, fs);

%% 分段平均平滑（第一次平滑）
% 对功率谱密度进行进一步的分段平均平滑
% 将频率轴分成若干段，对每段内的PSD值进行平均

num_segments = 50;  % 分段数量
segment_length = floor(length(psd_smooth) / num_segments);
psd_smoothed = zeros(num_segments, 1);
freq_smoothed = zeros(num_segments, 1);

for i = 1:num_segments
    start_idx = (i-1) * segment_length + 1;
    end_idx = min(i * segment_length, length(psd_smooth));
    
    % 对每段内的PSD值进行平均
    psd_smoothed(i) = mean(psd_smooth(start_idx:end_idx));
    freq_smoothed(i) = mean(freq_smooth(start_idx:end_idx));
end

%% 结果可视化
figure('Position', [100, 100, 1200, 800]);

% 子图1：OFDM时域信号
subplot(2, 2, 1);
t_ms = (0:length(ofdm_signal_noisy)-1) / fs * 1000;  % 时间轴（毫秒）
plot(t_ms(1:min(1000, length(t_ms))), real(ofdm_signal_noisy(1:min(1000, length(ofdm_signal_noisy)))));
xlabel('时间 (ms)');
ylabel('幅度');
title('OFDM时域信号（实部，前1000个采样点）');
grid on;

% 子图2：Welch方法计算的PSD（使用汉明窗，表3-1参数）
subplot(2, 2, 2);
semilogy(freq_smooth / 1e6, psd_smooth);
xlabel('频率 (MHz)');
ylabel('功率谱密度 (dB/Hz)');
title(sprintf('Welch周期图法（汉明窗%d点，nfft=%d）', hamming_window_length, welch_nfft));
grid on;

% 子图3：OFDM信号频谱（FFT）
subplot(2, 2, 3);
N_fft = 2^nextpow2(length(ofdm_signal_noisy));
fft_signal = fft(ofdm_signal_noisy, N_fft);
f_fft = (0:N_fft/2) * fs / N_fft;
psd_fft = abs(fft_signal(1:N_fft/2+1)).^2 / N_fft / fs;
semilogy(f_fft / 1e6, psd_fft);
xlabel('频率 (MHz)');
ylabel('功率谱密度');
title('OFDM信号频谱（FFT）');
grid on;

% 子图4：分段平均平滑后的PSD
subplot(2, 2, 4);
semilogy(freq_smoothed / 1e6, psd_smoothed);
xlabel('频率 (MHz)');
ylabel('功率谱密度 (dB/Hz)');
title('分段平均平滑后的PSD');
grid on;

%% 对比图：原始PSD vs 平滑后PSD
figure('Position', [200, 200, 1000, 600]);
semilogy(freq_smooth / 1e6, psd_smooth, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Welch PSD（原始）');
hold on;
semilogy(freq_smoothed / 1e6, psd_smoothed, 'r-', 'LineWidth', 2, 'DisplayName', '分段平均平滑后');
xlabel('频率 (MHz)');
ylabel('功率谱密度 (dB/Hz)');
title('功率谱密度对比：原始 vs 分段平均平滑');
legend('Location', 'best');
grid on;
hold off;

%% 输出统计信息
fprintf('=== OFDM信号功率谱密度分析结果（表3-1参数）===\n');
fprintf('总子载波个数: %d\n', N);
fprintf('OFDM符号数: %d\n', num_symbols);
fprintf('CP长度: %d\n', N_cp);
fprintf('调制方式: QPSK\n');
fprintf('采样率: %.2f MHz\n', fs / 1e6);
fprintf('理论带宽: %.2f MHz\n', bandwidth / 1e6);
fprintf('频偏: %.2f MHz\n', freq_offset / 1e6);
fprintf('信噪比: %.1f dB\n', SNR_dB);
fprintf('信号长度: %d 个采样点\n', length(ofdm_signal_noisy));
fprintf('Welch PSD点数: %d (nfft=%d)\n', length(psd_smooth), welch_nfft);
fprintf('汉明窗长: %d\n', hamming_window_length);
fprintf('分段平均后点数: %d\n', length(psd_smoothed));
fprintf('平均功率: %.2f dB\n', 10*log10(mean(psd_smooth)));
fprintf('峰值功率: %.2f dB\n', 10*log10(max(psd_smooth)));
fprintf('===================================\n');

%% ==================== 基于小波分解的OFDM带宽估计算法 ====================
% 算法步骤：
% 1) 利用Welch周期图法计算OFDM信号的功率谱，通过Welch算法对信号做第一次平滑处理
% 2) 选择与Welch谱点数相对应的小波分解层数，将平滑处理后的OFDM功率谱做小波分解，并以近似系数重构
%    注意：是在平滑处理后的功率谱（Welch平滑后的psd_smooth）上做小波分解
% 3) 计算重构信号各点间的差分值，统计两个最大差分值的位置以及与这两点相邻的次大值位置
% 4) 将这四个位置分为两组，分别代表带宽的起始点和终止点，做线性内插得到最终带宽

fprintf('\n========== 基于小波分解的OFDM带宽估计 ==========\n');

% 步骤1：使用Welch周期图法计算功率谱（已完成）
% Welch周期图法本身已经对信号进行了第一次平滑处理
% 在平滑处理后的OFDM功率谱上做小波分解
psd_input = psd_smooth;  % 使用Welch平滑后的功率谱（平滑处理后的功率谱）
freq_input = freq_smooth;
N_psd = length(psd_input);  % Welch谱点数

fprintf('步骤1完成：Welch功率谱点数 = %d（已平滑处理）\n', N_psd);

% 步骤2：在平滑处理后的功率谱上做小波分解与重构
% 根据表3-1：小波分解层数固定为5
% 注意：这里是在Welch平滑处理后的功率谱（psd_input）上做小波分解
wavelet_name = 'db4';  % 使用Daubechies 4小波
decomposition_level = wavelet_decomp_level;  % 小波分解层数：5（根据表3-1）

fprintf('步骤2：在平滑后的功率谱上做小波分解，层数 = %d (小波基: %s)\n', ...
    decomposition_level, wavelet_name);

% 执行小波分解（在平滑处理后的功率谱上）
[c, l] = wavedec(psd_input, decomposition_level, wavelet_name);

% 提取近似系数（最低频成分）
approximation = appcoef(c, l, wavelet_name);

% 使用近似系数重构信号（只保留近似系数，去除细节系数）
% 方法：使用wrcoef函数，只使用近似系数重构到指定层数
psd_reconstructed = wrcoef('a', c, l, wavelet_name, decomposition_level);

% 确保长度一致（如果长度不匹配，进行插值）
if length(psd_reconstructed) ~= N_psd
    % 使用插值使长度匹配
    psd_reconstructed = interp1(1:length(psd_reconstructed), psd_reconstructed, ...
        linspace(1, length(psd_reconstructed), N_psd), 'linear', 'extrap');
end

fprintf('重构信号长度 = %d\n', length(psd_reconstructed));

% 步骤3：计算差分值并找到关键位置
% 计算重构信号各点间的差分值
% 公式：r(k) = d(k+1) - d(k), k = 1, 2, ..., N-1
diff_values = diff(psd_reconstructed);  % 差分值
N_diff = length(diff_values);

% 找到两个最大差分值的位置（正的最大值和负的最大值，即最小值的绝对值）
% 正的最大值对应带宽起始点（功率上升），负的最大值对应带宽终止点（功率下降）
[~, idx_max_pos] = max(diff_values);  % 最大正值位置（起始点）
[~, idx_max_neg] = min(diff_values);  % 最小负值位置（终止点）

% 找到与这两个点相邻的次大值位置
% 对于起始点：找到idx_max_pos相邻位置的次大正值（前一个或后一个位置）
% 对于终止点：找到idx_max_neg相邻位置的次小负值（前一个或后一个位置）

% 起始点组：寻找idx_max_pos相邻的次大正值
% 检查前一个和后一个位置
candidates_start = [];
if idx_max_pos > 1 && diff_values(idx_max_pos-1) > 0
    candidates_start = [candidates_start, idx_max_pos-1];
end
if idx_max_pos < N_diff && diff_values(idx_max_pos+1) > 0
    candidates_start = [candidates_start, idx_max_pos+1];
end

if ~isempty(candidates_start)
    [~, local_idx] = max(diff_values(candidates_start));
    idx_second_max_pos = candidates_start(local_idx);
else
    % 如果相邻位置没有正值，扩大搜索范围
    search_range_start = max(1, idx_max_pos-3):min(N_diff, idx_max_pos+3);
    search_range_start(search_range_start == idx_max_pos) = [];  % 排除最大值位置本身
    diff_start_region = diff_values(search_range_start);
    diff_start_region(diff_start_region <= 0) = -inf;  % 只考虑正值
    [~, local_idx] = max(diff_start_region);
    if ~isinf(diff_start_region(local_idx))
        idx_second_max_pos = search_range_start(local_idx);
    else
        idx_second_max_pos = idx_max_pos;  % 如果找不到，使用最大值位置
    end
end

% 终止点组：寻找idx_max_neg相邻的次小负值
% 检查前一个和后一个位置
candidates_end = [];
if idx_max_neg > 1 && diff_values(idx_max_neg-1) < 0
    candidates_end = [candidates_end, idx_max_neg-1];
end
if idx_max_neg < N_diff && diff_values(idx_max_neg+1) < 0
    candidates_end = [candidates_end, idx_max_neg+1];
end

if ~isempty(candidates_end)
    [~, local_idx] = min(diff_values(candidates_end));
    idx_second_max_neg = candidates_end(local_idx);
else
    % 如果相邻位置没有负值，扩大搜索范围
    search_range_end = max(1, idx_max_neg-3):min(N_diff, idx_max_neg+3);
    search_range_end(search_range_end == idx_max_neg) = [];  % 排除最小值位置本身
    diff_end_region = diff_values(search_range_end);
    diff_end_region(diff_end_region >= 0) = inf;  % 只考虑负值
    [~, local_idx] = min(diff_end_region);
    if ~isinf(diff_end_region(local_idx))
        idx_second_max_neg = search_range_end(local_idx);
    else
        idx_second_max_neg = idx_max_neg;  % 如果找不到，使用最小值位置
    end
end

% 四个关键位置
positions = [idx_max_pos, idx_second_max_pos, idx_max_neg, idx_second_max_neg];
positions = sort(positions);  % 按位置排序

fprintf('步骤3完成：找到四个关键位置\n');
fprintf('  位置1 (起始点最大值): %d, 差分值 = %.6f\n', idx_max_pos, diff_values(idx_max_pos));
fprintf('  位置2 (起始点次大值): %d, 差分值 = %.6f\n', idx_second_max_pos, diff_values(idx_second_max_pos));
fprintf('  位置3 (终止点次小值): %d, 差分值 = %.6f\n', idx_second_max_neg, diff_values(idx_second_max_neg));
fprintf('  位置4 (终止点最小值): %d, 差分值 = %.6f\n', idx_max_neg, diff_values(idx_max_neg));

% 步骤4：线性内插得到最终带宽边界
% 将四个位置分为两组：起始点组和终止点组
% 起始点组：较小的两个位置（对应功率上升）
% 终止点组：较大的两个位置（对应功率下降）

% 确定起始点组和终止点组
start_group = [idx_max_pos, idx_second_max_pos];
start_group = sort(start_group);
end_group = [idx_max_neg, idx_second_max_neg];
end_group = sort(end_group);

% 如果起始点和终止点位置重叠或顺序错误，需要调整
if max(start_group) >= min(end_group)
    % 如果重叠，根据差分值的大小重新分组
    if abs(diff_values(idx_max_pos)) > abs(diff_values(idx_max_neg))
        % 起始点在前
        start_group = [min(positions(1:2))];
        end_group = [max(positions(3:4))];
    else
        % 终止点在前（异常情况，通常不会发生）
        start_group = [min(positions(3:4))];
        end_group = [max(positions(1:2))];
    end
end

% 线性内插
% 对于起始点组：根据最大值和次大值的差别做线性内插
if length(start_group) == 2 && start_group(1) ~= start_group(2)
    % 计算权重（根据差分值的差别）
    diff_start = [diff_values(start_group(1)), diff_values(start_group(2))];
    if abs(diff_start(1) - diff_start(2)) > eps
        weight1 = abs(diff_start(1)) / (abs(diff_start(1)) + abs(diff_start(2)));
        weight2 = abs(diff_start(2)) / (abs(diff_start(1)) + abs(diff_start(2)));
        freq_start = freq_input(start_group(1)) * weight1 + freq_input(start_group(2)) * weight2;
    else
        freq_start = mean(freq_input(start_group));
    end
else
    freq_start = freq_input(start_group(1));
end

% 对于终止点组：根据最小值和次小值的差别做线性内插
if length(end_group) == 2 && end_group(1) ~= end_group(2)
    % 计算权重（根据差分值的差别）
    diff_end = [diff_values(end_group(1)), diff_values(end_group(2))];
    if abs(diff_end(1) - diff_end(2)) > eps
        weight1 = abs(diff_end(1)) / (abs(diff_end(1)) + abs(diff_end(2)));
        weight2 = abs(diff_end(2)) / (abs(diff_end(1)) + abs(diff_end(2)));
        freq_end = freq_input(end_group(1)) * weight1 + freq_input(end_group(2)) * weight2;
    else
        freq_end = mean(freq_input(end_group));
    end
else
    freq_end = freq_input(end_group(1));
end

% 计算估计带宽
estimated_bandwidth = freq_end - freq_start;

fprintf('步骤4完成：线性内插得到带宽边界\n');
fprintf('  起始频率 (a点): %.4f MHz\n', freq_start / 1e6);
fprintf('  终止频率 (b点): %.4f MHz\n', freq_end / 1e6);
fprintf('  估计带宽: %.4f MHz\n', estimated_bandwidth / 1e6);

% 理论带宽（用于对比）
% 根据表3-1：带宽 = 12MHz
theoretical_bandwidth = bandwidth;  % 理论带宽：12MHz（根据表3-1）

fprintf('\n理论带宽: %.4f MHz\n', theoretical_bandwidth / 1e6);
fprintf('估计误差: %.2f%%\n', abs(estimated_bandwidth - theoretical_bandwidth) / theoretical_bandwidth * 100);
fprintf('==========================================\n\n');

%% 可视化带宽估计结果
figure('Position', [300, 300, 1400, 900]);

% 子图1：原始功率谱和重构功率谱
subplot(3, 2, 1);
semilogy(freq_input / 1e6, psd_input, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Welch功率谱');
hold on;
semilogy(freq_input / 1e6, psd_reconstructed, 'r--', 'LineWidth', 2, 'DisplayName', '小波重构功率谱');
xlabel('频率 (MHz)');
ylabel('功率谱密度');
title('步骤2：小波分解重构结果');
legend('Location', 'best');
grid on;
hold off;

% 子图2：差分值
subplot(3, 2, 2);
plot(freq_input(1:end-1) / 1e6, diff_values, 'b-', 'LineWidth', 1.5);
hold on;
plot(freq_input(idx_max_pos) / 1e6, diff_values(idx_max_pos), 'ro', ...
    'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', '起始点最大值');
plot(freq_input(idx_second_max_pos) / 1e6, diff_values(idx_second_max_pos), 'r^', ...
    'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', '起始点次大值');
plot(freq_input(idx_max_neg) / 1e6, diff_values(idx_max_neg), 'go', ...
    'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', '终止点最小值');
plot(freq_input(idx_second_max_neg) / 1e6, diff_values(idx_second_max_neg), 'g^', ...
    'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', '终止点次小值');
xlabel('频率 (MHz)');
ylabel('差分值');
title('步骤3：差分值及关键位置');
legend('Location', 'best');
grid on;
hold off;

% 子图3：带宽估计结果
subplot(3, 2, 3);
semilogy(freq_input / 1e6, psd_input, 'b-', 'LineWidth', 1.5, 'DisplayName', '功率谱');
hold on;
plot([freq_start, freq_start] / 1e6, [min(psd_input), max(psd_input)], ...
    'r--', 'LineWidth', 2, 'DisplayName', sprintf('起始点: %.4f MHz', freq_start/1e6));
plot([freq_end, freq_end] / 1e6, [min(psd_input), max(psd_input)], ...
    'g--', 'LineWidth', 2, 'DisplayName', sprintf('终止点: %.4f MHz', freq_end/1e6));
fill([freq_start, freq_end, freq_end, freq_start] / 1e6, ...
    [min(psd_input), min(psd_input), max(psd_input), max(psd_input)], ...
    'y', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', '估计带宽');
xlabel('频率 (MHz)');
ylabel('功率谱密度');
title(sprintf('步骤4：带宽估计结果 (%.4f MHz)', estimated_bandwidth/1e6));
legend('Location', 'best');
grid on;
hold off;

% 子图4：小波分解系数
subplot(3, 2, 4);
details = cell(decomposition_level, 1);
for i = 1:decomposition_level
    details{i} = detcoef(c, l, i);
end
hold on;
plot(approximation, 'r-', 'LineWidth', 2, 'DisplayName', '近似系数');
for i = 1:min(3, decomposition_level)
    plot(details{i}, 'LineWidth', 1, 'DisplayName', sprintf('细节系数 D%d', i));
end
xlabel('系数索引');
ylabel('系数值');
title('小波分解系数');
legend('Location', 'best');
grid on;
hold off;

% 子图5：对比：理论带宽 vs 估计带宽
subplot(3, 2, 5);
bar([theoretical_bandwidth, estimated_bandwidth] / 1e6, 'FaceColor', [0.3 0.6 0.9]);
set(gca, 'XTickLabel', {'理论带宽', '估计带宽'});
ylabel('带宽 (MHz)');
title('带宽对比');
grid on;
text(1, theoretical_bandwidth/1e6, sprintf('%.4f MHz', theoretical_bandwidth/1e6), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10);
text(2, estimated_bandwidth/1e6, sprintf('%.4f MHz', estimated_bandwidth/1e6), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10);

% 子图6：算法流程说明
subplot(3, 2, 6);
axis off;
text_str = {
    '算法流程：';
    '';
    '1. Welch周期图法计算功率谱';
    '   → 平滑处理';
    '';
    '2. 小波分解功率谱';
    '   → 使用近似系数重构';
    '';
    '3. 计算差分值';
    '   → 找到4个关键位置';
    '';
    '4. 线性内插';
    '   → 得到带宽边界';
    '';
    sprintf('结果：%.4f MHz', estimated_bandwidth/1e6);
    sprintf('误差：%.2f%%', abs(estimated_bandwidth - theoretical_bandwidth) / theoretical_bandwidth * 100);
};
text(0.1, 0.5, text_str, 'FontSize', 11, 'VerticalAlignment', 'middle', ...
    'FontName', 'FixedWidth');

sgtitle('基于小波分解的OFDM信号带宽估计', 'FontSize', 14, 'FontWeight', 'bold');

