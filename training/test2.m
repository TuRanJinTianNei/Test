%% OFDM信号生成与功率谱密度分析
% 使用 Welch 周期图法计算 OFDM 信号的功率谱密度（PSD）
% 通过分段平均对信号进行平滑，降低噪声影响

clear all;
close all;
clc;

%% 参数设置
N = 64;              % 子载波数量
N_cp = 16;           % 循环前缀长度
M = 4;               % QAM调制阶数（4-QAM）
num_symbols = 100;   % OFDM符号数量
fs = 1e6;            % 采样频率 (Hz)

%% 生成OFDM信号
% 1. 生成随机数据并进行QAM调制
data = randi([0 M-1], N, num_symbols);
modulated_data = qammod(data, M, 'gray');

% 2. IFFT变换（将频域信号转换为时域）
ofdm_symbols = ifft(modulated_data, N);

% 3. 添加循环前缀
ofdm_symbols_cp = [ofdm_symbols(end-N_cp+1:end, :); ofdm_symbols];

% 4. 串行化，生成完整的OFDM信号
ofdm_signal = reshape(ofdm_symbols_cp, [], 1);

% 添加噪声（可选，模拟实际信道）
SNR_dB = 20;  % 信噪比（dB）
ofdm_signal_noisy = awgn(ofdm_signal, SNR_dB, 'measured');

%% 使用Welch周期图法计算功率谱密度
% Welch方法通过分段平均实现平滑，降低方差
% 参数说明：
%   - 窗口长度：使用默认值（信号长度的1/8）
%   - 重叠率：默认50%
%   - FFT点数：默认值
%   - 采样频率：fs

% 方法1：使用默认参数
[psd_default, freq_default] = pwelch(ofdm_signal_noisy, [], [], [], fs);

% 方法2：自定义窗口长度和重叠率（更平滑）
window_length = length(ofdm_signal_noisy) / 8;  % 窗口长度为信号长度的1/8
overlap_ratio = 0.5;  % 50%重叠
noverlap = round(window_length * overlap_ratio);
nfft = 2^nextpow2(window_length);  % FFT点数，取2的幂次

[psd_smooth, freq_smooth] = pwelch(ofdm_signal_noisy, window_length, noverlap, nfft, fs);

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
t = (0:length(ofdm_signal_noisy)-1) / fs * 1000;  % 时间轴（毫秒）
plot(t(1:min(1000, length(t))), real(ofdm_signal_noisy(1:min(1000, length(ofdm_signal_noisy)))));
xlabel('时间 (ms)');
ylabel('幅度');
title('OFDM时域信号（实部，前1000个采样点）');
grid on;

% 子图2：Welch方法计算的PSD（默认参数）
subplot(2, 2, 2);
semilogy(freq_default / 1e6, psd_default);
xlabel('频率 (MHz)');
ylabel('功率谱密度 (dB/Hz)');
title('Welch周期图法 - 默认参数');
grid on;

% 子图3：Welch方法计算的PSD（自定义参数，更平滑）
subplot(2, 2, 3);
semilogy(freq_smooth / 1e6, psd_smooth);
xlabel('频率 (MHz)');
ylabel('功率谱密度 (dB/Hz)');
title('Welch周期图法 - 自定义参数（平滑）');
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
fprintf('=== OFDM信号功率谱密度分析结果 ===\n');
fprintf('信号长度: %d 个采样点\n', length(ofdm_signal_noisy));
fprintf('采样频率: %.2f MHz\n', fs / 1e6);
fprintf('Welch PSD点数: %d\n', length(psd_smooth));
fprintf('分段平均后点数: %d\n', length(psd_smoothed));
fprintf('平均功率: %.2f dB\n', 10*log10(mean(psd_smooth)));
fprintf('峰值功率: %.2f dB\n', 10*log10(max(psd_smooth)));
fprintf('===================================\n');

