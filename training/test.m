%% 详细仿真短时傅里叶变换（STFT, Short-Time Fourier Transform）
% 本脚本详细演示 STFT 的原理、实现和应用

clear; close all; clc;

%% ==================== 第一部分：STFT 基本原理 ====================
fprintf('========== STFT 基本原理 ==========\n');
fprintf('STFT 通过在时域加窗，对每个窗口内的信号进行 FFT，得到时频表示\n');
fprintf('公式：STFT(t, f) = ∫ x(τ)w(τ-t)e^(-j2πfτ)dτ\n');
fprintf('====================================\n\n');

%% ==================== 第二部分：生成测试信号 ====================
fs = 1000;          % 采样频率 (Hz)
duration = 4;       % 信号时长 (秒)
t = 0:1/fs:duration-1/fs;
N = length(t);

% 创建多频率成分的测试信号
% 信号1：频率随时间变化的线性调频信号（chirp）
f1_start = 10;      % 起始频率 (Hz)
f1_end = 100;       % 结束频率 (Hz)
chirp_signal = chirp(t, f1_start, duration, f1_end, 'linear');

% 信号2：分段恒定频率信号
f2_1 = 50;          % 第一段频率 (Hz)
f2_2 = 150;         % 第二段频率 (Hz)
segment_signal = zeros(size(t));
segment_idx = t < duration/2;
segment_signal(segment_idx) = sin(2*pi*f2_1*t(segment_idx));
segment_signal(~segment_idx) = sin(2*pi*f2_2*t(~segment_idx));

% 信号3：多个瞬时频率成分
multi_freq_signal = sin(2*pi*30*t) .* (t < duration/2) + ...
                    sin(2*pi*80*t) .* (t >= duration/2 & t < 3*duration/4) + ...
                    sin(2*pi*120*t) .* (t >= 3*duration/4);

% 信号4：频率跳变信号
jump_signal = sin(2*pi*40*t) .* (t < duration/3) + ...
              sin(2*pi*100*t) .* (t >= duration/3 & t < 2*duration/3) + ...
              sin(2*pi*60*t) .* (t >= 2*duration/3);

% 添加噪声（可选）
noise_level = 0.1;
chirp_signal = chirp_signal + noise_level*randn(size(chirp_signal));
segment_signal = segment_signal + noise_level*randn(size(segment_signal));

%% ==================== 第三部分：不同窗函数对比 ====================
% 窗函数参数
win_length = 256;       % 窗长度
noverlap = win_length/2; % 重叠长度（50%重叠）
nfft = 512;             % FFT 点数

% 不同窗函数
windows = {
    @(n) ones(n,1), '矩形窗 (Rectangular)';
    @(n) hann(n), '汉宁窗 (Hann)';
    @(n) hamming(n), '汉明窗 (Hamming)';
    @(n) blackman(n), '布莱克曼窗 (Blackman)';
    @(n) kaiser(n, 5), '凯泽窗 (Kaiser, β=5)';
};

%% ==================== 第五部分：可视化结果 ====================

% 选择要分析的信号
test_signal = segment_signal;  % 可以改为其他信号
signal_name = '分段恒定频率信号';

figure('Position', [50, 50, 1600, 1000]);

% 子图1：原始时域信号
subplot(3, 3, 1);
plot(t, test_signal);
xlabel('时间 (秒)', 'FontSize', 11);
ylabel('幅度', 'FontSize', 11);
title(['原始信号: ', signal_name], 'FontSize', 12, 'FontWeight', 'bold');
grid on;
xlim([0, duration]);

% 子图2-6：不同窗函数的 STFT
for i = 1:min(5, size(windows, 1))
    window_func = windows{i, 1};
    window_name = windows{i, 2};
    win = window_func(win_length);
    
    [S, f, t_stft] = my_stft(test_signal, fs, win, noverlap, nfft);
    
    subplot(3, 3, i+1);
    imagesc(t_stft, f, 20*log10(abs(S) + eps));
    axis xy;
    colorbar;
    xlabel('时间 (秒)', 'FontSize', 11);
    ylabel('频率 (Hz)', 'FontSize', 11);
    title(['STFT - ', window_name], 'FontSize', 12, 'FontWeight', 'bold');
    colormap('jet');
    caxis([-60, 0]);  % 动态范围 -60dB
    ylim([0, 200]);
end

% 子图7：使用 MATLAB 内置函数对比
subplot(3, 3, 7);
[S_matlab, f_matlab, t_matlab] = spectrogram(test_signal, win_length, noverlap, nfft, fs);
imagesc(t_matlab, f_matlab, 20*log10(abs(S_matlab) + eps));
axis xy;
colorbar;
xlabel('时间 (秒)', 'FontSize', 11);
ylabel('频率 (Hz)', 'FontSize', 11);
title('MATLAB spectrogram 函数', 'FontSize', 12, 'FontWeight', 'bold');
colormap('jet');
caxis([-60, 0]);
ylim([0, 200]);

% 子图8：窗函数形状对比
subplot(3, 3, 8);
hold on;
n = 0:win_length-1;
plot(n, ones(win_length,1), 'LineWidth', 2, 'DisplayName', '矩形窗');
plot(n, hann(win_length), 'LineWidth', 2, 'DisplayName', '汉宁窗');
plot(n, hamming(win_length), 'LineWidth', 2, 'DisplayName', '汉明窗');
plot(n, blackman(win_length), 'LineWidth', 2, 'DisplayName', '布莱克曼窗');
xlabel('样本点', 'FontSize', 11);
ylabel('幅度', 'FontSize', 11);
title('窗函数形状对比', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best');
grid on;
hold off;

% 子图9：时频分辨率权衡
subplot(3, 3, 9);
% 短窗：时间分辨率高，频率分辨率低
win_short = 128;
[S_short, f_short, t_short] = my_stft(test_signal, fs, hann(win_short), win_short/2, nfft);
imagesc(t_short, f_short, 20*log10(abs(S_short) + eps));
axis xy;
colorbar;
xlabel('时间 (秒)', 'FontSize', 11);
ylabel('频率 (Hz)', 'FontSize', 11);
title('短窗 (高时间分辨率)', 'FontSize', 12, 'FontWeight', 'bold');
colormap('jet');
caxis([-60, 0]);
ylim([0, 200]);

sgtitle('短时傅里叶变换 (STFT) 详细仿真', 'FontSize', 16, 'FontWeight', 'bold');

%% ==================== 第六部分：窗长度对分辨率的影响 ====================
figure('Position', [100, 100, 1600, 600]);

win_lengths = [64, 128, 256, 512];
n_wins = length(win_lengths);

for i = 1:n_wins
    wl = win_lengths(i);
    noverlap_wl = wl/2;
    
    [S_res, f_res, t_res] = my_stft(test_signal, fs, hann(wl), noverlap_wl, nfft);
    
    subplot(2, n_wins, i);
    imagesc(t_res, f_res, 20*log10(abs(S_res) + eps));
    axis xy;
    colorbar;
    xlabel('时间 (秒)', 'FontSize', 10);
    ylabel('频率 (Hz)', 'FontSize', 10);
    title(sprintf('窗长度 = %d (时间分辨率)', wl), 'FontSize', 11, 'FontWeight', 'bold');
    colormap('jet');
    caxis([-60, 0]);
    ylim([0, 200]);
    
    % 显示窗函数
    subplot(2, n_wins, i+n_wins);
    plot(0:wl-1, hann(wl), 'LineWidth', 2);
    xlabel('样本点', 'FontSize', 10);
    ylabel('幅度', 'FontSize', 10);
    title(sprintf('窗函数 (长度=%d)', wl), 'FontSize', 11);
    grid on;
    xlim([0, wl-1]);
end

sgtitle('窗长度对时频分辨率的影响', 'FontSize', 14, 'FontWeight', 'bold');

%% ==================== 第七部分：不同信号的 STFT 对比 ====================
signals = {
    chirp_signal, '线性调频信号 (Chirp)';
    segment_signal, '分段恒定频率信号';
    multi_freq_signal, '多频率成分信号';
    jump_signal, '频率跳变信号';
};

figure('Position', [150, 150, 1600, 800]);

for i = 1:size(signals, 1)
    sig = signals{i, 1};
    sig_name = signals{i, 2};
    
    % 时域信号
    subplot(4, 2, 2*i-1);
    plot(t, sig);
    xlabel('时间 (秒)', 'FontSize', 10);
    ylabel('幅度', 'FontSize', 10);
    title(['时域: ', sig_name], 'FontSize', 11, 'FontWeight', 'bold');
    grid on;
    xlim([0, duration]);
    
    % STFT
    [S_sig, f_sig, t_sig] = my_stft(sig, fs, hann(win_length), noverlap, nfft);
    subplot(4, 2, 2*i);
    imagesc(t_sig, f_sig, 20*log10(abs(S_sig) + eps));
    axis xy;
    colorbar;
    xlabel('时间 (秒)', 'FontSize', 10);
    ylabel('频率 (Hz)', 'FontSize', 10);
    title(['STFT: ', sig_name], 'FontSize', 11, 'FontWeight', 'bold');
    colormap('jet');
    caxis([-60, 0]);
    ylim([0, 200]);
end

sgtitle('不同信号的 STFT 时频分析', 'FontSize', 14, 'FontWeight', 'bold');

%% ==================== 第八部分：STFT 参数说明 ====================
fprintf('\n========== STFT 参数说明 ==========\n');
fprintf('1. 窗长度 (Window Length):\n');
fprintf('   - 短窗：时间分辨率高，频率分辨率低\n');
fprintf('   - 长窗：时间分辨率低，频率分辨率高\n');
fprintf('   - 当前设置: %d 样本 (%.3f 秒)\n', win_length, win_length/fs);

fprintf('\n2. 重叠长度 (Overlap):\n');
fprintf('   - 增加重叠可以平滑时频图，减少边界效应\n');
fprintf('   - 当前设置: %d 样本 (%.1f%% 重叠)\n', noverlap, noverlap/win_length*100);

fprintf('\n3. FFT 点数 (NFFT):\n');
fprintf('   - 影响频率分辨率\n');
fprintf('   - 当前设置: %d 点 (频率分辨率: %.2f Hz)\n', nfft, fs/nfft);

fprintf('\n4. 窗函数选择:\n');
fprintf('   - 矩形窗：主瓣窄，旁瓣高（频谱泄漏大）\n');
fprintf('   - 汉宁窗：旁瓣低，主瓣较宽（常用）\n');
fprintf('   - 汉明窗：类似汉宁窗，但旁瓣更低\n');
fprintf('   - 布莱克曼窗：旁瓣最低，主瓣最宽\n');

fprintf('\n5. 时频分辨率权衡:\n');
fprintf('   - 不确定性原理：Δt × Δf ≥ 1/(4π)\n');
fprintf('   - 无法同时获得高时间分辨率和高频率分辨率\n');
fprintf('====================================\n\n');

%% ==================== 第九部分：STFT 逆变换 (ISTFT) ====================
fprintf('========== STFT 逆变换 (ISTFT) ==========\n');
fprintf('STFT 是可逆的，可以通过 ISTFT 重建原始信号\n');
fprintf('实现方法：重叠相加法 (Overlap-Add)\n');
fprintf('========================================\n\n');

% 计算 STFT
[S_recon, f_recon, t_recon] = my_stft(test_signal, fs, hann(win_length), noverlap, nfft);

% 使用自定义 ISTFT 函数进行逆变换
reconstructed = my_istft(S_recon, win_length, noverlap, nfft, fs);

% 可视化重建结果
figure('Position', [200, 200, 1400, 500]);

subplot(1, 3, 1);
plot(t, test_signal);
xlabel('时间 (秒)', 'FontSize', 11);
ylabel('幅度', 'FontSize', 11);
title('原始信号', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
xlim([0, duration]);

subplot(1, 3, 2);
plot(t(1:length(reconstructed)), reconstructed);
xlabel('时间 (秒)', 'FontSize', 11);
ylabel('幅度', 'FontSize', 11);
title('ISTFT 重建信号', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
xlim([0, duration]);

subplot(1, 3, 3);
error_signal = test_signal(1:length(reconstructed)) - reconstructed;
plot(t(1:length(error_signal)), error_signal);
xlabel('时间 (秒)', 'FontSize', 11);
ylabel('误差', 'FontSize', 11);
title('重建误差', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
xlim([0, duration]);

fprintf('重建误差统计:\n');
fprintf('  最大误差: %.6f\n', max(abs(error_signal)));
fprintf('  均方误差: %.6f\n', mean(error_signal.^2));
fprintf('  信噪比: %.2f dB\n\n', 10*log10(var(test_signal(1:length(reconstructed)))/var(error_signal)));

sgtitle('STFT 逆变换重建验证', 'FontSize', 14, 'FontWeight', 'bold');

fprintf('========== 仿真完成 ==========\n');
fprintf('所有图形已生成，请查看各个窗口了解 STFT 的详细特性。\n\n');

%% ==================== 函数定义 ====================

function [S, f, t_stft] = my_stft(x, fs, window, noverlap, nfft)
    % 自定义 STFT 实现
    % 输入：
    %   x - 输入信号
    %   fs - 采样频率
    %   window - 窗函数或窗长度
    %   noverlap - 重叠样本数
    %   nfft - FFT 点数
    % 输出：
    %   S - STFT 系数矩阵 (频率 x 时间)
    %   f - 频率向量
    %   t_stft - 时间向量
    
    if isscalar(window)
        win_length = window;
        window = hann(win_length);  % 默认使用汉宁窗
    else
        win_length = length(window);
    end
    
    hop_size = win_length - noverlap;
    n_frames = floor((length(x) - noverlap) / hop_size);
    
    % 初始化输出
    S = zeros(nfft/2+1, n_frames);
    t_stft = zeros(1, n_frames);
    
    % 对每一帧进行 FFT
    for i = 1:n_frames
        start_idx = (i-1)*hop_size + 1;
        end_idx = start_idx + win_length - 1;
        
        if end_idx > length(x)
            break;
        end
        
        % 加窗
        frame = x(start_idx:end_idx) .* window(:);
        
        % FFT
        X = fft(frame, nfft);
        
        % 只保留正频率部分
        S(:, i) = X(1:nfft/2+1);
        
        % 时间中心点
        t_stft(i) = (start_idx + end_idx) / 2 / fs;
    end
    
    % 频率向量
    f = (0:nfft/2) * fs / nfft;
end

function x_recon = my_istft(S, win_length, noverlap, nfft, fs)
    % 自定义 ISTFT 实现（重叠相加法）
    % 输入：
    %   S - STFT 系数矩阵 (频率 x 时间)
    %   win_length - 窗长度
    %   noverlap - 重叠样本数
    %   nfft - FFT 点数
    %   fs - 采样频率
    % 输出：
    %   x_recon - 重建的时域信号
    
    hop_size = win_length - noverlap;
    n_frames = size(S, 2);
    
    % 重建完整频谱（对称性）
    S_full = zeros(nfft, n_frames);
    S_full(1:nfft/2+1, :) = S;
    % 利用共轭对称性填充负频率部分
    for k = 2:nfft/2
        S_full(nfft-k+2, :) = conj(S(k, :));
    end
    
    % 初始化输出信号
    signal_length = (n_frames - 1) * hop_size + win_length;
    x_recon = zeros(signal_length, 1);
    window_sum = zeros(signal_length, 1);
    
    % 使用汉宁窗（与 STFT 中使用的窗相同）
    window = hann(win_length);
    
    % 重叠相加法重建
    for i = 1:n_frames
        % IFFT
        frame = ifft(S_full(:, i), nfft);
        frame = real(frame(1:win_length));  % 取实部
        
        % 加窗
        frame_windowed = frame .* window;
        
        % 计算位置
        start_idx = (i-1)*hop_size + 1;
        end_idx = start_idx + win_length - 1;
        
        if end_idx > signal_length
            end_idx = signal_length;
            frame_windowed = frame_windowed(1:end_idx-start_idx+1);
        end
        
        % 重叠相加
        x_recon(start_idx:end_idx) = x_recon(start_idx:end_idx) + frame_windowed;
        window_sum(start_idx:end_idx) = window_sum(start_idx:end_idx) + window(1:end_idx-start_idx+1).^2;
    end
    
    % 归一化（补偿窗函数的影响）
    window_sum(window_sum < eps) = 1;  % 避免除零
    x_recon = x_recon ./ window_sum;
end

