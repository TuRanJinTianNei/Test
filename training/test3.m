%% 详细仿真小波变换（Wavelet Transform）
% 本脚本详细演示连续小波变换（CWT）和离散小波变换（DWT）的原理、实现和应用

clear; close all; clc;

%% ==================== 第一部分：小波变换基本原理 ====================
fprintf('========== 小波变换基本原理 ==========\n');
fprintf('小波变换通过平移和缩放母小波函数来分析信号的时频特性\n');
fprintf('连续小波变换（CWT）公式：\n');
fprintf('  W(a,b) = (1/√a) ∫ x(t) ψ*((t-b)/a) dt\n');
fprintf('  其中 a 是尺度参数，b 是平移参数\n');
fprintf('离散小波变换（DWT）使用二进制的尺度和平移\n');
fprintf('====================================\n\n');

%% ==================== 第二部分：生成测试信号 ====================
fs = 1000;          % 采样频率 (Hz)
duration = 4;       % 信号时长 (秒)
t = 0:1/fs:duration-1/fs;
N = length(t);

% 信号1：多频率成分信号
f1 = 10; f2 = 50; f3 = 100;
signal1 = sin(2*pi*f1*t) + 0.5*sin(2*pi*f2*t) + 0.3*sin(2*pi*f3*t);

% 信号2：频率随时间变化的信号（chirp）
f_start = 5; f_end = 200;
signal2 = chirp(t, f_start, duration, f_end, 'linear');

% 信号3：包含突变点的信号
signal3 = sin(2*pi*30*t);
signal3(1000:1100) = signal3(1000:1100) + 2;  % 添加突变
signal3(2000:2100) = signal3(2000:2100) - 1.5;  % 添加突变

% 信号4：分段不同频率的信号
signal4 = zeros(size(t));
signal4(t < duration/4) = sin(2*pi*20*t(t < duration/4));
signal4(t >= duration/4 & t < duration/2) = sin(2*pi*80*t(t >= duration/4 & t < duration/2));
signal4(t >= duration/2 & t < 3*duration/4) = sin(2*pi*150*t(t >= duration/2 & t < 3*duration/4));
signal4(t >= 3*duration/4) = sin(2*pi*30*t(t >= 3*duration/4));

% 信号5：包含噪声的信号
signal5 = signal1 + 0.2*randn(size(signal1));

% 选择主要测试信号
test_signal = signal3;  % 可以改为其他信号
signal_name = '包含突变点的信号';

%% ==================== 第三部分：连续小波变换（CWT） ====================
fprintf('========== 连续小波变换（CWT）分析 ==========\n');

% 可用的小波基（CWT专用）
% 注意：
%   - 'morl': Morlet小波（最常用）
%   - 'cmor1-1.5': 复Morlet小波，格式为 'cmor带宽-中心频率'
%   - 'cgau1-1': 复高斯小波，格式为 'cgau阶数-参数'
%   - 'mexh': Mexican Hat小波
%   - 'gaus1': 高斯小波，格式为 'gaus阶数'（1-8）
wavelets = {'morl', 'cmor1-1.5', 'cgau1-1', 'mexh', 'gaus1'};
wavelet_names = {'Morlet', '复Morlet', '复高斯', 'Mexican Hat', '高斯'};

% 选择小波基
selected_wavelet = 'morl';  % Morlet小波，最常用
selected_idx = find(strcmp(wavelets, selected_wavelet));

% 计算CWT
scales = 1:1:128;  % 尺度范围（对应不同的频率）
coefs = cwt(test_signal, scales, selected_wavelet, 1/fs);

% 将尺度转换为频率
freqs = scal2frq(scales, selected_wavelet, 1/fs);

fprintf('使用小波基: %s\n', wavelet_names{selected_idx});
fprintf('尺度范围: %d - %d\n', min(scales), max(scales));
fprintf('对应频率范围: %.2f - %.2f Hz\n', min(freqs), max(freqs));
fprintf('==========================================\n\n');

%% ==================== 第四部分：可视化CWT结果 ====================
% 优化窗口布局
screen_size = get(0, 'ScreenSize');
fig_width = min(1600, screen_size(3) - 100);
fig_height = min(1000, screen_size(4) - 150);
figure('Position', [50, 50, fig_width, fig_height], 'Color', 'white');

% 子图1：原始时域信号
subplot(3, 3, 1);
plot(t, test_signal, 'LineWidth', 1.5, 'Color', [0.2 0.4 0.8]);
xlabel('时间 (秒)', 'FontSize', 11);
ylabel('幅度', 'FontSize', 11);
title(['原始信号: ', signal_name], 'FontSize', 12, 'FontWeight', 'bold');
grid on;
xlim([0, duration]);
set(gca, 'FontSize', 10);

% 子图2：CWT时频图（幅度）
subplot(3, 3, 2);
imagesc(t, freqs, abs(coefs));
axis xy;
cb = colorbar;
cb.Label.String = '幅度';
cb.Label.FontSize = 10;
xlabel('时间 (秒)', 'FontSize', 11);
ylabel('频率 (Hz)', 'FontSize', 11);
title(['CWT时频图 - ', wavelet_names{selected_idx}], 'FontSize', 12, 'FontWeight', 'bold');
colormap('jet');
ylim([0, 200]);
set(gca, 'FontSize', 10);

% 子图3：CWT时频图（相位）
subplot(3, 3, 3);
imagesc(t, freqs, angle(coefs));
axis xy;
cb = colorbar;
cb.Label.String = '相位 (rad)';
cb.Label.FontSize = 10;
xlabel('时间 (秒)', 'FontSize', 11);
ylabel('频率 (Hz)', 'FontSize', 11);
title('CWT相位图', 'FontSize', 12, 'FontWeight', 'bold');
colormap('hsv');
ylim([0, 200]);
set(gca, 'FontSize', 10);

% 子图4-6：不同小波基的CWT对比
for i = 1:min(3, length(wavelets))
    wav = wavelets{i};
    wav_name = wavelet_names{i};
    
    try
        coefs_wav = cwt(test_signal, scales, wav, 1/fs);
        freqs_wav = scal2frq(scales, wav, 1/fs);
        
        subplot(3, 3, 3+i);
        imagesc(t, freqs_wav, abs(coefs_wav));
        axis xy;
        cb = colorbar;
        cb.Label.String = '幅度';
        cb.Label.FontSize = 9;
        xlabel('时间 (秒)', 'FontSize', 10);
        ylabel('频率 (Hz)', 'FontSize', 10);
        title(['CWT - ', wav_name], 'FontSize', 11, 'FontWeight', 'bold');
        colormap('jet');
        ylim([0, 200]);
        set(gca, 'FontSize', 9);
    catch ME
        subplot(3, 3, 3+i);
        text(0.5, 0.5, sprintf('无法使用 %s', wav_name), ...
            'HorizontalAlignment', 'center', 'FontSize', 12);
        title(['CWT - ', wav_name, ' (错误)'], 'FontSize', 11);
        fprintf('警告: 无法使用小波基 %s (%s)\n', wav_name, ME.message);
    end
end

% 子图7：不同尺度的小波系数
subplot(3, 3, 7);
selected_scales = [10, 30, 50, 80];
colors_scale = lines(length(selected_scales));  % 使用不同颜色
hold on;
for i = 1:length(selected_scales)
    scale_idx = find(scales == selected_scales(i), 1);
    if ~isempty(scale_idx)
        plot(t, real(coefs(scale_idx, :)), 'LineWidth', 1.5, ...
            'Color', colors_scale(i, :), ...
            'DisplayName', sprintf('尺度=%d (%.1f Hz)', selected_scales(i), freqs(scale_idx)));
    end
end
xlabel('时间 (秒)', 'FontSize', 11);
ylabel('小波系数', 'FontSize', 11);
title('不同尺度的小波系数', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
grid on;
hold off;
set(gca, 'FontSize', 10);

% 子图8：小波能量分布
subplot(3, 3, 8);
energy = sum(abs(coefs).^2, 2);  % 每个尺度的总能量
plot(freqs, energy, 'LineWidth', 2, 'Color', [0.8 0.2 0.2]);
xlabel('频率 (Hz)', 'FontSize', 11);
ylabel('能量', 'FontSize', 11);
title('小波能量分布', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
xlim([0, 200]);
set(gca, 'FontSize', 10);

% 子图9：特定时刻的频谱
subplot(3, 3, 9);
time_idx = round(length(t) * 0.5);  % 选择中间时刻
plot(freqs, abs(coefs(:, time_idx)), 'LineWidth', 2, 'Color', [0.2 0.6 0.8]);
xlabel('频率 (Hz)', 'FontSize', 11);
ylabel('幅度', 'FontSize', 11);
title(sprintf('t=%.2f秒时的频谱', t(time_idx)), 'FontSize', 12, 'FontWeight', 'bold');
grid on;
xlim([0, 200]);
set(gca, 'FontSize', 10);

sgtitle('连续小波变换（CWT）详细分析', 'FontSize', 16, 'FontWeight', 'bold');
% 优化子图间距
set(gcf, 'Units', 'normalized');
set(gcf, 'Position', [0.02, 0.05, 0.96, 0.88]);

%% ==================== 第五部分：离散小波变换（DWT） ====================
fprintf('========== 离散小波变换（DWT）分析 ==========\n');

% 选择小波基
dwt_wavelet = 'db4';  % Daubechies 4小波
decomposition_level = 5;  % 分解层数

% 说明：频率范围是二进制的（每层减半）
% 这是因为DWT每层分解时：
% 1. 使用低通和高通滤波器分离信号
% 2. 对滤波后的信号进行降采样（采样率减半）
% 3. 降采样后，有效频率范围也减半（奈奎斯特定理）
% 4. 因此频率范围是 2的幂次：fs/2, fs/4, fs/8, fs/16, ...

% 执行DWT分解
[c, l] = wavedec(test_signal, decomposition_level, dwt_wavelet);

% 提取各层细节系数和近似系数
details = cell(decomposition_level, 1);
for i = 1:decomposition_level
    details{i} = detcoef(c, l, i);
end
approximation = appcoef(c, l, dwt_wavelet);

fprintf('使用小波基: %s\n', dwt_wavelet);
fprintf('分解层数: %d\n', decomposition_level);
fprintf('原始信号长度: %d 个采样点 (%.2f 秒)\n', length(test_signal), duration);
fprintf('近似系数长度: %d 个采样点 (覆盖整个时间范围 0-%.2f秒)\n', length(approximation), duration);
for i = 1:decomposition_level
    fprintf('第%d层细节系数长度: %d 个采样点 (覆盖整个时间范围 0-%.2f秒)\n', ...
        i, length(details{i}), duration);
end
fprintf('\n重要说明：\n');
fprintf('  • 所有系数（近似+细节）都是完整的时域信号\n');
fprintf('  • 它们在整个时间范围内都有值（0-%.2f秒）\n', duration);
fprintf('  • 长度不同是因为降采样（每层采样率减半）\n');
fprintf('  • 不是"某个时间点"的值，而是"整个时间序列"\n');
fprintf('\n频带分解说明：\n');
fs_nyquist = fs / 2;
fprintf('  • 近似系数 A%d: 0 - %.2f Hz (最低频)\n', decomposition_level, fs_nyquist / 2^decomposition_level);
for i = 1:decomposition_level
    level = decomposition_level - i + 1;
    f_low = fs_nyquist / 2^(level+1);
    f_high = fs_nyquist / 2^level;
    fprintf('  • 细节系数 D%d: %.2f - %.2f Hz (高频频带 %d)\n', level, f_low, f_high, i);
end
fprintf('  • 每个细节系数对应不同的高频频带\n');
fprintf('  • D1频率最高，D%d频率最低（但仍高于近似系数）\n', decomposition_level);
fprintf('==========================================\n\n');

%% ==================== 第六部分：可视化DWT结果 ====================
fig_width = min(1600, screen_size(3) - 100);
fig_height = min(1000, screen_size(4) - 150);
figure('Position', [100, 100, fig_width, fig_height], 'Color', 'white');

% 子图1：原始信号
subplot(decomposition_level+2, 1, 1);
plot(t, test_signal, 'LineWidth', 1.5, 'Color', [0.2 0.4 0.8]);
ylabel('幅度', 'FontSize', 11);
title('原始信号', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
xlim([0, duration]);
set(gca, 'FontSize', 10);
set(gca, 'XTickLabel', []);  % 隐藏x轴标签，最后一个子图显示

% 子图2：近似系数（最低频成分）
subplot(decomposition_level+2, 1, 2);
t_approx = linspace(0, duration, length(approximation));  % 覆盖整个时间范围
plot(t_approx, approximation, 'LineWidth', 1.5, 'Color', [0.8 0.2 0.2]);
ylabel('幅度', 'FontSize', 11);
title(sprintf('近似系数 A%d (长度=%d, 覆盖0-%.2f秒)', ...
    decomposition_level, length(approximation), duration), ...
    'FontSize', 12, 'FontWeight', 'bold');
grid on;
xlim([0, duration]);
set(gca, 'FontSize', 10);
set(gca, 'XTickLabel', []);  % 隐藏x轴标签

% 子图3-7：各层细节系数
colors_detail = lines(decomposition_level);  % 使用不同颜色
fs_nyquist = fs / 2;  % 奈奎斯特频率
for i = 1:decomposition_level
    subplot(decomposition_level+2, 1, i+2);
    t_detail = linspace(0, duration, length(details{i}));  % 覆盖整个时间范围
    plot(t_detail, details{i}, 'LineWidth', 1.5, 'Color', colors_detail(i, :));
    ylabel('幅度', 'FontSize', 11);
    level = decomposition_level - i + 1;
    % 计算该层对应的频率范围
    f_low = fs_nyquist / 2^(level+1);
    f_high = fs_nyquist / 2^level;
    title(sprintf('细节系数 D%d (第%d层, 频率范围: %.1f-%.1f Hz, 长度=%d)', ...
        level, i, f_low, f_high, length(details{i})), ...
        'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    xlim([0, duration]);
    set(gca, 'FontSize', 10);
    if i < decomposition_level
        set(gca, 'XTickLabel', []);  % 隐藏x轴标签，最后一个子图显示
    else
        xlabel('时间 (秒)', 'FontSize', 11);
    end
end

sgtitle(sprintf('离散小波变换（DWT）多分辨率分解 - %s', dwt_wavelet), ...
    'FontSize', 16, 'FontWeight', 'bold');
% 优化子图间距
set(gcf, 'Units', 'normalized');
set(gcf, 'Position', [0.02, 0.05, 0.96, 0.88]);

%% ==================== DWT频带分解可视化 ====================
% 详细展示每个细节系数对应的频率范围
fprintf('========== DWT频带分解可视化 ==========\n');
fprintf('展示每个细节系数对应的频率范围\n');
fprintf('==========================================\n\n');

figure('Position', [250, 250, fig_width, fig_height], 'Color', 'white');

% 计算各频带的频率范围
fs_nyquist = fs / 2;
freq_bands = zeros(decomposition_level + 1, 2);
band_names = cell(decomposition_level + 1, 1);

% 近似系数频率范围
freq_bands(1, 1) = 0;
freq_bands(1, 2) = fs_nyquist / 2^decomposition_level;
band_names{1} = sprintf('A%d', decomposition_level);

% 各层细节系数频率范围
for i = 1:decomposition_level
    level = decomposition_level - i + 1;
    freq_bands(i+1, 1) = fs_nyquist / 2^(level+1);
    freq_bands(i+1, 2) = fs_nyquist / 2^level;
    band_names{i+1} = sprintf('D%d', level);
end

% 子图1：频带分解示意图
subplot(2, 2, 1);
hold on;
colors_band = [0.8 0.2 0.2; lines(decomposition_level)];
for i = 1:size(freq_bands, 1)
    y_pos = size(freq_bands, 1) - i + 1;
    rectangle('Position', [freq_bands(i,1), y_pos-0.4, ...
        freq_bands(i,2)-freq_bands(i,1), 0.8], ...
        'FaceColor', colors_band(i, :), 'EdgeColor', 'k', 'LineWidth', 1.5);
    text((freq_bands(i,1)+freq_bands(i,2))/2, y_pos, ...
        sprintf('%s: %.1f-%.1f Hz', band_names{i}, freq_bands(i,1), freq_bands(i,2)), ...
        'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
end
xlabel('频率 (Hz)', 'FontSize', 11);
ylabel('频带', 'FontSize', 11);
title('DWT频带分解示意图', 'FontSize', 12, 'FontWeight', 'bold');
xlim([0, fs_nyquist]);
ylim([0.5, size(freq_bands, 1) + 0.5]);
set(gca, 'YTick', 1:size(freq_bands, 1));
set(gca, 'YTickLabel', flipud(band_names));
grid on;
set(gca, 'FontSize', 10);

% 子图2：各频带的能量分布
subplot(2, 2, 2);
energy_bands = zeros(size(freq_bands, 1), 1);
energy_bands(1) = sum(approximation.^2);
for i = 1:decomposition_level
    energy_bands(i+1) = sum(details{i}.^2);
end
bar(energy_bands, 'FaceColor', [0.3 0.6 0.9]);
set(gca, 'XTickLabel', band_names);
xlabel('系数类型', 'FontSize', 11);
ylabel('能量', 'FontSize', 11);
title('各频带的能量分布', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 10);

% 子图3：原始信号频谱（FFT）与DWT频带对比
subplot(2, 2, 3);
% 计算原始信号的FFT
N_fft = 2^nextpow2(length(test_signal));
fft_signal = fft(test_signal, N_fft);
f_fft = (0:N_fft/2) * fs / N_fft;
psd = abs(fft_signal(1:N_fft/2+1)).^2;

semilogy(f_fft, psd, 'LineWidth', 1.5, 'Color', [0.2 0.4 0.8], 'DisplayName', '原始信号频谱');
hold on;
% 标记各频带范围
for i = 1:size(freq_bands, 1)
    if i == 1
        color_band = [0.8 0.2 0.2];
    else
        color_band = colors_detail(i-1, :);
    end
    plot([freq_bands(i,1), freq_bands(i,2)], [max(psd)*0.1, max(psd)*0.1], ...
        'LineWidth', 8, 'Color', color_band, 'DisplayName', band_names{i});
end
xlabel('频率 (Hz)', 'FontSize', 11);
ylabel('功率谱密度', 'FontSize', 11);
title('原始信号频谱与DWT频带对比', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
grid on;
xlim([0, min(200, fs_nyquist)]);
set(gca, 'FontSize', 10);

% 子图4：说明文字
subplot(2, 2, 4);
axis off;
text_str = {
    'DWT频带分解说明：';
    '';
    '✓ 每个细节系数对应不同的高频频带';
    '✓ 频率范围是二进制的（每层减半）';
    '';
    '频带划分（fs=1000 Hz）：';
    sprintf('  • A%d: 0 - %.1f Hz (最低频)', decomposition_level, freq_bands(1,2));
};
for i = 1:decomposition_level
    level = decomposition_level - i + 1;
    text_str{end+1} = sprintf('  • D%d: %.1f - %.1f Hz (高频频带 %d)', ...
        level, freq_bands(i+1,1), freq_bands(i+1,2), i);
end
text_str{end+1} = '';
text_str{end+1} = '特点：';
text_str{end+1} = '  • D1: 最高频 (250-500 Hz)';
text_str{end+1} = sprintf('  • D%d: 最低高频 (%.1f-%.1f Hz)', ...
    decomposition_level, freq_bands(end,1), freq_bands(end,2));
text_str{end+1} = '  • 每个频带互不重叠';
text_str{end+1} = '  • 所有频带覆盖完整频谱';
text(0.1, 0.5, text_str, 'FontSize', 11, 'VerticalAlignment', 'middle', ...
    'FontName', 'FixedWidth');

sgtitle('DWT频带分解详细说明', 'FontSize', 14, 'FontWeight', 'bold');
set(gcf, 'Units', 'normalized');
set(gcf, 'Position', [0.02, 0.1, 0.96, 0.8]);

%% ==================== 二进制频率范围详解 ====================
% 详细解释"频率范围是二进制的（每层减半）"的含义
fprintf('========== 二进制频率范围详解 ==========\n');
fprintf('解释"每层减半"的含义和原理\n');
fprintf('==========================================\n\n');

figure('Position', [300, 300, fig_width, fig_height], 'Color', 'white');

fs_nyquist = fs / 2;  % 奈奎斯特频率 = 500 Hz

% 子图1：频率范围减半过程示意图
subplot(2, 3, 1);
hold on;
% 原始信号频率范围
rectangle('Position', [0, 4.5, fs_nyquist, 0.8], 'FaceColor', [0.8 0.8 0.8], ...
    'EdgeColor', 'k', 'LineWidth', 2);
text(fs_nyquist/2, 4.9, sprintf('原始信号: 0 - %.0f Hz', fs_nyquist), ...
    'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');

% 第1层分解
rectangle('Position', [0, 3.5, fs_nyquist/2, 0.8], 'FaceColor', [0.8 0.2 0.2], ...
    'EdgeColor', 'k', 'LineWidth', 1.5);
text(fs_nyquist/4, 3.9, sprintf('A1: 0 - %.0f Hz', fs_nyquist/2), ...
    'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', 'w', 'FontWeight', 'bold');
rectangle('Position', [fs_nyquist/2, 3.5, fs_nyquist/2, 0.8], 'FaceColor', [0.2 0.6 0.8], ...
    'EdgeColor', 'k', 'LineWidth', 1.5);
text(3*fs_nyquist/4, 3.9, sprintf('D1: %.0f - %.0f Hz', fs_nyquist/2, fs_nyquist), ...
    'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', 'w', 'FontWeight', 'bold');

% 第2层分解（对A1继续分解）
rectangle('Position', [0, 2.5, fs_nyquist/4, 0.8], 'FaceColor', [0.8 0.2 0.2], ...
    'EdgeColor', 'k', 'LineWidth', 1.5);
text(fs_nyquist/8, 2.9, sprintf('A2: 0 - %.0f Hz', fs_nyquist/4), ...
    'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', 'w', 'FontWeight', 'bold');
rectangle('Position', [fs_nyquist/4, 2.5, fs_nyquist/4, 0.8], 'FaceColor', [0.2 0.8 0.4], ...
    'EdgeColor', 'k', 'LineWidth', 1.5);
text(3*fs_nyquist/8, 2.9, sprintf('D2: %.0f - %.0f Hz', fs_nyquist/4, fs_nyquist/2), ...
    'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', 'w', 'FontWeight', 'bold');

% 第3层分解（对A2继续分解）
rectangle('Position', [0, 1.5, fs_nyquist/8, 0.8], 'FaceColor', [0.8 0.2 0.2], ...
    'EdgeColor', 'k', 'LineWidth', 1.5);
text(fs_nyquist/16, 1.9, sprintf('A3: 0 - %.0f Hz', fs_nyquist/8), ...
    'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', 'w', 'FontWeight', 'bold');
rectangle('Position', [fs_nyquist/8, 1.5, fs_nyquist/8, 0.8], 'FaceColor', [0.8 0.6 0.2], ...
    'EdgeColor', 'k', 'LineWidth', 1.5);
text(3*fs_nyquist/16, 1.9, sprintf('D3: %.0f - %.0f Hz', fs_nyquist/8, fs_nyquist/4), ...
    'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', 'w', 'FontWeight', 'bold');

% 标注箭头（使用annotation函数标注分解过程）
annotation('arrow', [0.3, 0.3], [0.7, 0.65], 'LineWidth', 2, 'Color', [0.5 0.5 0.5]);
text(fs_nyquist/2 + 20, 3.95, '第1层分解', 'FontSize', 9, 'Rotation', -90);

annotation('arrow', [0.15, 0.15], [0.55, 0.5], 'LineWidth', 2, 'Color', [0.5 0.5 0.5]);
text(fs_nyquist/4 + 15, 2.95, '第2层分解', 'FontSize', 9, 'Rotation', -90);

xlabel('频率 (Hz)', 'FontSize', 11);
ylabel('分解层数', 'FontSize', 11);
title('频率范围减半过程示意图', 'FontSize', 12, 'FontWeight', 'bold');
xlim([0, fs_nyquist]);
ylim([1, 5.5]);
set(gca, 'YTick', []);
grid on;
set(gca, 'FontSize', 10);

% 子图2：频率范围的数学关系
subplot(2, 3, 2);
freq_values = zeros(decomposition_level + 1, 1);
freq_values(1) = fs_nyquist / 2^decomposition_level;  % A5的上限
for i = 1:decomposition_level
    freq_values(i+1) = fs_nyquist / 2^(decomposition_level - i + 1);  % 各层D的上限
end

% 绘制频率值的对数关系
semilogy(0:decomposition_level, [freq_values(1); freq_values(2:end)], ...
    'o-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', [0.3 0.6 0.9]);
xlabel('分解层数', 'FontSize', 11);
ylabel('频率上限 (Hz, 对数刻度)', 'FontSize', 11);
title('频率范围的指数衰减（每层减半）', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
set(gca, 'XTick', 0:decomposition_level);
set(gca, 'XTickLabel', {'A5', 'D5', 'D4', 'D3', 'D2', 'D1'});
set(gca, 'FontSize', 10);

% 子图3：频率范围数值表
subplot(2, 3, 3);
axis off;
table_data = cell(decomposition_level + 2, 3);
table_data{1, 1} = '系数';
table_data{1, 2} = '频率范围 (Hz)';
table_data{1, 3} = '频率比';
table_data{2, 1} = sprintf('A%d', decomposition_level);
table_data{2, 2} = sprintf('0 - %.2f', fs_nyquist / 2^decomposition_level);
table_data{2, 3} = sprintf('1/%d', 2^decomposition_level);
for i = 1:decomposition_level
    level = decomposition_level - i + 1;
    f_low = fs_nyquist / 2^(level+1);
    f_high = fs_nyquist / 2^level;
    table_data{i+2, 1} = sprintf('D%d', level);
    table_data{i+2, 2} = sprintf('%.2f - %.2f', f_low, f_high);
    table_data{i+2, 3} = sprintf('1/%d - 1/%d', 2^(level+1), 2^level);
end
text(0.1, 0.5, table_data, 'FontSize', 10, 'VerticalAlignment', 'middle', ...
    'FontName', 'FixedWidth', 'FontWeight', 'bold');
title('频率范围数值表', 'FontSize', 12, 'FontWeight', 'bold', 'Position', [0.5, 0.95]);

% 子图4：二进制关系可视化
subplot(2, 3, 4);
% 绘制频率范围的二进制关系
freq_ranges = zeros(decomposition_level + 1, 2);
freq_ranges(1, 2) = fs_nyquist / 2^decomposition_level;
for i = 1:decomposition_level
    level = decomposition_level - i + 1;
    freq_ranges(i+1, 1) = fs_nyquist / 2^(level+1);
    freq_ranges(i+1, 2) = fs_nyquist / 2^level;
end

hold on;
for i = 1:size(freq_ranges, 1)
    if i == 1
        color = [0.8 0.2 0.2];
        name = sprintf('A%d', decomposition_level);
    else
        color = colors_detail(i-1, :);
        name = sprintf('D%d', decomposition_level - i + 2);
    end
    % 使用rectangle绘制条形图（兼容性更好）
    width = freq_ranges(i, 2) - freq_ranges(i, 1);
    rectangle('Position', [freq_ranges(i, 1), i-0.4, width, 0.8], ...
        'FaceColor', color, 'EdgeColor', 'k', 'LineWidth', 1.5);
    text((freq_ranges(i,1)+freq_ranges(i,2))/2, i, name, ...
        'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'w');
end
xlabel('频率 (Hz)', 'FontSize', 11);
ylabel('系数', 'FontSize', 11);
title('各系数频率范围（二进制划分）', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'YTick', 1:size(freq_ranges, 1));
set(gca, 'YTickLabel', flipud(band_names));
xlim([0, fs_nyquist]);
grid on;
set(gca, 'FontSize', 10);

% 子图5：每层减半的数学证明
subplot(2, 3, 5);
axis off;
proof_text = {
    '二进制频率范围的数学原理：';
    '';
    '1. 奈奎斯特定理：';
    '   有效频率范围 = 采样率 / 2';
    '';
    '2. DWT每层操作：';
    '   • 低通/高通滤波';
    '   • 降采样（采样率 ÷ 2）';
    '';
    '3. 降采样后的频率范围：';
    '   新采样率 = 原采样率 / 2';
    '   新频率范围 = 新采样率 / 2';
    '              = (原采样率/2) / 2';
    '              = 原频率范围 / 2';
    '';
    '4. 因此频率范围每层减半：';
    sprintf('   第1层: fs/2 = %.0f Hz', fs_nyquist);
    sprintf('   第2层: fs/4 = %.0f Hz', fs_nyquist/2);
    sprintf('   第3层: fs/8 = %.0f Hz', fs_nyquist/4);
    sprintf('   第4层: fs/16 = %.0f Hz', fs_nyquist/8);
    sprintf('   第5层: fs/32 = %.0f Hz', fs_nyquist/16);
    '';
    '5. 这就是"二进制"的含义：';
    '   频率范围 = fs / 2^n';
    '   (n = 1, 2, 3, ...)';
};
text(0.05, 0.5, proof_text, 'FontSize', 10, 'VerticalAlignment', 'middle', ...
    'FontName', 'FixedWidth');

% 子图6：与CWT的对比
subplot(2, 3, 6);
% DWT的二进制频率点
dwt_freq_points = [fs_nyquist / 2^decomposition_level];
for i = 1:decomposition_level
    level = decomposition_level - i + 1;
    dwt_freq_points(end+1) = fs_nyquist / 2^level;
end

% CWT的连续频率范围
cwt_freqs = linspace(1, fs_nyquist, 100);

hold on;
% 绘制CWT的连续频率
plot(cwt_freqs, ones(size(cwt_freqs)), 'b-', 'LineWidth', 3, 'DisplayName', 'CWT: 连续频率');
% 绘制DWT的二进制频率点
plot(dwt_freq_points, ones(size(dwt_freq_points)), 'ro', 'MarkerSize', 10, ...
    'MarkerFaceColor', 'r', 'LineWidth', 2, 'DisplayName', 'DWT: 二进制频率点');
xlabel('频率 (Hz)', 'FontSize', 11);
ylabel('', 'FontSize', 11);
title('DWT vs CWT 频率分辨率对比', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
xlim([0, fs_nyquist]);
ylim([0.5, 1.5]);
grid on;
set(gca, 'YTick', []);
set(gca, 'FontSize', 10);

sgtitle('"频率范围是二进制的（每层减半）"详解', 'FontSize', 14, 'FontWeight', 'bold');
set(gcf, 'Units', 'normalized');
set(gcf, 'Position', [0.02, 0.05, 0.96, 0.88]);

%% ==================== DWT系数时间覆盖说明 ====================
% 添加一个可视化来说明DWT系数覆盖整个时间范围
fprintf('========== DWT系数时间覆盖说明 ==========\n');
fprintf('可视化：展示DWT系数如何覆盖整个时间范围\n');
fprintf('==========================================\n\n');

figure('Position', [200, 200, fig_width, fig_height], 'Color', 'white');

% 选择一个特定的时间点来展示
selected_time = duration / 2;  % 选择中间时刻
time_idx = round(selected_time * fs);

% 子图1：原始信号，标记选定的时间点
subplot(2, 2, 1);
plot(t, test_signal, 'LineWidth', 1.5, 'Color', [0.2 0.4 0.8]);
hold on;
plot([selected_time, selected_time], [min(test_signal), max(test_signal)], ...
    'r--', 'LineWidth', 2, 'DisplayName', sprintf('t=%.2f秒', selected_time));
xlabel('时间 (秒)', 'FontSize', 11);
ylabel('幅度', 'FontSize', 11);
title('原始信号（标记选定时间点）', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best');
grid on;
xlim([0, duration]);
set(gca, 'FontSize', 10);

% 子图2：展示各层系数在选定时间点的值
subplot(2, 2, 2);
% 计算各层系数在选定时间点的值（通过插值）
coeff_values = zeros(decomposition_level + 1, 1);
coeff_names = cell(decomposition_level + 1, 1);

% 近似系数
t_approx = linspace(0, duration, length(approximation));
approx_val = interp1(t_approx, approximation, selected_time, 'linear', 'extrap');
coeff_values(1) = abs(approx_val);
coeff_names{1} = sprintf('A%d', decomposition_level);

% 各层细节系数
for i = 1:decomposition_level
    t_detail = linspace(0, duration, length(details{i}));
    detail_val = interp1(t_detail, details{i}, selected_time, 'linear', 'extrap');
    level = decomposition_level - i + 1;
    coeff_values(i+1) = abs(detail_val);
    coeff_names{i+1} = sprintf('D%d', level);
end

bar(coeff_values, 'FaceColor', [0.3 0.6 0.9]);
set(gca, 'XTickLabel', coeff_names);
xlabel('系数类型', 'FontSize', 11);
ylabel('系数幅度', 'FontSize', 11);
title(sprintf('t=%.2f秒时各层系数的值', selected_time), 'FontSize', 12, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 10);

% 子图3：展示各层系数在整个时间范围内的分布
subplot(2, 2, 3);
hold on;
% 近似系数
t_approx = linspace(0, duration, length(approximation));
plot(t_approx, approximation, 'LineWidth', 2, 'Color', [0.8 0.2 0.2], ...
    'DisplayName', sprintf('A%d (长度=%d)', decomposition_level, length(approximation)));

% 各层细节系数（为了可视化，将它们垂直偏移）
offset = 0;
for i = 1:decomposition_level
    t_detail = linspace(0, duration, length(details{i}));
    level = decomposition_level - i + 1;
    plot(t_detail, details{i} + offset, 'LineWidth', 1.5, ...
        'Color', colors_detail(i, :), ...
        'DisplayName', sprintf('D%d (长度=%d)', level, length(details{i})));
    offset = offset + (max(details{i}) - min(details{i})) * 1.5;
end
xlabel('时间 (秒)', 'FontSize', 11);
ylabel('系数值（已偏移）', 'FontSize', 11);
title('各层系数覆盖整个时间范围', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
grid on;
xlim([0, duration]);
set(gca, 'FontSize', 10);

% 子图4：说明文字
subplot(2, 2, 4);
axis off;
text_str = {
    'DWT系数的时间特性：';
    '';
    '✓ 所有系数都是完整的时域信号';
    '✓ 它们覆盖整个时间范围 (0-4秒)';
    '✓ 不是"某个时间点"的值';
    '';
    '长度不同的原因：';
    '  • 每层分解后降采样（采样率减半）';
    '  • A5: 最低采样率，长度最短';
    '  • D1: 最高采样率，长度最长';
    '';
    '时间对应关系：';
    '  • t_approx = linspace(0, duration, length(approx))';
    '  • t_detail = linspace(0, duration, length(detail))';
    '  • 都映射到 0-4秒 的完整时间范围';
};
text(0.1, 0.5, text_str, 'FontSize', 11, 'VerticalAlignment', 'middle', ...
    'FontName', 'FixedWidth');

sgtitle('DWT系数时间覆盖范围说明', 'FontSize', 14, 'FontWeight', 'bold');
set(gcf, 'Units', 'normalized');
set(gcf, 'Position', [0.02, 0.1, 0.96, 0.8]);

%% ==================== DWT时频表示（小波尺度图）====================
% 说明：DWT的输出虽然看起来像多个频谱，但实际上每个系数都对应时域位置
% 这里构建DWT的时频表示（类似CWT的时频图）
fprintf('========== DWT时频表示说明 ==========\n');
fprintf('DWT输出的是多分辨率系数，不是时频矩阵\n');
fprintf('但可以通过构建小波尺度图来可视化时频关系\n');
fprintf('====================================\n\n');

figure('Position', [150, 150, fig_width, fig_height], 'Color', 'white');

% 构建DWT的时频表示
% 方法：将各层系数按频率范围映射到时频矩阵中
n_time = length(test_signal);
n_freq_bands = decomposition_level + 1;  % 近似系数 + 各层细节系数

% 创建时频矩阵
dwt_tf_matrix = zeros(n_freq_bands, n_time);

% 计算各频带对应的频率范围
fs_nyquist = fs / 2;
freq_bands = cell(n_freq_bands, 1);
freq_bands{1} = sprintf('A%d: 0 - %.1f Hz', decomposition_level, fs_nyquist / 2^decomposition_level);

% 填充近似系数（最低频）
approx_interp = interp1(linspace(0, duration, length(approximation)), approximation, t, 'linear', 'extrap');
dwt_tf_matrix(1, :) = abs(approx_interp);

% 填充各层细节系数
for i = 1:decomposition_level
    level = decomposition_level - i + 1;  % 从高到低
    detail = details{i};
    detail_interp = interp1(linspace(0, duration, length(detail)), detail, t, 'linear', 'extrap');
    dwt_tf_matrix(i+1, :) = abs(detail_interp);
    
    % 计算频率范围
    f_low = fs_nyquist / 2^(level+1);
    f_high = fs_nyquist / 2^level;
    freq_bands{i+1} = sprintf('D%d: %.1f - %.1f Hz', level, f_low, f_high);
end

% 子图1：原始信号
subplot(2, 2, 1);
plot(t, test_signal, 'LineWidth', 1.5, 'Color', [0.2 0.4 0.8]);
xlabel('时间 (秒)', 'FontSize', 11);
ylabel('幅度', 'FontSize', 11);
title('原始信号', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
xlim([0, duration]);
set(gca, 'FontSize', 10);

% 子图2：DWT时频表示（类似CWT）
subplot(2, 2, 2);
freq_axis = 1:n_freq_bands;
imagesc(t, freq_axis, dwt_tf_matrix);
axis xy;
cb = colorbar;
cb.Label.String = '系数幅度';
cb.Label.FontSize = 10;
xlabel('时间 (秒)', 'FontSize', 11);
ylabel('频带', 'FontSize', 11);
title('DWT时频表示（小波尺度图）', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'YTick', 1:n_freq_bands);
set(gca, 'YTickLabel', freq_bands);
colormap('jet');
set(gca, 'FontSize', 10);

% 子图3：CWT时频图（对比）
subplot(2, 2, 3);
imagesc(t, freqs, abs(coefs));
axis xy;
cb2 = colorbar;
cb2.Label.String = '幅度';
cb2.Label.FontSize = 10;
xlabel('时间 (秒)', 'FontSize', 11);
ylabel('频率 (Hz)', 'FontSize', 11);
title('CWT时频图（对比）', 'FontSize', 12, 'FontWeight', 'bold');
colormap('jet');
ylim([0, 200]);
set(gca, 'FontSize', 10);

% 子图4：说明
subplot(2, 2, 4);
axis off;
text_str = {
    'DWT vs CWT 输出形式对比：';
    '';
    'CWT输出：';
    '  • 时频联合矩阵 (时间 × 频率)';
    '  • 连续尺度，频率分辨率高';
    '  • 冗余表示';
    '';
    'DWT输出：';
    '  • 多分辨率系数（近似+细节）';
    '  • 每个系数对应时域位置';
    '  • 非冗余，紧凑表示';
    '';
    'DWT时频表示：';
    '  • 通过插值构建时频矩阵';
    '  • 频率分辨率是二进制的';
    '  • 高频时间分辨率高，低频频率分辨率高';
};
text(0.1, 0.5, text_str, 'FontSize', 11, 'VerticalAlignment', 'middle', ...
    'FontName', 'FixedWidth');

sgtitle('DWT时频表示 vs CWT时频图', 'FontSize', 14, 'FontWeight', 'bold');
set(gcf, 'Units', 'normalized');
set(gcf, 'Position', [0.02, 0.1, 0.96, 0.8]);

%% ==================== 第七部分：DWT重构 ====================
% 使用所有系数重构
reconstructed_full = waverec(c, l, dwt_wavelet);

% 只使用近似系数重构（去除细节）
c_approx = c;
c_approx(1:length(c)-length(approximation)) = 0;
reconstructed_approx = waverec(c_approx, l, dwt_wavelet);

% 只使用第1层细节系数重构
c_detail1 = zeros(size(c));
detail1_start = length(c) - length(details{1}) - length(approximation) + 1;
detail1_end = length(c) - length(approximation);
c_detail1(detail1_start:detail1_end) = c(detail1_start:detail1_end);
reconstructed_detail1 = waverec(c_detail1, l, dwt_wavelet);

fig_width = min(1400, screen_size(3) - 100);
fig_height = min(800, screen_size(4) - 150);
figure('Position', [150, 150, fig_width, fig_height], 'Color', 'white');

% 原始信号
subplot(2, 2, 1);
plot(t, test_signal, 'LineWidth', 1.5, 'Color', [0.2 0.4 0.8]);
xlabel('时间 (秒)', 'FontSize', 11);
ylabel('幅度', 'FontSize', 11);
title('原始信号', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
xlim([0, duration]);
set(gca, 'FontSize', 10);

% 完整重构
subplot(2, 2, 2);
plot(t(1:length(reconstructed_full)), reconstructed_full, 'LineWidth', 1.5, 'Color', [0.2 0.6 0.2]);
xlabel('时间 (秒)', 'FontSize', 11);
ylabel('幅度', 'FontSize', 11);
title('完整重构信号', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
xlim([0, duration]);
set(gca, 'FontSize', 10);

% 只使用近似系数重构
subplot(2, 2, 3);
plot(t(1:length(reconstructed_approx)), reconstructed_approx, 'LineWidth', 1.5, 'Color', [0.8 0.2 0.2]);
xlabel('时间 (秒)', 'FontSize', 11);
ylabel('幅度', 'FontSize', 11);
title('仅使用近似系数重构', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
xlim([0, duration]);
set(gca, 'FontSize', 10);

% 重构误差
subplot(2, 2, 4);
error_signal = test_signal(1:length(reconstructed_full)) - reconstructed_full;
plot(t(1:length(error_signal)), error_signal, 'LineWidth', 1.5, 'Color', [0.9 0.6 0.1]);
xlabel('时间 (秒)', 'FontSize', 11);
ylabel('误差', 'FontSize', 11);
title('重构误差', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
xlim([0, duration]);
set(gca, 'FontSize', 10);

sgtitle('DWT重构验证', 'FontSize', 14, 'FontWeight', 'bold');
% 优化子图间距
set(gcf, 'Units', 'normalized');
set(gcf, 'Position', [0.02, 0.1, 0.96, 0.8]);

fprintf('DWT重构误差统计:\n');
fprintf('  最大误差: %.6f\n', max(abs(error_signal)));
fprintf('  均方误差: %.6f\n', mean(error_signal.^2));
fprintf('  信噪比: %.2f dB\n\n', 10*log10(var(test_signal(1:length(reconstructed_full)))/var(error_signal)));

%% ==================== 第八部分：不同小波基的DWT对比 ====================
dwt_wavelets = {'haar', 'db4', 'db8', 'coif2', 'sym4', 'bior2.2'};
dwt_names = {'Haar', 'Daubechies 4', 'Daubechies 8', 'Coiflet 2', 'Symlet 4', 'Biorthogonal 2.2'};

fig_width = min(1600, screen_size(3) - 100);
fig_height = min(1000, screen_size(4) - 150);
figure('Position', [200, 200, fig_width, fig_height], 'Color', 'white');

for i = 1:length(dwt_wavelets)
    wav = dwt_wavelets{i};
    wav_name = dwt_names{i};
    
    try
        [c_wav, l_wav] = wavedec(test_signal, 4, wav);
        recon_wav = waverec(c_wav, l_wav, wav);
        
        subplot(3, 3, i);
        plot(t(1:length(recon_wav)), recon_wav, 'LineWidth', 1.5, 'Color', [0.2 0.4 0.8]);
        xlabel('时间 (秒)', 'FontSize', 10);
        ylabel('幅度', 'FontSize', 10);
        title(sprintf('DWT重构 - %s', wav_name), 'FontSize', 11, 'FontWeight', 'bold');
        grid on;
        xlim([0, duration]);
        set(gca, 'FontSize', 9);
    catch
        subplot(3, 3, i);
        text(0.5, 0.5, sprintf('%s 不支持', wav_name), 'HorizontalAlignment', 'center');
        title(sprintf('DWT - %s', wav_name), 'FontSize', 11);
    end
end

sgtitle('不同小波基的DWT重构对比', 'FontSize', 14, 'FontWeight', 'bold');
% 优化子图间距
set(gcf, 'Units', 'normalized');
set(gcf, 'Position', [0.02, 0.05, 0.96, 0.88]);

%% ==================== 第九部分：不同信号的CWT和DWT对比 ====================
signals = {
    signal1, '多频率成分信号';
    signal2, '线性调频信号';
    signal3, '包含突变点的信号';
    signal4, '分段不同频率信号';
    signal5, '含噪声信号';
};

fig_width = min(1600, screen_size(3) - 100);
fig_height = min(1200, screen_size(4) - 150);
figure('Position', [250, 250, fig_width, fig_height], 'Color', 'white');

for i = 1:size(signals, 1)
    sig = signals{i, 1};
    sig_name = signals{i, 2};
    
    % 时域信号
    subplot(size(signals, 1), 3, 3*i-2);
    plot(t, sig, 'LineWidth', 1.5, 'Color', [0.2 0.4 0.8]);
    xlabel('时间 (秒)', 'FontSize', 10);
    ylabel('幅度', 'FontSize', 10);
    title(['时域: ', sig_name], 'FontSize', 11, 'FontWeight', 'bold');
    grid on;
    xlim([0, duration]);
    set(gca, 'FontSize', 9);
    
    % CWT
    subplot(size(signals, 1), 3, 3*i-1);
    coefs_sig = cwt(sig, scales, selected_wavelet, 1/fs);
    freqs_sig = scal2frq(scales, selected_wavelet, 1/fs);
    imagesc(t, freqs_sig, abs(coefs_sig));
    axis xy;
    cb = colorbar;
    cb.Label.String = '幅度';
    cb.Label.FontSize = 9;
    xlabel('时间 (秒)', 'FontSize', 10);
    ylabel('频率 (Hz)', 'FontSize', 10);
    title(['CWT: ', sig_name], 'FontSize', 11, 'FontWeight', 'bold');
    colormap('jet');
    ylim([0, 200]);
    set(gca, 'FontSize', 9);
    
    % DWT细节系数
    subplot(size(signals, 1), 3, 3*i);
    [c_sig, l_sig] = wavedec(sig, 3, dwt_wavelet);
    details_sig = cell(3, 1);
    for j = 1:3
        details_sig{j} = detcoef(c_sig, l_sig, j);
    end
    % 显示各层细节系数
    hold on;
    colors = {'r', 'g', 'b'};
    for j = 1:3
        t_detail = linspace(0, duration, length(details_sig{j}));
        plot(t_detail, details_sig{j} + (4-j)*2, 'Color', colors{j}, 'LineWidth', 1);
    end
    hold off;
    xlabel('时间 (秒)', 'FontSize', 10);
    ylabel('细节系数', 'FontSize', 10);
    title(['DWT细节: ', sig_name], 'FontSize', 11, 'FontWeight', 'bold');
    legend('D1', 'D2', 'D3', 'Location', 'best', 'FontSize', 8);
    grid on;
    xlim([0, duration]);
    set(gca, 'FontSize', 9);
end

sgtitle('不同信号的小波变换对比', 'FontSize', 14, 'FontWeight', 'bold');
% 优化子图间距
set(gcf, 'Units', 'normalized');
set(gcf, 'Position', [0.02, 0.05, 0.96, 0.88]);

%% ==================== 第十部分：小波去噪应用 ====================
fprintf('========== 小波去噪应用 ==========\n');

% 使用含噪声信号
noisy_signal = signal5;
noise_level = 0.2;

% 方法1：使用DWT去噪
[c_denoise, l_denoise] = wavedec(noisy_signal, 4, 'db4');
% 使用软阈值去噪
sigma = median(abs(c_denoise)) / 0.6745;  % 估计噪声标准差
threshold = sigma * sqrt(2*log(length(noisy_signal)));
c_denoise_thresh = wthresh(c_denoise, 's', threshold);
denoised_dwt = waverec(c_denoise_thresh, l_denoise, 'db4');

% 方法2：使用MATLAB内置去噪函数
denoised_matlab = wdenoise(noisy_signal, 4, 'Wavelet', 'db4', 'DenoisingMethod', 'Bayes');

fig_width = min(1400, screen_size(3) - 100);
fig_height = min(600, screen_size(4) - 150);
figure('Position', [300, 300, fig_width, fig_height], 'Color', 'white');

subplot(2, 2, 1);
plot(t, signal1, 'LineWidth', 1.5, 'Color', [0.2 0.4 0.8]);
xlabel('时间 (秒)', 'FontSize', 11);
ylabel('幅度', 'FontSize', 11);
title('原始干净信号', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
xlim([0, duration]);
set(gca, 'FontSize', 10);

subplot(2, 2, 2);
plot(t, noisy_signal, 'LineWidth', 1.5, 'Color', [0.8 0.2 0.2]);
xlabel('时间 (秒)', 'FontSize', 11);
ylabel('幅度', 'FontSize', 11);
title(sprintf('含噪声信号 (SNR=%.1f dB)', 10*log10(var(signal1)/var(noise_level^2))), ...
    'FontSize', 12, 'FontWeight', 'bold');
grid on;
xlim([0, duration]);
set(gca, 'FontSize', 10);

subplot(2, 2, 3);
plot(t(1:length(denoised_dwt)), denoised_dwt, 'LineWidth', 1.5, 'Color', [0.2 0.8 0.2]);
xlabel('时间 (秒)', 'FontSize', 11);
ylabel('幅度', 'FontSize', 11);
title('DWT去噪结果', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
xlim([0, duration]);
set(gca, 'FontSize', 10);

subplot(2, 2, 4);
plot(t(1:length(denoised_matlab)), denoised_matlab, 'LineWidth', 1.5, 'Color', [0.8 0.2 0.8]);
xlabel('时间 (秒)', 'FontSize', 11);
ylabel('幅度', 'FontSize', 11);
title('MATLAB内置去噪结果', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
xlim([0, duration]);
set(gca, 'FontSize', 10);

sgtitle('小波去噪应用', 'FontSize', 14, 'FontWeight', 'bold');
% 优化子图间距
set(gcf, 'Units', 'normalized');
set(gcf, 'Position', [0.02, 0.15, 0.96, 0.75]);

% 计算去噪性能
snr_original = 10*log10(var(signal1(1:length(denoised_dwt)))/var(noisy_signal(1:length(denoised_dwt))-signal1(1:length(denoised_dwt))));
snr_dwt = 10*log10(var(signal1(1:length(denoised_dwt)))/var(denoised_dwt(1:length(signal1))-signal1(1:length(denoised_dwt))));
snr_matlab = 10*log10(var(signal1(1:length(denoised_matlab)))/var(denoised_matlab(1:length(signal1))-signal1(1:length(denoised_matlab))));

fprintf('去噪性能:\n');
fprintf('  原始信号SNR: %.2f dB\n', snr_original);
fprintf('  DWT去噪后SNR: %.2f dB\n', snr_dwt);
fprintf('  MATLAB去噪后SNR: %.2f dB\n', snr_matlab);
fprintf('====================================\n\n');

%% ==================== 第十一部分：小波变换参数说明 ====================
fprintf('========== 小波变换参数说明 ==========\n');
fprintf('1. 连续小波变换（CWT）:\n');
fprintf('   - 优点：时频分辨率可调，适合分析非平稳信号\n');
fprintf('   - 缺点：计算量大，冗余度高\n');
fprintf('   - 常用小波：Morlet, Mexican Hat, Gaussian\n');

fprintf('\n2. 离散小波变换（DWT）:\n');
fprintf('   - 优点：计算效率高，无冗余，适合压缩和去噪\n');
fprintf('   - 缺点：时频分辨率固定（二进制的）\n');
fprintf('   - 常用小波：Daubechies (dbN), Haar, Coiflet, Symlet\n');

fprintf('\n3. 小波基选择:\n');
fprintf('   - Haar: 最简单，但频域特性差\n');
fprintf('   - Daubechies (dbN): 平衡时频特性，N越大越平滑\n');
fprintf('   - Coiflet: 近似系数和细节系数都有消失矩\n');
fprintf('   - Biorthogonal: 线性相位，适合图像处理\n');

fprintf('\n4. 分解层数选择:\n');
fprintf('   - 层数越多，频率分辨率越高，但计算量增大\n');
fprintf('   - 通常选择3-5层\n');

fprintf('\n5. 小波变换 vs 傅里叶变换:\n');
fprintf('   - 傅里叶变换：全局频率信息，无时间定位\n');
fprintf('   - 短时傅里叶变换：固定时频分辨率\n');
fprintf('   - 小波变换：多分辨率分析，高频时间分辨率高，低频频率分辨率高\n');
fprintf('====================================\n\n');

fprintf('========== 仿真完成 ==========\n');
fprintf('所有图形已生成，请查看各个窗口了解小波变换的详细特性。\n');
fprintf('小波变换在信号处理、图像压缩、去噪等领域有广泛应用。\n\n');

