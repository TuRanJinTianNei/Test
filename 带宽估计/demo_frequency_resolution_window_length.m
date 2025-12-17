%===============================================================================
% demo_frequency_resolution_window_length.m - 演示窗长度对频率分辨率的影响
% 
% 功能说明:
%   通过仿真演示周期图算法中频率分辨率受窗长度限制的原理
%   核心概念：
%   1. 频率分辨率 Δf = fs / N（N为窗长度）
%   2. 窗长度越大，频率分辨率越高（能区分更接近的频率）
%   3. 窗函数的频域主瓣宽度决定频率分辨率
%
% 仿真内容:
%   - 生成两个频率接近的正弦信号
%   - 使用不同窗长度进行周期图估计
%   - 可视化窗函数的频域特性
%   - 验证频率分辨率公式
%
% 使用示例:
%   demo_frequency_resolution_window_length
%
% 创建日期: 2025.12.16
%===============================================================================

clc;
clear all;
close all;

fprintf('========================================\n');
fprintf('周期图频率分辨率与窗长度关系仿真\n');
fprintf('========================================\n\n');

%===============================================================================
% 参数设置
%===============================================================================
fs = 1000;              % 采样频率 (Hz)
duration = 1.0;          % 信号持续时间 (秒)
t = 0:1/fs:duration-1/fs;  % 时间向量
N_total = length(t);    % 总样本数

% 两个频率接近的正弦信号
f1 = 100;               % 第一个信号频率 (Hz)
f2 = 105;               % 第二个信号频率 (Hz)，频率差为5Hz
amplitude1 = 1.0;       % 第一个信号幅度
amplitude2 = 1.0;       % 第二个信号幅度

fprintf('信号参数：\n');
fprintf('  采样频率: %.1f Hz\n', fs);
fprintf('  信号1频率: %.1f Hz\n', f1);
fprintf('  信号2频率: %.1f Hz\n', f2);
fprintf('  频率差: %.1f Hz\n', abs(f2-f1));
fprintf('  信号长度: %d 个样本\n\n', N_total);

%===============================================================================
% 生成测试信号
%===============================================================================
signal1 = amplitude1 * sin(2*pi*f1*t);
signal2 = amplitude2 * sin(2*pi*f2*t);
signal_combined = signal1 + signal2;  % 两个信号的叠加

fprintf('生成测试信号完成！\n\n');

%===============================================================================
% 测试不同的窗长度
%===============================================================================
window_lengths = [50, 100, 200, 500, 1000];  % 不同的窗长度
n_windows = length(window_lengths);

% 存储结果
results = struct();
results.window_lengths = window_lengths;
results.frequency_resolution = zeros(1, n_windows);
results.mainlobe_width = zeros(1, n_windows);
results.can_resolve = false(1, n_windows);  % 是否能区分两个频率

fprintf('========================================\n');
fprintf('测试不同窗长度的频率分辨率\n');
fprintf('========================================\n');

for i = 1:n_windows
    L = window_lengths(i);
    
    fprintf('\n窗长度 %d: ', L);
    
    % 计算理论频率分辨率
    delta_f = fs / L;
    results.frequency_resolution(i) = delta_f;
    
    % 计算矩形窗的主瓣宽度（第一个零点之间的宽度）
    mainlobe_width = 2 * fs / L;
    results.mainlobe_width(i) = mainlobe_width;
    
    % 判断是否能区分两个频率（频率差应大于频率分辨率）
    frequency_separation = abs(f2 - f1);
    can_resolve = frequency_separation >= delta_f;
    results.can_resolve(i) = can_resolve;
    
    fprintf('频率分辨率 = %.2f Hz, ', delta_f);
    fprintf('主瓣宽度 = %.2f Hz, ', mainlobe_width);
    if can_resolve
        fprintf('✓ 能区分 (频率差 %.1f Hz >= 分辨率 %.2f Hz)', ...
            frequency_separation, delta_f);
    else
        fprintf('✗ 不能区分 (频率差 %.1f Hz < 分辨率 %.2f Hz)', ...
            frequency_separation, delta_f);
    end
end

fprintf('\n\n');

%===============================================================================
% 计算周期图（使用不同窗长度）
%===============================================================================
fprintf('计算周期图...\n');

% 存储周期图结果
Pxx_results = cell(1, n_windows);
f_results = cell(1, n_windows);

for i = 1:n_windows
    L = window_lengths(i);
    
    % 截取信号（加矩形窗）
    signal_windowed = signal_combined(1:L);
    
    % 计算周期图（使用FFT）
    NFFT = max(8192, 2^nextpow2(L));  % FFT点数，至少8192以获得足够的频率采样
    X = fft(signal_windowed, NFFT);
    Pxx = (abs(X).^2) / (L * fs);  % 归一化功率谱密度
    
    % 频率向量（单边谱：0到fs/2）
    f = (0:NFFT/2) * fs / NFFT;
    Pxx_single = Pxx(1:NFFT/2+1);
    Pxx_single(2:end-1) = 2 * Pxx_single(2:end-1);  % 转换为单边谱
    
    % 转换为dB
    Pxx_db = 10*log10(Pxx_single);
    Pxx_db = Pxx_db - max(Pxx_db);  % 归一化到峰值0dB
    
    Pxx_results{i} = Pxx_db;
    f_results{i} = f;
end

fprintf('周期图计算完成！\n\n');

%===============================================================================
% 绘制窗函数的频域特性
%===============================================================================
fprintf('绘制窗函数频域特性...\n');

figure('Name', '窗函数频域特性与频率分辨率', 'Position', [100, 100, 1400, 800]);

% 子图1：矩形窗的频域特性（不同长度）
subplot(2, 3, 1);
hold on;
colors = lines(n_windows);
for i = 1:n_windows
    L = window_lengths(i);
    % 矩形窗的频域响应
    NFFT_window = 8192;
    window_rect = ones(1, L);
    W = fft(window_rect, NFFT_window);
    f_window = (0:NFFT_window/2) * fs / NFFT_window;
    W_single = abs(W(1:NFFT_window/2+1));
    W_single(2:end-1) = 2 * W_single(2:end-1);
    W_db = 20*log10(W_single);
    W_db = W_db - max(W_db);
    
    plot(f_window, W_db, 'Color', colors(i,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('L=%d', L));
end
xlim([0, 200]);
ylim([-60, 5]);
xlabel('频率 (Hz)');
ylabel('幅度 (dB)');
title('矩形窗频域特性（不同长度）');
legend('Location', 'best', 'FontSize', 9);
grid on;
hold off;

% 子图2：主瓣宽度与窗长度的关系
subplot(2, 3, 2);
plot(window_lengths, results.mainlobe_width, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
% 理论曲线
L_theory = 50:10:1000;
mainlobe_theory = 2 * fs ./ L_theory;
plot(L_theory, mainlobe_theory, 'r--', 'LineWidth', 1.5, 'DisplayName', '理论值');
xlabel('窗长度 L');
ylabel('主瓣宽度 (Hz)');
title('主瓣宽度 vs 窗长度');
legend('仿真值', '理论值 (2fs/L)', 'Location', 'best');
grid on;
hold off;

% 子图3：频率分辨率与窗长度的关系
subplot(2, 3, 3);
plot(window_lengths, results.frequency_resolution, 'go-', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
% 理论曲线
delta_f_theory = fs ./ L_theory;
plot(L_theory, delta_f_theory, 'r--', 'LineWidth', 1.5, 'DisplayName', '理论值');
xlabel('窗长度 L');
ylabel('频率分辨率 Δf (Hz)');
title('频率分辨率 vs 窗长度');
legend('仿真值', '理论值 (fs/L)', 'Location', 'best');
grid on;
hold off;

% 子图4-6：不同窗长度下的周期图
for i = 1:min(3, n_windows)
    subplot(2, 3, 3+i);
    L = window_lengths(i);
    plot(f_results{i}, Pxx_results{i}, 'b-', 'LineWidth', 1.5);
    hold on;
    % 获取y轴范围
    y_limits = ylim;
    % 标记两个信号的频率
    plot([f1, f1], y_limits, 'r--', 'LineWidth', 1.5, 'DisplayName', sprintf('f1=%.1f Hz', f1));
    plot([f2, f2], y_limits, 'g--', 'LineWidth', 1.5, 'DisplayName', sprintf('f2=%.1f Hz', f2));
    % 标记频率分辨率
    delta_f = results.frequency_resolution(i);
    plot([f1-delta_f/2, f1+delta_f/2], [y_limits(2)-5, y_limits(2)-5], ...
        'k-', 'LineWidth', 3, 'DisplayName', sprintf('分辨率=%.2f Hz', delta_f));
    
    xlim([80, 130]);
    xlabel('频率 (Hz)');
    ylabel('功率谱密度 (dB)');
    if results.can_resolve(i)
        title(sprintf('周期图 (L=%d, Δf=%.2f Hz) ✓可区分', L, delta_f));
    else
        title(sprintf('周期图 (L=%d, Δf=%.2f Hz) ✗不可区分', L, delta_f));
    end
    legend('Location', 'best', 'FontSize', 9);
    grid on;
    hold off;
end

fprintf('图形绘制完成！\n\n');

%===============================================================================
% 绘制详细对比图
%===============================================================================
fprintf('绘制详细对比图...\n');

figure('Name', '频率分辨率对比：不同窗长度', 'Position', [150, 150, 1600, 1000]);

% 绘制所有窗长度的周期图对比
for i = 1:n_windows
    subplot(2, 3, i);
    L = window_lengths(i);
    plot(f_results{i}, Pxx_results{i}, 'b-', 'LineWidth', 1.5);
    hold on;
    
    % 获取y轴范围
    y_limits = ylim;
    
    % 标记两个信号的频率
    plot([f1, f1], y_limits, 'r--', 'LineWidth', 2, 'DisplayName', sprintf('f1=%.1f Hz', f1));
    plot([f2, f2], y_limits, 'g--', 'LineWidth', 2, 'DisplayName', sprintf('f2=%.1f Hz', f2));
    
    % 标记频率分辨率范围
    delta_f = results.frequency_resolution(i);
    y_pos = y_limits(2) - 3;
    plot([f1-delta_f/2, f1+delta_f/2], [y_pos, y_pos], ...
        'k-', 'LineWidth', 4, 'DisplayName', sprintf('Δf=%.2f Hz', delta_f));
    
    % 标记主瓣宽度
    mainlobe = results.mainlobe_width(i);
    y_pos2 = y_limits(2) - 6;
    plot([f1-mainlobe/2, f1+mainlobe/2], [y_pos2, y_pos2], ...
        'm-', 'LineWidth', 2, 'LineStyle', ':', 'DisplayName', sprintf('主瓣=%.2f Hz', mainlobe));
    
    xlim([80, 130]);
    xlabel('频率 (Hz)', 'FontSize', 11);
    ylabel('功率谱密度 (dB)', 'FontSize', 11);
    
    if results.can_resolve(i)
        title(sprintf('窗长度 L=%d\n频率分辨率 Δf=%.2f Hz ✓可区分', L, delta_f), ...
            'FontSize', 12, 'Color', 'green');
    else
        title(sprintf('窗长度 L=%d\n频率分辨率 Δf=%.2f Hz ✗不可区分', L, delta_f), ...
            'FontSize', 12, 'Color', 'red');
    end
    
    legend('Location', 'best', 'FontSize', 9);
    grid on;
    hold off;
end

% 第6个子图：总结表格
subplot(2, 3, 6);
axis off;
table_text = {'窗长度', '频率分辨率', '主瓣宽度', '频率差', '能否区分'; ...
    'L', 'Δf (Hz)', '主瓣 (Hz)', '|f2-f1|', '结果'};
for i = 1:n_windows
    L = window_lengths(i);
    delta_f = results.frequency_resolution(i);
    mainlobe = results.mainlobe_width(i);
    freq_diff = abs(f2 - f1);
    if results.can_resolve(i)
        result_str = '✓ 能';
    else
        result_str = '✗ 不能';
    end
    table_text{end+1, 1} = sprintf('%d', L);
    table_text{end, 2} = sprintf('%.2f', delta_f);
    table_text{end, 3} = sprintf('%.2f', mainlobe);
    table_text{end, 4} = sprintf('%.1f', freq_diff);
    table_text{end, 5} = result_str;
end

% 显示表格
text(0.1, 0.9, '频率分辨率对比表', 'FontSize', 14, 'FontWeight', 'bold');
y_pos = 0.8;
for i = 1:size(table_text, 1)
    row_text = sprintf('%-8s  %-12s  %-12s  %-10s  %-8s', ...
        table_text{i, 1}, table_text{i, 2}, table_text{i, 3}, ...
        table_text{i, 4}, table_text{i, 5});
    if i == 1 || i == 2
        text(0.1, y_pos, row_text, 'FontSize', 11, 'FontWeight', 'bold');
    else
        text(0.1, y_pos, row_text, 'FontSize', 10);
    end
    y_pos = y_pos - 0.08;
end

% 添加关键公式
formula_text = {
    '关键公式：';
    '';
    '频率分辨率：Δf = fs / L';
    '主瓣宽度：主瓣 ≈ 2fs / L';
    '';
    '结论：';
    '• 窗长度L越大，频率分辨率越高';
    '• 频率差必须 ≥ Δf 才能区分';
    '• FFT点数不影响分辨率，只影响采样密度'
};
y_pos = y_pos - 0.1;
for i = 1:length(formula_text)
    text(0.1, y_pos, formula_text{i}, 'FontSize', 10, ...
        'FontName', 'Courier New');
    y_pos = y_pos - 0.06;
end

fprintf('详细对比图绘制完成！\n\n');

%===============================================================================
% 输出总结
%===============================================================================
fprintf('========================================\n');
fprintf('仿真总结\n');
fprintf('========================================\n');
fprintf('核心结论：\n');
fprintf('  1. 频率分辨率 Δf = fs / L（L为窗长度）\n');
fprintf('  2. 主瓣宽度 ≈ 2fs / L\n');
fprintf('  3. 窗长度越大，频率分辨率越高\n');
fprintf('  4. 能区分两个频率的条件：|f2-f1| ≥ Δf\n');
fprintf('\n');
fprintf('本次仿真参数：\n');
fprintf('  采样频率: %.1f Hz\n', fs);
fprintf('  信号频率: %.1f Hz 和 %.1f Hz（频率差 %.1f Hz）\n', ...
    f1, f2, abs(f2-f1));
fprintf('\n');
fprintf('测试结果：\n');
for i = 1:n_windows
    L = window_lengths(i);
    delta_f = results.frequency_resolution(i);
    if results.can_resolve(i)
        fprintf('  L=%d: Δf=%.2f Hz ✓ 能区分（%.1f Hz ≥ %.2f Hz）\n', ...
            L, delta_f, abs(f2-f1), delta_f);
    else
        fprintf('  L=%d: Δf=%.2f Hz ✗ 不能区分（%.1f Hz < %.2f Hz）\n', ...
            L, delta_f, abs(f2-f1), delta_f);
    end
end
fprintf('\n');
fprintf('========================================\n');
fprintf('仿真完成！\n');
fprintf('========================================\n');

