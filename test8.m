% ============================================================================
% 文件名: test8.m
% 功能: 功率谱密度（PSD）计算仿真 - 单侧PSD示例
% 描述: 
%   1. 生成音频信号（采样率 Fs = 44.1 kHz，长度 N = 1024）
%   2. 计算频率分辨率 Fres = Fs/N
%   3. 进行 FFT 得到 X(k)
%   4. 计算单侧功率谱密度（PSD）：
%      - DC 分量 (k=0): PSD(0) = (1/Fres) * |X(0)|^2
%      - 其他频率 (k=1 到 N/2): PSD(k) = (1/(2*Fres)) * |X(k)|^2
%   5. 可视化时域信号、FFT结果和单侧PSD
% ============================================================================

tic
clc;
clear all;
close all;

%===============================================================================
% 【参数配置】
%===============================================================================
Fs = 44.1e3;                % 采样率：44.1 kHz（音频标准采样率）
N = 1024;                   % 信号长度：1024 个样本
Fres = Fs / N;              % 频率分辨率：44.1 kHz / 1024 ≈ 43.07 Hz

fprintf('==== 功率谱密度（PSD）计算仿真 ====\n');
fprintf('采样率 Fs = %.1f kHz\n', Fs/1000);
fprintf('信号长度 N = %d\n', N);
fprintf('频率分辨率 Fres = %.2f Hz\n', Fres);
fprintf('==========================================\n\n');

%===============================================================================
% 【生成测试信号】
%===============================================================================
% 创建时间向量
t = (0:N-1) / Fs;           % 时间向量（秒）

% 生成多频率合成信号（模拟音频信号）
% 包含多个频率分量，便于观察PSD
f1 = 440;                   % 440 Hz（A4音符）
f2 = 880;                   % 880 Hz（A5音符）
f3 = 1320;                  % 1320 Hz
f4 = 5000;                  % 5000 Hz

% 生成信号：多个正弦波的叠加，加上一些噪声
signal = 2.0 * sin(2*pi*f1*t) + ...
         1.5 * sin(2*pi*f2*t) + ...
         1.0 * sin(2*pi*f3*t) + ...
         0.8 * sin(2*pi*f4*t) + ...
         0.3 * randn(size(t));  % 添加少量噪声

% 归一化信号（可选）
signal = signal / max(abs(signal));

%===============================================================================
% 【步骤1：进行 FFT 得到 X(k)】
%===============================================================================
X = fft(signal, N);         % FFT 变换，得到 X(k)

% 频率向量（双边谱，0 到 Fs）
f_double = (0:N-1) * Fs / N;

% 频率向量（单边谱，0 到 Fs/2）
f_single = (0:N/2) * Fs / N;

%===============================================================================
% 【步骤2：计算单侧功率谱密度（PSD）】
%===============================================================================
% 初始化单侧 PSD 向量（长度为 N/2+1，包含 DC 和正频率）
PSD_single = zeros(1, N/2+1);

% DC 分量 (k = 0): PSD(0) = (1 / Fres) * |X(0)|^2
PSD_single(1) = (1 / Fres) * abs(X(1))^2;  % MATLAB索引从1开始，X(1)对应k=0

% 其他频率 (k = 1 到 N/2): PSD(k) = (1 / (2 * Fres)) * |X(k)|^2
for k = 2:(N/2+1)
    PSD_single(k) = (1 / (2 * Fres)) * abs(X(k))^2;
end

% 注意：对于实信号，负频率分量是正频率的共轭，所以单侧PSD已经包含了所有功率

%===============================================================================
% 【对比：双侧PSD（用于对比）】
%===============================================================================
% 双侧PSD：PSD(k) = (1 / (N * Fres)) * |X(k)|^2
PSD_double = (1 / (N * Fres)) * abs(X).^2;

%===============================================================================
% 【可视化结果】
%===============================================================================
figure(1);
set(gcf, 'Position', [100, 100, 1400, 900]);

% 子图1：时域信号
subplot(3,2,1);
plot(t*1000, signal, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('时间 t (ms)', 'FontSize', 12);
ylabel('幅度', 'FontSize', 12);
title('时域信号（多频率合成）', 'FontSize', 14);
axis([0 max(t)*1000 -1.2 1.2]);
legend('信号 x(t)', 'Location', 'best');

% 子图2：FFT 幅度谱（双边）
subplot(3,2,2);
plot(f_double/1000, abs(X), 'r-', 'LineWidth', 1.5);
grid on;
xlabel('频率 f (kHz)', 'FontSize', 12);
ylabel('|X(k)|', 'FontSize', 12);
title('FFT 幅度谱（双边谱）', 'FontSize', 14);
axis([0 Fs/2000 0 max(abs(X))*1.1]);
legend('|X(k)|', 'Location', 'best');

% 子图3：FFT 幅度谱（单边，0 到 Fs/2）
subplot(3,2,3);
plot(f_single/1000, abs(X(1:N/2+1)), 'g-', 'LineWidth', 1.5);
grid on;
xlabel('频率 f (kHz)', 'FontSize', 12);
ylabel('|X(k)|', 'FontSize', 12);
title('FFT 幅度谱（单边谱，0 到 Fs/2）', 'FontSize', 14);
axis([0 Fs/2000 0 max(abs(X(1:N/2+1)))*1.1]);
legend('|X(k)| (单边)', 'Location', 'best');

% 子图4：单侧功率谱密度（PSD）- 线性刻度
subplot(3,2,4);
plot(f_single/1000, PSD_single, 'b-', 'LineWidth', 2);
grid on;
xlabel('频率 f (kHz)', 'FontSize', 12);
ylabel('PSD (功率/Hz)', 'FontSize', 12);
title('单侧功率谱密度（PSD）- 线性刻度', 'FontSize', 14);
axis([0 Fs/2000 0 max(PSD_single)*1.1]);
legend('单侧 PSD', 'Location', 'best');

% 标记主要频率分量
hold on;
[~, idx1] = min(abs(f_single - f1));
[~, idx2] = min(abs(f_single - f2));
[~, idx3] = min(abs(f_single - f3));
[~, idx4] = min(abs(f_single - f4));
plot(f_single(idx1)/1000, PSD_single(idx1), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
plot(f_single(idx2)/1000, PSD_single(idx2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
plot(f_single(idx3)/1000, PSD_single(idx3), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
plot(f_single(idx4)/1000, PSD_single(idx4), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
hold off;

% 子图5：单侧功率谱密度（PSD）- 对数刻度（dB）
subplot(3,2,5);
plot(f_single/1000, 10*log10(PSD_single + eps), 'b-', 'LineWidth', 2);
grid on;
xlabel('频率 f (kHz)', 'FontSize', 12);
ylabel('PSD (dB/Hz)', 'FontSize', 12);
title('单侧功率谱密度（PSD）- 对数刻度（dB）', 'FontSize', 14);
axis([0 Fs/2000 -20 max(10*log10(PSD_single + eps))*1.1]);
legend('单侧 PSD (dB)', 'Location', 'best');

% 标记主要频率分量
hold on;
plot(f_single(idx1)/1000, 10*log10(PSD_single(idx1) + eps), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
plot(f_single(idx2)/1000, 10*log10(PSD_single(idx2) + eps), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
plot(f_single(idx3)/1000, 10*log10(PSD_single(idx3) + eps), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
plot(f_single(idx4)/1000, 10*log10(PSD_single(idx4) + eps), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
hold off;

% 子图6：双侧PSD vs 单侧PSD对比（仅正频率部分）
subplot(3,2,6);
plot(f_single/1000, PSD_double(1:N/2+1), 'r--', 'LineWidth', 1.5);
hold on;
plot(f_single/1000, PSD_single, 'b-', 'LineWidth', 2);
grid on;
xlabel('频率 f (kHz)', 'FontSize', 12);
ylabel('PSD (功率/Hz)', 'FontSize', 12);
title('双侧PSD vs 单侧PSD对比（正频率部分）', 'FontSize', 14);
legend('双侧 PSD', '单侧 PSD', 'Location', 'best');
axis([0 Fs/2000 0 max([max(PSD_double(1:N/2+1)), max(PSD_single)])*1.1]);
hold off;

%===============================================================================
% 【验证：功率计算】
%===============================================================================
% 时域功率
power_time = sum(signal.^2) / N;

% 频域功率（使用单侧PSD）
% 注意：DC分量 + 2倍正频率分量（因为单侧PSD已经考虑了负频率的功率）
power_freq_single = PSD_single(1) * Fres + 2 * sum(PSD_single(2:end)) * Fres;

% 频域功率（使用双侧PSD）
power_freq_double = sum(PSD_double) * Fres;

%===============================================================================
% 【详细输出：主要频率分量的PSD值】
%===============================================================================
fprintf('==== 主要频率分量的PSD值 ====\n');
fprintf('频率 (Hz)    |X(k)|      PSD (功率/Hz)    PSD (dB/Hz)\n');
fprintf('------------------------------------------------------------\n');
fprintf('DC (0 Hz)    %.4f      %.6f          %.2f\n', ...
    abs(X(1)), PSD_single(1), 10*log10(PSD_single(1) + eps));
fprintf('%.1f Hz      %.4f      %.6f          %.2f\n', ...
    f1, abs(X(idx1)), PSD_single(idx1), 10*log10(PSD_single(idx1) + eps));
fprintf('%.1f Hz      %.4f      %.6f          %.2f\n', ...
    f2, abs(X(idx2)), PSD_single(idx2), 10*log10(PSD_single(idx2) + eps));
fprintf('%.1f Hz      %.4f      %.6f          %.2f\n', ...
    f3, abs(X(idx3)), PSD_single(idx3), 10*log10(PSD_single(idx3) + eps));
fprintf('%.1f Hz      %.4f      %.6f          %.2f\n', ...
    f4, abs(X(idx4)), PSD_single(idx4), 10*log10(PSD_single(idx4) + eps));
fprintf('============================================================\n\n');

%===============================================================================
% 【验证：功率守恒（帕塞瓦尔定理）】
%===============================================================================
fprintf('==== 功率守恒验证（帕塞瓦尔定理）====\n');
fprintf('时域功率：%.6f\n', power_time);
fprintf('频域功率（单侧PSD）：%.6f\n', power_freq_single);
fprintf('频域功率（双侧PSD）：%.6f\n', power_freq_double);
fprintf('相对误差（单侧）：%.4f%%\n', abs(power_time - power_freq_single)/power_time*100);
fprintf('相对误差（双侧）：%.4f%%\n', abs(power_time - power_freq_double)/power_time*100);
fprintf('==========================================\n\n');

%===============================================================================
% 【公式说明】
%===============================================================================
fprintf('==== 单侧PSD计算公式说明 ====\n');
fprintf('频率分辨率：Fres = Fs / N = %.2f Hz\n', Fres);
fprintf('\n单侧PSD计算公式：\n');
fprintf('  DC 分量 (k = 0):\n');
fprintf('    PSD(0) = (1 / Fres) * |X(0)|^2\n');
fprintf('    计算值：PSD(0) = (1 / %.2f) * |X(0)|^2 = %.6f\n', ...
    Fres, PSD_single(1));
fprintf('\n  其他频率 (k = 1 到 %d):\n', N/2);
fprintf('    PSD(k) = (1 / (2 * Fres)) * |X(k)|^2\n');
fprintf('    计算值：PSD(1) = (1 / (2 * %.2f)) * |X(1)|^2 = %.6f\n', ...
    Fres, PSD_single(2));
fprintf('\n注意：\n');
fprintf('  - 单侧PSD只计算正频率部分（0 到 Fs/2）\n');
fprintf('  - DC分量使用不同的归一化因子（1/Fres）\n');
fprintf('  - 其他频率使用归一化因子（1/(2*Fres)），因为负频率功率已包含\n');
fprintf('  - 单侧PSD的总功率 = DC功率 + 2×正频率功率\n');
fprintf('==========================================\n\n');

%===============================================================================
% 【额外可视化：PSD公式验证】
%===============================================================================
figure(2);
set(gcf, 'Position', [200, 200, 1400, 600]);

% 子图1：展示前几个频率bin的PSD计算
k_range = 0:10;  % 显示前11个频率bin（包括DC）
f_range = k_range * Fres;
PSD_range = [PSD_single(1), PSD_single(2:min(11, N/2+1))];

subplot(1,2,1);
bar(k_range, PSD_range, 'b', 'LineWidth', 1.5);
grid on;
xlabel('频率索引 k', 'FontSize', 12);
ylabel('PSD (功率/Hz)', 'FontSize', 12);
title('前11个频率bin的单侧PSD值', 'FontSize', 14);
legend('单侧 PSD', 'Location', 'best');
for k = 1:length(k_range)
    text(k_range(k), PSD_range(k), sprintf('%.4f', PSD_range(k)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 9);
end

% 子图2：展示PSD计算公式的验证
% 手动计算前几个bin的PSD，验证公式
PSD_manual = zeros(size(k_range));
for k = 1:length(k_range)
    if k_range(k) == 0
        % DC分量
        PSD_manual(k) = (1 / Fres) * abs(X(1))^2;
    else
        % 其他频率
        PSD_manual(k) = (1 / (2 * Fres)) * abs(X(k_range(k)+1))^2;
    end
end

subplot(1,2,2);
plot(k_range, PSD_range, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
plot(k_range, PSD_manual, 'r--s', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on;
xlabel('频率索引 k', 'FontSize', 12);
ylabel('PSD (功率/Hz)', 'FontSize', 12);
title('PSD公式验证（前11个频率bin）', 'FontSize', 14);
legend('程序计算', '手动计算', 'Location', 'best');
hold off;

%===============================================================================
% 【采样过程与带宽限制的仿真解释】
%===============================================================================
fprintf('\n==== 采样过程与带宽限制的仿真解释 =====\n');
fprintf('==========================================\n\n');

% 创建连续时间信号（更高采样率，模拟连续信号）
Fs_continuous = 10 * Fs;   % 连续信号的"采样率"（用于模拟连续信号）
N_continuous = 10 * N;      % 连续信号的样本数
t_continuous = (0:N_continuous-1) / Fs_continuous;
T_continuous = max(t_continuous);

% 生成连续时间信号（带宽受限信号）
f_max_signal = Fs / 2 - 1000;  % 信号最大频率（略小于奈奎斯特频率）
signal_continuous = sin(2*pi*f1*t_continuous) + ...
                    0.5*sin(2*pi*f_max_signal*t_continuous);

% 采样过程：在时域用冲击函数采样
% 采样函数：s(t) = Σ δ(t - n*Ts)，其中 Ts = 1/Fs
Ts = 1 / Fs;                % 采样间隔
n_samples = 0:floor(T_continuous/Ts);  % 采样时刻索引
t_samples = n_samples * Ts;  % 采样时刻

% 创建采样冲击函数序列（在采样时刻有值，其他时刻为0）
impulse_train = zeros(size(t_continuous));
sample_indices = round(t_samples * Fs_continuous) + 1;
sample_indices = sample_indices(sample_indices <= length(t_continuous));
impulse_train(sample_indices) = 1;

% 采样后的信号：x_sampled(t) = x(t) * s(t) = x(t) * Σ δ(t - n*Ts)
signal_sampled = signal_continuous .* impulse_train;

% 计算连续信号的频谱
X_continuous = fftshift(fft(signal_continuous, N_continuous)) * (1/Fs_continuous);
f_continuous = (-N_continuous/2:N_continuous/2-1) * Fs_continuous / N_continuous;

% 计算采样后信号的频谱
X_sampled = fftshift(fft(signal_sampled, N_continuous)) * (1/Fs_continuous);

% 计算采样函数的频谱（冲击函数序列的频谱）
S_impulse = fftshift(fft(impulse_train, N_continuous)) * (1/Fs_continuous);

%===============================================================================
% 【可视化：采样过程与时域表示】
%===============================================================================
figure(3);
set(gcf, 'Position', [300, 300, 1400, 1000]);

% 子图1：连续时间信号
subplot(3,2,1);
plot(t_continuous*1000, signal_continuous, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('时间 t (ms)', 'FontSize', 12);
ylabel('幅度', 'FontSize', 12);
title('连续时间信号 x(t)', 'FontSize', 14);
axis([0 min(5, max(t_continuous)*1000) -2 2]);
legend('连续信号', 'Location', 'best');

% 子图2：采样冲击函数序列
subplot(3,2,2);
stem(t_samples(1:min(50, length(t_samples)))*1000, ...
     impulse_train(sample_indices(1:min(50, length(sample_indices)))), ...
     'r', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('时间 t (ms)', 'FontSize', 12);
ylabel('幅度', 'FontSize', 12);
title(sprintf('采样冲击函数序列 s(t) = Σ δ(t - n*Ts), Ts = %.3f ms', Ts*1000), ...
    'FontSize', 14);
axis([0 min(5, max(t_samples)*1000) -0.2 1.2]);
legend('冲击函数', 'Location', 'best');

% 子图3：采样后的信号（时域）
subplot(3,2,3);
plot(t_continuous*1000, signal_continuous, 'b--', 'LineWidth', 1);
hold on;
stem(t_samples(1:min(50, length(t_samples)))*1000, ...
     signal_sampled(sample_indices(1:min(50, length(sample_indices)))), ...
     'r', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('时间 t (ms)', 'FontSize', 12);
ylabel('幅度', 'FontSize', 12);
title('采样后的信号 x_s(t) = x(t) * s(t)', 'FontSize', 14);
axis([0 min(5, max(t_continuous)*1000) -2 2]);
legend('原始连续信号', '采样点', 'Location', 'best');
hold off;

% 子图4：连续信号的频谱（单边）
f_continuous_single = f_continuous(f_continuous >= 0);
X_continuous_single = X_continuous(f_continuous >= 0);
subplot(3,2,4);
plot(f_continuous_single/1000, abs(X_continuous_single), 'b-', 'LineWidth', 2);
grid on;
xlabel('频率 f (kHz)', 'FontSize', 12);
ylabel('|X(f)|', 'FontSize', 12);
title('连续信号的频谱（单边）', 'FontSize', 14);
axis([0 Fs/1000 0 max(abs(X_continuous_single))*1.1]);
xline(Fs/2000, 'r--', 'LineWidth', 2, 'DisplayName', sprintf('奈奎斯特频率 = %.1f kHz', Fs/2000));
legend('连续信号频谱', '奈奎斯特频率', 'Location', 'best');

% 子图5：采样后信号的频谱（周期延拓）
f_sampled_single = f_continuous(f_continuous >= 0);
X_sampled_single = X_sampled(f_continuous >= 0);
subplot(3,2,5);
plot(f_sampled_single/1000, abs(X_sampled_single), 'r-', 'LineWidth', 2);
grid on;
xlabel('频率 f (kHz)', 'FontSize', 12);
ylabel('|X_s(f)|', 'FontSize', 12);
title('采样后信号的频谱（周期延拓，每Fs重复）', 'FontSize', 14);
axis([0 Fs/1000 0 max(abs(X_sampled_single))*1.1]);
xline(Fs/2000, 'g--', 'LineWidth', 2, 'DisplayName', sprintf('奈奎斯特频率 = %.1f kHz', Fs/2000));
xline(Fs/1000, 'm--', 'LineWidth', 2, 'DisplayName', sprintf('采样率 = %.1f kHz', Fs/1000));
legend('采样信号频谱', '奈奎斯特频率', '采样率', 'Location', 'best');

% 子图6：采样函数的频谱（冲击函数序列的频谱）
S_impulse_single = S_impulse(f_continuous >= 0);
subplot(3,2,6);
plot(f_continuous_single/1000, abs(S_impulse_single), 'g-', 'LineWidth', 2);
grid on;
xlabel('频率 f (kHz)', 'FontSize', 12);
ylabel('|S(f)|', 'FontSize', 12);
title('采样函数 s(t) 的频谱（冲击函数序列）', 'FontSize', 14);
axis([0 Fs/1000 0 max(abs(S_impulse_single))*1.1]);
xline(Fs/1000, 'm--', 'LineWidth', 2, 'DisplayName', sprintf('采样率 = %.1f kHz', Fs/1000));
legend('采样函数频谱', '采样率', 'Location', 'best');

%===============================================================================
% 【可视化：带宽限制与混叠现象】
%===============================================================================
figure(4);
set(gcf, 'Position', [350, 350, 1400, 800]);

% 创建超过奈奎斯特频率的信号（会产生混叠）
f_alias = Fs / 2 + 2000;  % 超过奈奎斯特频率的信号
signal_alias = sin(2*pi*f_alias*t_continuous);
signal_alias_sampled = signal_alias .* impulse_train;

% 计算混叠信号的频谱
X_alias = fftshift(fft(signal_alias_sampled, N_continuous)) * (1/Fs_continuous);

% 计算混叠后的频率（根据采样定理）
f_alias_result = f_alias - Fs;  % 混叠到负频率，然后映射到正频率

subplot(2,2,1);
plot(t_continuous(1:min(1000, length(t_continuous)))*1000, ...
     signal_alias(1:min(1000, length(signal_alias))), 'b-', 'LineWidth', 1.5);
hold on;
stem(t_samples(1:min(20, length(t_samples)))*1000, ...
     signal_alias_sampled(sample_indices(1:min(20, length(sample_indices)))), ...
     'r', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('时间 t (ms)', 'FontSize', 12);
ylabel('幅度', 'FontSize', 12);
title(sprintf('高频信号（f = %.1f Hz > 奈奎斯特频率）的采样', f_alias), 'FontSize', 14);
axis([0 min(2, max(t_continuous))*1000 -1.5 1.5]);
legend('原始高频信号', '采样点', 'Location', 'best');
hold off;

subplot(2,2,2);
plot(f_continuous_single/1000, abs(X_alias(f_continuous >= 0)), 'r-', 'LineWidth', 2);
grid on;
xlabel('频率 f (kHz)', 'FontSize', 12);
ylabel('|X(f)|', 'FontSize', 12);
title('采样后信号的频谱（显示混叠）', 'FontSize', 14);
axis([0 Fs/1000 0 max(abs(X_alias(f_continuous >= 0)))*1.1]);
xline(Fs/2000, 'g--', 'LineWidth', 2, 'DisplayName', sprintf('奈奎斯特频率 = %.1f kHz', Fs/2000));
xline(f_alias/1000, 'b--', 'LineWidth', 1.5, 'DisplayName', sprintf('原始频率 = %.1f kHz', f_alias/1000));
xline(abs(f_alias_result)/1000, 'r--', 'LineWidth', 1.5, 'DisplayName', sprintf('混叠频率 = %.1f kHz', abs(f_alias_result)/1000));
legend('采样信号频谱', '奈奎斯特频率', '原始频率', '混叠频率', 'Location', 'best');

% 展示频谱周期延拓的完整视图
subplot(2,2,3);
f_extended = (-2*Fs:Fs/1000:2*Fs) / 1000;  % 扩展到多个周期
% 简化展示：只显示主要频率分量
plot([-Fs/1000, -Fs/2000, 0, Fs/2000, Fs/1000, 1.5*Fs/1000], ...
     [0, 0.5, 1, 0.5, 0, 0.5], 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('频率 f (kHz)', 'FontSize', 12);
ylabel('|X(f)| (归一化)', 'FontSize', 12);
title('采样后信号的频谱周期延拓（多个周期）', 'FontSize', 14);
axis([-2*Fs/1000 2*Fs/1000 -0.1 1.1]);
xline(-Fs/2000, 'g--', 'LineWidth', 1.5);
xline(Fs/2000, 'g--', 'LineWidth', 1.5);
xline(-Fs/1000, 'm--', 'LineWidth', 1.5);
xline(Fs/1000, 'm--', 'LineWidth', 1.5);
xline(0, 'k-', 'LineWidth', 0.5);
text(0, 1.05, 'DC', 'HorizontalAlignment', 'center', 'FontSize', 10);
text(Fs/2000, 0.6, 'Fs/2', 'HorizontalAlignment', 'center', 'FontSize', 10);
text(Fs/1000, 0.05, 'Fs', 'HorizontalAlignment', 'center', 'FontSize', 10);
legend('频谱周期延拓', 'Location', 'best');

% 解释图：为什么带宽限制在采样率的一半
subplot(2,2,4);
% 绘制理想低通滤波器响应
f_ideal = 0:Fs/1000/100:Fs/1000;
H_ideal = ones(size(f_ideal));
H_ideal(f_ideal > Fs/2000) = 0;
plot(f_ideal, H_ideal, 'b-', 'LineWidth', 3);
hold on;
fill([0, Fs/2000, Fs/2000, 0], [0, 0, 1, 1], 'b', 'FaceAlpha', 0.3);
grid on;
xlabel('频率 f (kHz)', 'FontSize', 12);
ylabel('幅度', 'FontSize', 12);
title('理想低通滤波器（带宽 = Fs/2）', 'FontSize', 14);
axis([0 Fs/1000 -0.1 1.2]);
xline(Fs/2000, 'r--', 'LineWidth', 2, 'DisplayName', sprintf('截止频率 = Fs/2 = %.1f kHz', Fs/2000));
text(Fs/4000, 0.5, '通带\n(0 到 Fs/2)', 'HorizontalAlignment', 'center', ...
    'FontSize', 12, 'BackgroundColor', 'white', 'EdgeColor', 'blue');
text(3*Fs/4000, 0.5, '阻带\n(> Fs/2)', 'HorizontalAlignment', 'center', ...
    'FontSize', 12, 'BackgroundColor', 'white', 'EdgeColor', 'red');
legend('理想低通滤波器', '截止频率', 'Location', 'best');
hold off;

%===============================================================================
% 【命令行输出：解释说明】
%===============================================================================
fprintf('==== 采样过程与时域表示的解释 =====\n');
fprintf('\n1. 采样过程在时域的表示：\n');
fprintf('   采样过程可以表示为：x_s(t) = x(t) * s(t)\n');
fprintf('   其中 s(t) = Σ δ(t - n*Ts) 是采样冲击函数序列\n');
fprintf('   Ts = 1/Fs = %.6f 秒 是采样间隔\n', Ts);
fprintf('   采样率 Fs = %.1f Hz 表示每秒采样 %.1f 次\n', Fs, Fs);
fprintf('   在时域，采样就是在离散时刻 t = n*Ts 对信号进行"冲击采样"\n');
fprintf('\n2. 采样函数的频谱特性：\n');
fprintf('   采样函数 s(t) = Σ δ(t - n*Ts) 的频谱是周期性的\n');
fprintf('   频谱在频率轴上每 Fs Hz 重复一次\n');
fprintf('   这导致采样后信号的频谱也是周期性的（周期 = Fs）\n');
fprintf('\n3. 为什么信号带宽限制在采样率的一半（奈奎斯特频率）：\n');
fprintf('   奈奎斯特频率 f_Nyquist = Fs/2 = %.1f Hz\n', Fs/2);
fprintf('   如果信号带宽 > Fs/2，频谱会发生混叠（Aliasing）\n');
fprintf('   混叠：高频分量会"折叠"到低频区域，无法区分\n');
fprintf('   因此，要无失真地恢复信号，信号带宽必须 ≤ Fs/2\n');
fprintf('\n4. 采样定理（奈奎斯特定理）：\n');
fprintf('   如果信号的最大频率为 f_max，则采样率必须满足：\n');
fprintf('   Fs ≥ 2 * f_max\n');
fprintf('   或者：f_max ≤ Fs/2（奈奎斯特频率）\n');
fprintf('   这样才能从采样信号中无失真地恢复原始连续信号\n');
fprintf('\n5. 实际例子（本仿真）：\n');
fprintf('   采样率 Fs = %.1f kHz\n', Fs/1000);
fprintf('   奈奎斯特频率 = %.1f kHz\n', Fs/2000);
fprintf('   信号带宽应限制在 ≤ %.1f kHz\n', Fs/2000);
fprintf('   如果信号频率 > %.1f kHz，会发生混叠\n', Fs/2000);
fprintf('   例如：%.1f Hz 的信号会混叠到 %.1f Hz\n', ...
    f_alias, abs(f_alias_result));
fprintf('==========================================\n\n');

toc

%===============================================================================
% 文件结束
%===============================================================================

