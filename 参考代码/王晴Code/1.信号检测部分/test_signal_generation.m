% 测试OFDM信号生成函数，并展示时域和频域特性
% 
% 本脚本演示如何使用 generate_OFDM_signal 函数生成OFDM信号，
% 并绘制信号的时域波形和频域功率谱

clc;
clear all;
close all;

%% 参数设置
SubCarryNN = 256;      % 有效子载波数
SubCarryN = 1024;      % 总子载波数
SymbN = 20;            % OFDM符号个数
ratio = 1/4;           % 循环前缀比例
SNR = 10;              % 信噪比（dB），设为空[]则不添加噪声

%% 生成OFDM信号
fprintf('正在生成OFDM信号...\n');
fprintf('参数设置：\n');
fprintf('  有效子载波数：%d\n', SubCarryNN);
fprintf('  总子载波数：%d\n', SubCarryN);
fprintf('  OFDM符号数：%d\n', SymbN);
fprintf('  循环前缀比例：%.2f\n', ratio);
if ~isempty(SNR)
    fprintf('  信噪比：%.1f dB\n', SNR);
else
    fprintf('  信噪比：无噪声\n');
end

[TrData, ReData] = generate_OFDM_signal(SubCarryNN, SubCarryN, SymbN, ratio, SNR);

fprintf('信号生成完成！\n');
fprintf('信号长度：%d 个采样点\n', length(TrData));
fprintf('信号持续时间：%.4f 秒（假设采样率为1024 Hz）\n', length(TrData)/1024);

%% 计算采样率
fs = 1024;  % 采样率（Hz）
Ts = 1/fs;  % 采样间隔
t = (0:length(ReData)-1) * Ts;  % 时间轴

%% ========== 时域分析 ==========
figure('Position', [100, 100, 1200, 800]);

% 子图1：时域波形（实部）
subplot(3, 2, 1);
plot(t(1:min(5000, length(t))), real(ReData(1:min(5000, length(ReData)))), 'b-', 'LineWidth', 1);
xlabel('时间 (s)');
ylabel('幅度');
title('时域波形 - 实部（前5000个采样点）');
grid on;
xlim([t(1), t(min(5000, length(t)))]);

% 子图2：时域波形（虚部）
subplot(3, 2, 2);
plot(t(1:min(5000, length(t))), imag(ReData(1:min(5000, length(ReData)))), 'r-', 'LineWidth', 1);
xlabel('时间 (s)');
ylabel('幅度');
title('时域波形 - 虚部（前5000个采样点）');
grid on;
xlim([t(1), t(min(5000, length(t)))]);

% 子图3：时域幅度
subplot(3, 2, 3);
plot(t(1:min(5000, length(t))), abs(ReData(1:min(5000, length(ReData)))), 'g-', 'LineWidth', 1);
xlabel('时间 (s)');
ylabel('幅度');
title('时域波形 - 幅度（前5000个采样点）');
grid on;
xlim([t(1), t(min(5000, length(t)))]);

% 子图4：时域相位
subplot(3, 2, 4);
plot(t(1:min(5000, length(t))), angle(ReData(1:min(5000, length(ReData)))), 'm-', 'LineWidth', 1);
xlabel('时间 (s)');
ylabel('相位 (rad)');
title('时域波形 - 相位（前5000个采样点）');
grid on;
xlim([t(1), t(min(5000, length(t)))]);

% 子图5：一个OFDM符号的时域波形（包含循环前缀）
symbol_length = SubCarryN + SubCarryN * ratio;  % 一个符号的长度
if length(ReData) >= symbol_length
    subplot(3, 2, 5);
    t_symbol = (0:symbol_length-1) * Ts;
    plot(t_symbol, real(ReData(1:symbol_length)), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(t_symbol, imag(ReData(1:symbol_length)), 'r--', 'LineWidth', 1.5);
    xlabel('时间 (s)');
    ylabel('幅度');
    title(sprintf('单个OFDM符号时域波形（长度=%d）', symbol_length));
    legend('实部', '虚部', 'Location', 'best');
    grid on;
    
    % 标记循环前缀
    cp_length = SubCarryN * ratio;
    xline(t_symbol(cp_length), 'k--', 'LineWidth', 1, 'Label', '循环前缀结束');
end

% 子图6：信号功率统计
subplot(3, 2, 6);
histogram(abs(ReData), 50, 'FaceColor', 'c', 'EdgeColor', 'none');
xlabel('幅度');
ylabel('频数');
title('信号幅度分布直方图');
grid on;

sgtitle('OFDM信号时域特性分析', 'FontSize', 14, 'FontWeight', 'bold');

%% ========== 频域分析 ==========
figure('Position', [150, 150, 1200, 800]);

% 计算功率谱密度（使用Welch方法）
window = hanning(150);
noverlap = 50;
[Pxx1, f_welch] = pwelch(ReData, window, noverlap, length(ReData), fs);
Pxx_welch = 10 * log10(fftshift(Pxx1));
f_welch_shifted = fftshift(f_welch);

% 子图1：Welch法功率谱密度
subplot(3, 2, 1);
plot(f_welch_shifted, Pxx_welch, 'b-', 'LineWidth', 1.5);
xlabel('频率 (Hz)');
ylabel('功率谱密度 (dB/Hz)');
title('Welch法功率谱密度估计');
grid on;

% 子图2：FFT频谱（线性刻度）
subplot(3, 2, 2);
N_fft = 2^nextpow2(length(ReData));  % FFT长度（2的幂次）
Y = fft(ReData, N_fft);
P = abs(Y).^2 / N_fft;  % 功率谱
f_fft = (0:N_fft/2-1) * fs / N_fft;
plot(f_fft, P(1:N_fft/2), 'r-', 'LineWidth', 1.5);
xlabel('频率 (Hz)');
ylabel('功率');
title('FFT功率谱（线性刻度）');
grid on;
xlim([0, fs/2]);

% 子图3：FFT频谱（对数刻度）
subplot(3, 2, 3);
plot(f_fft, 10*log10(P(1:N_fft/2)), 'g-', 'LineWidth', 1.5);
xlabel('频率 (Hz)');
ylabel('功率 (dB)');
title('FFT功率谱（对数刻度）');
grid on;
xlim([0, fs/2]);

% 子图4：频域幅度谱
subplot(3, 2, 4);
Y_shifted = fftshift(Y);
f_shifted = (-N_fft/2:N_fft/2-1) * fs / N_fft;
plot(f_shifted, abs(Y_shifted), 'm-', 'LineWidth', 1.5);
xlabel('频率 (Hz)');
ylabel('幅度');
title('频域幅度谱（FFT）');
grid on;
xlim([-fs/2, fs/2]);

% 子图5：频域相位谱
subplot(3, 2, 5);
plot(f_shifted, angle(Y_shifted), 'c-', 'LineWidth', 1.5);
xlabel('频率 (Hz)');
ylabel('相位 (rad)');
title('频域相位谱（FFT）');
grid on;
xlim([-fs/2, fs/2]);

% 子图6：子载波功率分布（仅显示有效子载波）
subplot(3, 2, 6);
% 提取有效子载波位置的功率
valid_start = N_fft/2 - SubCarryNN/2 + 1;
valid_end = N_fft/2 + SubCarryNN/2;
valid_indices = valid_start:valid_end;
subcarrier_power = P(valid_indices);
subcarrier_idx = 1:length(subcarrier_power);
bar(subcarrier_idx, 10*log10(subcarrier_power), 'FaceColor', 'y', 'EdgeColor', 'k');
xlabel('子载波索引');
ylabel('功率 (dB)');
title(sprintf('有效子载波功率分布（共%d个子载波）', SubCarryNN));
grid on;

sgtitle('OFDM信号频域特性分析', 'FontSize', 14, 'FontWeight', 'bold');

%% ========== 信号统计信息 ==========
fprintf('\n========== 信号统计信息 ==========\n');
fprintf('时域统计：\n');
fprintf('  均值（实部）：%.6f\n', mean(real(ReData)));
fprintf('  均值（虚部）：%.6f\n', mean(imag(ReData)));
fprintf('  标准差（实部）：%.6f\n', std(real(ReData)));
fprintf('  标准差（虚部）：%.6f\n', std(imag(ReData)));
fprintf('  平均功率：%.6f\n', mean(abs(ReData).^2));
fprintf('  峰值功率：%.6f\n', max(abs(ReData).^2));
fprintf('  峰均功率比（PAPR）：%.2f dB\n', 10*log10(max(abs(ReData).^2) / mean(abs(ReData).^2)));

fprintf('\n频域统计：\n');
fprintf('  总带宽：%.2f Hz\n', fs);
fprintf('  子载波间隔：%.2f Hz\n', fs / SubCarryN);
fprintf('  有效带宽：%.2f Hz\n', fs * SubCarryNN / SubCarryN);
fprintf('  频谱利用率：%.2f%%\n', SubCarryNN / SubCarryN * 100);

fprintf('\n信号参数：\n');
fprintf('  采样点数：%d\n', length(ReData));
fprintf('  采样率：%.2f Hz\n', fs);
fprintf('  信号时长：%.4f 秒\n', length(ReData) / fs);
fprintf('  每个OFDM符号长度：%d 采样点\n', symbol_length);
fprintf('  循环前缀长度：%d 采样点\n', SubCarryN * ratio);
fprintf('  有效符号长度：%d 采样点\n', SubCarryN);

fprintf('\n========== 分析完成 ==========\n');

