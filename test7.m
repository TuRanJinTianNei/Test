% ============================================================================
% 文件名: test7.m
% 功能: 模拟带通复包络信号并展示
% 描述:
%   1. 生成基带复包络信号（包含幅度和相位调制）
%   2. 将复包络调制到带通频率
%   3. 生成多个Figure展示：
%      - Figure 1: 复包络的I/Q分量（实部/虚部）
%      - Figure 2: 复包络的幅度和相位
%      - Figure 3: 带通信号的时域波形
%      - Figure 4: 复包络在复平面上的轨迹（星座图）
%      - Figure 5: 带通信号的频谱
%      - Figure 6: 复包络的频谱（基带）
% 原理: 
%   带通信号：s(t) = Re{ṡ(t) * exp(j*2π*fc*t)}
%   其中 ṡ(t) 是复包络：ṡ(t) = A(t)*exp(j*φ(t)) = I(t) + j*Q(t)
%   - A(t): 包络（幅度）
%   - φ(t): 相位
%   - I(t): 同相分量（实部）
%   - Q(t): 正交分量（虚部）
% ============================================================================

clc; clear; close all;

%===============================================================================
% 【系统参数配置】
%===============================================================================
fs = 10000;                      % 采样频率 (Hz)
T = 1;                           % 信号持续时间 (秒)
t = 0:1/fs:T-1/fs;              % 时间向量
N = length(t);                   % 采样点数

fc = 1000;                       % 载波频率 (Hz)，带通中心频率
f1 = 50;                         % 基带信号频率1 (Hz)
f2 = 100;                        % 基带信号频率2 (Hz)

%===============================================================================
% 【步骤1: 生成基带复包络信号】
%===============================================================================
% 方法1: 使用幅度和相位调制
% 幅度调制：A(t) = 1 + 0.5*cos(2π*f1*t)
A_t = 1 + 0.5*cos(2*pi*f1*t);

% 相位调制：φ(t) = 2π*f2*t + 0.3*sin(2π*f1*t)
phi_t = 2*pi*f2*t + 0.3*sin(2*pi*f1*t);

% 复包络：ṡ(t) = A(t)*exp(j*φ(t))
s_complex_envelope = A_t .* exp(1j*phi_t);

% 提取I/Q分量
I_t = real(s_complex_envelope);  % 同相分量（实部）
Q_t = imag(s_complex_envelope);  % 正交分量（虚部）

% 从I/Q分量计算幅度和相位（验证）
A_calc = abs(s_complex_envelope);
phi_calc = angle(s_complex_envelope);

%===============================================================================
% 【步骤2: 生成带通信号】
%===============================================================================
% 带通信号：s_bp(t) = Re{ṡ(t) * exp(j*2π*fc*t)}
s_bandpass = real(s_complex_envelope .* exp(1j*2*pi*fc*t));

% 也可以表示为：s_bp(t) = I(t)*cos(2π*fc*t) - Q(t)*sin(2π*fc*t)
s_bandpass_alt = I_t.*cos(2*pi*fc*t) - Q_t.*sin(2*pi*fc*t);

%===============================================================================
% 【可视化分析】
%===============================================================================

%------------------------------------------------------------------------------
% Figure 1: 复包络的I/Q分量（实部/虚部）
%------------------------------------------------------------------------------
figure(1);
subplot(2,1,1);
plot(t(1:1000), I_t(1:1000), 'b-', 'LineWidth', 1.5);
grid on;
xlabel('时间 t (秒)');
ylabel('I(t) - 同相分量（实部）');
title('复包络的I分量（前0.1秒）');
legend('I(t) = Re{ṡ(t)}', 'Location', 'best');

subplot(2,1,2);
plot(t(1:1000), Q_t(1:1000), 'r-', 'LineWidth', 1.5);
grid on;
xlabel('时间 t (秒)');
ylabel('Q(t) - 正交分量（虚部）');
title('复包络的Q分量（前0.1秒）');
legend('Q(t) = Im{ṡ(t)}', 'Location', 'best');

%------------------------------------------------------------------------------
% Figure 2: 复包络的幅度和相位
%------------------------------------------------------------------------------
figure(2);
subplot(2,1,1);
plot(t(1:1000), A_calc(1:1000), 'g-', 'LineWidth', 1.5);
hold on;
plot(t(1:1000), A_t(1:1000), 'k--', 'LineWidth', 1);
grid on;
xlabel('时间 t (秒)');
ylabel('幅度 A(t)');
title('复包络的幅度（前0.1秒）');
legend('计算值 |ṡ(t)|', '理论值 A(t)', 'Location', 'best');

subplot(2,1,2);
plot(t(1:1000), phi_calc(1:1000), 'm-', 'LineWidth', 1.5);
hold on;
plot(t(1:1000), phi_t(1:1000), 'k--', 'LineWidth', 1);
grid on;
xlabel('时间 t (秒)');
ylabel('相位 φ(t) (弧度)');
title('复包络的相位（前0.1秒）');
legend('计算值 arg{ṡ(t)}', '理论值 φ(t)', 'Location', 'best');

%------------------------------------------------------------------------------
% Figure 3: 带通信号的时域波形
%------------------------------------------------------------------------------
figure(3);
subplot(3,1,1);
plot(t(1:500), s_bandpass(1:500), 'b-', 'LineWidth', 1);
grid on;
xlabel('时间 t (秒)');
ylabel('幅度');
title(sprintf('带通信号 s(t) = Re{ṡ(t)*exp(j*2π*fc*t)} (fc=%d Hz, 前0.05秒)', fc));
legend('s_{bp}(t)', 'Location', 'best');

subplot(3,1,2);
plot(t(1:500), I_t(1:500), 'r-', 'LineWidth', 1.5);
hold on;
plot(t(1:500), Q_t(1:500), 'g-', 'LineWidth', 1.5);
grid on;
xlabel('时间 t (秒)');
ylabel('幅度');
title('复包络的I/Q分量（前0.05秒）');
legend('I(t)', 'Q(t)', 'Location', 'best');

subplot(3,1,3);
plot(t(1:500), A_calc(1:500), 'k-', 'LineWidth', 1.5);
grid on;
xlabel('时间 t (秒)');
ylabel('幅度');
title('复包络的幅度（前0.05秒）');
legend('|ṡ(t)|', 'Location', 'best');

%------------------------------------------------------------------------------
% Figure 4: 复包络在复平面上的轨迹（星座图）
%------------------------------------------------------------------------------
figure(4);
subplot(2,2,1);
plot(I_t(1:500), Q_t(1:500), 'b-', 'LineWidth', 1);
hold on;
plot(I_t(1), Q_t(1), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
plot(I_t(500), Q_t(500), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
grid on;
axis equal;
xlabel('I(t) - 实部');
ylabel('Q(t) - 虚部');
title('复包络在复平面上的轨迹（前0.05秒）');
legend('轨迹', '起点', '终点', 'Location', 'best');

subplot(2,2,2);
scatter(I_t(1:100:end), Q_t(1:100:end), 20, t(1:100:end), 'filled');
colorbar;
colormap('jet');
grid on;
axis equal;
xlabel('I(t) - 实部');
ylabel('Q(t) - 虚部');
title('复包络轨迹（颜色表示时间）');

subplot(2,2,3);
polarplot(phi_calc(1:500), A_calc(1:500), 'b-', 'LineWidth', 1);
title('复包络的极坐标表示（幅度-相位）');

subplot(2,2,4);
histogram2(I_t, Q_t, 30, 'DisplayStyle', 'tile', 'ShowEmptyBins', 'off');
xlabel('I(t)');
ylabel('Q(t)');
title('复包络的2D分布直方图');
colorbar;
axis equal;

%------------------------------------------------------------------------------
% Figure 5: 带通信号的频谱
%------------------------------------------------------------------------------
figure(5);
% 计算FFT
N_fft = 2^nextpow2(N);
S_bandpass_fft = fft(s_bandpass, N_fft);
f_fft = (0:N_fft-1) * fs / N_fft;

% 转换为单边频谱
S_bandpass_onesided = S_bandpass_fft(1:N_fft/2+1);
S_bandpass_onesided(2:end-1) = 2*S_bandpass_onesided(2:end-1);
f_onesided = f_fft(1:N_fft/2+1);

subplot(2,1,1);
plot(f_onesided, abs(S_bandpass_onesided), 'b-', 'LineWidth', 1.5);
grid on;
xlabel('频率 f (Hz)');
ylabel('幅度谱 |S(f)|');
title('带通信号的幅度谱');
xlim([0, 2*fc]);
hold on;
plot([fc, fc], ylim, 'r--', 'LineWidth', 1.5);
legend('幅度谱', sprintf('载波频率 fc=%d Hz', fc), 'Location', 'best');

subplot(2,1,2);
plot(f_onesided, 20*log10(abs(S_bandpass_onesided) + eps), 'b-', 'LineWidth', 1.5);
grid on;
xlabel('频率 f (Hz)');
ylabel('幅度谱 (dB)');
title('带通信号的幅度谱（对数尺度）');
xlim([0, 2*fc]);
hold on;
plot([fc, fc], ylim, 'r--', 'LineWidth', 1.5);
legend('幅度谱 (dB)', sprintf('载波频率 fc=%d Hz', fc), 'Location', 'best');

%------------------------------------------------------------------------------
% Figure 6: 复包络的频谱（基带）
%------------------------------------------------------------------------------
figure(6);
% 计算复包络的FFT
S_envelope_fft = fft(s_complex_envelope, N_fft);
S_envelope_onesided = S_envelope_fft(1:N_fft/2+1);
S_envelope_onesided(2:end-1) = 2*S_envelope_onesided(2:end-1);

subplot(2,1,1);
plot(f_onesided, abs(S_envelope_onesided), 'g-', 'LineWidth', 1.5);
grid on;
xlabel('频率 f (Hz)');
ylabel('幅度谱 |Ṡ(f)|');
title('复包络的幅度谱（基带）');
xlim([0, 300]);
hold on;
plot([f1, f1], ylim, 'r--', 'LineWidth', 1);
plot([f2, f2], ylim, 'm--', 'LineWidth', 1);
legend('复包络幅度谱', sprintf('f1=%d Hz', f1), sprintf('f2=%d Hz', f2), 'Location', 'best');

subplot(2,1,2);
plot(f_onesided, angle(S_envelope_onesided)*180/pi, 'g-', 'LineWidth', 1.5);
grid on;
xlabel('频率 f (Hz)');
ylabel('相位谱 (度)');
title('复包络的相位谱（基带）');
xlim([0, 300]);

%===============================================================================
% 【命令行输出摘要】
%===============================================================================
fprintf('\n==== 带通复包络信号模拟总结 ====\n');
fprintf('采样频率: %.0f Hz\n', fs);
fprintf('信号持续时间: %.1f 秒\n', T);
fprintf('采样点数: %d\n', N);
fprintf('载波频率: %.0f Hz\n', fc);
fprintf('基带频率1 (幅度调制): %.0f Hz\n', f1);
fprintf('基带频率2 (相位调制): %.0f Hz\n', f2);
fprintf('\n【信号表达式】\n');
fprintf('复包络: ṡ(t) = A(t)*exp(j*φ(t)) = I(t) + j*Q(t)\n');
fprintf('  其中: A(t) = 1 + 0.5*cos(2π*%d*t)\n', f1);
fprintf('       φ(t) = 2π*%d*t + 0.3*sin(2π*%d*t)\n', f2, f1);
fprintf('带通信号: s(t) = Re{ṡ(t)*exp(j*2π*%d*t)}\n', fc);
fprintf('        = I(t)*cos(2π*%d*t) - Q(t)*sin(2π*%d*t)\n', fc, fc);
fprintf('\n【统计特性】\n');
fprintf('复包络平均功率: %.6f\n', mean(abs(s_complex_envelope).^2));
fprintf('复包络平均幅度: %.6f\n', mean(abs(s_complex_envelope)));
fprintf('I分量均值: %.6f, 标准差: %.6f\n', mean(I_t), std(I_t));
fprintf('Q分量均值: %.6f, 标准差: %.6f\n', mean(Q_t), std(Q_t));
fprintf('带通信号平均功率: %.6f\n', mean(s_bandpass.^2));
fprintf('==================================\n\n');

