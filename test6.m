% ============================================================================
% 文件名: test6.m
% 功能: Jakes信道数学特征与参数详解仿真
% 描述: 
%   1. 实现Jakes信道模型，展示其核心数学特征
%   2. 验证Jakes信道的统计特性：
%      - 幅度分布：瑞利分布
%      - 相位分布：均匀分布
%      - 功率谱：U型多普勒功率谱
%      - 自相关函数：J_0(2πf_d Δt)
%   3. 生成多个Figure详细展示Jakes信道特征
% ============================================================================

tic
clc;
clear all;
close all;

%===================== Figure 索引说明 =====================
%
% 【一】Jakes信道时域特征
%   Figure 1 : Jakes信道冲激响应（时域：幅度和相位）
%   Figure 2 : Jakes信道包络（幅度）随时间变化
%   Figure 3 : Jakes信道相位随时间变化
%   Figure 4 : Jakes信道实部和虚部随时间变化
%
% 【二】Jakes信道频域特征
%   Figure 5 : Jakes信道频率响应（幅度）
%   Figure 6 : Jakes信道频率响应（相位）
%   Figure 7 : Jakes信道多普勒功率谱密度（U型谱）
%
% 【三】Jakes信道统计特性验证
%   Figure 8 : Jakes信道包络分布（瑞利分布对比）
%   Figure 9 : Jakes信道相位分布（均匀分布对比）
%   Figure 10: Jakes信道自相关函数（J_0函数对比）
%
% 【四】Jakes信道参数影响分析
%   Figure 11: 不同fd下的多普勒功率谱对比
%   Figure 12: 不同N0下的包络分布对比
%
%==========================================================

%===============================================================================
% 【Jakes信道参数配置】
%===============================================================================
% 说明：以下参数定义了Jakes信道的核心配置

% 基本参数
fd = 926;                          % 最大多普勒频移（Hz）
% 示例：若载波频率fc=2.4 GHz，移动速度v=10 m/s，则：
%       fd = v*fc/c = 10*2.4e9/3e8 = 80 Hz
% 这里使用926 Hz作为示例（对应较高移动速度）

Ts = 1e-6;                         % 采样周期（s）
fs = 1/Ts;                         % 采样频率（Hz）
Ns = 500000;                       % 采样点数（信号长度）- 增加到500000以提高统计精度

% Jakes模型参数
N0 = 100;                          % 基础正弦波数量 - 增加到100以获得更好的相位均匀性
N = 4*N0 + 2;                      % 总正弦波数量（Jakes标准公式）= 402
% 注意：N0 >= 50才能获得较好的统计特性，N0 >= 100可以获得更好的相位均匀性
% 更大的N0（如200）可以获得更好的精度，但计算时间会增加

% 时间轴
t = (0:Ns-1) * Ts;                % 时间轴（秒）

fprintf('\n==== Jakes信道参数配置 ====\n');
fprintf('最大多普勒频移 fd     : %.2f Hz\n', fd);
fprintf('采样周期 Ts           : %.2e s\n', Ts);
fprintf('采样频率 fs           : %.2e Hz\n', fs);
fprintf('采样点数 Ns           : %d\n', Ns);
fprintf('基础正弦波数量 N0     : %d\n', N0);
fprintf('总正弦波数量 N        : %d\n', N);
fprintf('信号时长              : %.4f s\n', Ns*Ts);
fprintf('随机数生成器          : Mersenne Twister (高质量)\n');
fprintf('=====================================\n\n');

%===============================================================================
% 【Jakes信道生成】
%===============================================================================

%------------------------------------------------------------------------------
% Jakes信道数学模型
%------------------------------------------------------------------------------
% Jakes信道的复包络可表示为多个正弦波的叠加：
% h(t) = sum_{n=1}^{N} a_n * exp(j*(ω_n*t + φ_n))
%
% 其中：
% - a_n：第n条路径的幅度，服从瑞利分布
% - ω_n：第n条路径的多普勒频移，ω_n = 2π*fd*cos(α_n)
% - φ_n：第n条路径的随机相位
% - α_n：到达角度，在[0, 2π]均匀分布
%
% Jakes模型实现：
% h(t) = sqrt(2/N0) * [sum_{n=1}^{N0} cos(2π*fd*cos(β_n)*t + φ_n) + 
%                      j*sum_{n=1}^{N0} sin(2π*fd*cos(β_n)*t + φ_n)]
%
% 其中 β_n = π*n/N0，φ_n 是随机相位
%------------------------------------------------------------------------------

fprintf('\n==== 生成Jakes信道 ====\n');

% 初始化随机数生成器 - 使用Mersenne Twister（更高质量的随机数生成器）
rng('default');                    % 重置为默认状态
rng('shuffle', 'twister');         % 使用Mersenne Twister生成器，基于时间设置随机种子
% Mersenne Twister是MATLAB中质量最高的随机数生成器，周期为2^19937-1

% Jakes模型中的角度
beta_n = pi * (1:N0) / N0;         % β_n = π*n/N0，n=1,2,...,N0

% 随机初始相位（均匀分布在[0, 2π]）
% 为I和Q分量分别生成独立的随机相位，确保相位分布均匀
% 使用rand()函数，它基于当前设置的随机数生成器（Mersenne Twister）
phi_I_n = 2*pi*rand(1, N0);        % I分量的随机相位（高质量均匀分布）
phi_Q_n = 2*pi*rand(1, N0);        % Q分量的随机相位（独立于I分量，高质量均匀分布）

% 初始化同相分量（实部）和正交分量（虚部）
h_I = zeros(1, Ns);                % 同相分量（In-phase）
h_Q = zeros(1, Ns);                % 正交分量（Quadrature）

% 生成Jakes信道（改进的实现，确保I和Q分量独立）
for n = 1:N0
    % 多普勒频移：f_n = fd * cos(β_n)
    f_n = fd * cos(beta_n(n));
    
    % 同相分量（使用独立的相位）
    h_I = h_I + cos(2*pi*f_n*t + phi_I_n(n));
    
    % 正交分量（使用独立的相位）
    h_Q = h_Q + sin(2*pi*f_n*t + phi_Q_n(n));
end

% 归一化：sqrt(2/N0)
h_I = sqrt(2/N0) * h_I;
h_Q = sqrt(2/N0) * h_Q;

% 复信道冲激响应
h_jakes = h_I + 1j*h_Q;

% 归一化信道功率（使平均功率为1）
h_jakes = h_jakes / sqrt(mean(abs(h_jakes).^2));

fprintf('信道功率（归一化后）: %.6f\n', mean(abs(h_jakes).^2));
fprintf('信道包络均值        : %.6f\n', mean(abs(h_jakes)));
fprintf('信道包络标准差      : %.6f\n', std(abs(h_jakes)));
fprintf('=====================================\n\n');

%===============================================================================
% 【一】Jakes信道时域特征可视化】
%===============================================================================

%------------------------------------------------------------------------------
% Figure 1: Jakes信道冲激响应（时域：幅度和相位）
%------------------------------------------------------------------------------
figure(1);
subplot(2,1,1);
plot(t*1000, abs(h_jakes), 'b-', 'LineWidth', 1.5);
grid on;
xlabel('时间 (ms)');
ylabel('幅度 |h(t)|');
title(sprintf('Jakes信道冲激响应幅度（fd=%.0f Hz, N=%d）', fd, N));
legend('|h(t)|', 'Location', 'best');

subplot(2,1,2);
plot(t*1000, angle(h_jakes)*180/pi, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('时间 (ms)');
ylabel('相位 (度)');
title(sprintf('Jakes信道冲激响应相位（fd=%.0f Hz, N=%d）', fd, N));
legend('∠h(t)', 'Location', 'best');

%------------------------------------------------------------------------------
% Figure 2: Jakes信道包络（幅度）随时间变化
%------------------------------------------------------------------------------
figure(2);
plot(t*1000, 20*log10(abs(h_jakes) + eps), 'b-', 'LineWidth', 1.5);
grid on;
xlabel('时间 (ms)');
ylabel('幅度 (dB)');
title(sprintf('Jakes信道包络随时间变化（fd=%.0f Hz, N=%d）', fd, N));
legend('20log_{10}(|h(t)|)', 'Location', 'best');

%------------------------------------------------------------------------------
% Figure 3: Jakes信道相位随时间变化
%------------------------------------------------------------------------------
figure(3);
plot(t*1000, unwrap(angle(h_jakes))*180/pi, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('时间 (ms)');
ylabel('相位 (度)');
title(sprintf('Jakes信道相位随时间变化（fd=%.0f Hz, N=%d）', fd, N));
legend('∠h(t) (unwrapped)', 'Location', 'best');

%------------------------------------------------------------------------------
% Figure 4: Jakes信道实部和虚部随时间变化
%------------------------------------------------------------------------------
figure(4);
subplot(2,1,1);
plot(t*1000, real(h_jakes), 'b-', 'LineWidth', 1.5);
hold on;
plot(t*1000, imag(h_jakes), 'r-', 'LineWidth', 1.5);
grid on;
xlabel('时间 (ms)');
ylabel('幅度');
title(sprintf('Jakes信道实部和虚部（fd=%.0f Hz, N=%d）', fd, N));
legend('实部 Re[h(t)]', '虚部 Im[h(t)]', 'Location', 'best');
hold off;

subplot(2,1,2);
plot(real(h_jakes), imag(h_jakes), 'b.', 'MarkerSize', 1);
grid on;
xlabel('实部 Re[h(t)]');
ylabel('虚部 Im[h(t)]');
title('Jakes信道复平面轨迹');
axis equal;

%===============================================================================
% 【二】Jakes信道频域特征可视化
%===============================================================================

%------------------------------------------------------------------------------
% Figure 5: Jakes信道频率响应（幅度）
%------------------------------------------------------------------------------
Nfft = 2^18;                      % FFT点数 - 使用2^18=262144点以提高频率分辨率
% 如果Ns < Nfft，会进行零填充；如果Ns > Nfft，会截断
H_jakes = fft(h_jakes, Nfft);      % FFT变换
f_axis = (-Nfft/2:(Nfft/2-1)) * fs / Nfft;  % 频率轴（Hz）

figure(5);
subplot(2,1,1);
plot(f_axis/1000, 20*log10(abs(fftshift(H_jakes)) + eps), 'b-', 'LineWidth', 1.5);
grid on;
xlabel('频率 (kHz)');
ylabel('幅度 (dB)');
title('Jakes信道频率响应（幅度）');
legend('|H(f)|', 'Location', 'best');
xlim([-fd*2/1000, fd*2/1000]);

subplot(2,1,2);
plot(f_axis/1000, abs(fftshift(H_jakes)), 'b-', 'LineWidth', 1.5);
grid on;
xlabel('频率 (kHz)');
ylabel('幅度');
title('Jakes信道频率响应（幅度，线性刻度）');
legend('|H(f)|', 'Location', 'best');
xlim([-fd*2/1000, fd*2/1000]);

%------------------------------------------------------------------------------
% Figure 6: Jakes信道频率响应（相位）
%------------------------------------------------------------------------------
figure(6);
plot(f_axis/1000, angle(fftshift(H_jakes))*180/pi, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('频率 (kHz)');
ylabel('相位 (度)');
title('Jakes信道频率响应（相位）');
legend('∠H(f)', 'Location', 'best');
xlim([-fd*2/1000, fd*2/1000]);

%------------------------------------------------------------------------------
% Figure 7: Jakes信道多普勒功率谱密度（U型谱）
%------------------------------------------------------------------------------
% 计算多普勒功率谱密度（通过FFT的功率谱）
P_doppler = abs(fftshift(H_jakes)).^2;
f_doppler_normalized = f_axis / fd;  % 归一化到最大多普勒频移

% 理论Jakes多普勒功率谱密度
% S(f) = 1/(π*fd*sqrt(1-(f/fd)^2)) for |f| < fd
f_theory = linspace(-fd*1.1, fd*1.1, 1000);
S_theory = zeros(size(f_theory));
idx_valid = abs(f_theory) < fd;
S_theory(idx_valid) = 1./(pi*fd*sqrt(1 - (f_theory(idx_valid)/fd).^2));
S_theory = S_theory / max(S_theory);  % 归一化

figure(7);
subplot(2,1,1);
plot(f_doppler_normalized, 10*log10(P_doppler/max(P_doppler) + eps), 'b-', 'LineWidth', 1.5);
grid on;
xlabel('归一化频率 f/f_d');
ylabel('功率谱密度 (dB)');
title(sprintf('Jakes信道多普勒功率谱密度（fd=%.0f Hz）', fd));
legend('仿真 S(f)', 'Location', 'best');
axis([-2 2 -40 5]);

subplot(2,1,2);
plot(f_theory/fd, 10*log10(S_theory + eps), 'r--', 'LineWidth', 2);
hold on;
plot(f_doppler_normalized, 10*log10(P_doppler/max(P_doppler) + eps), 'b-', 'LineWidth', 1.5);
grid on;
xlabel('归一化频率 f/f_d');
ylabel('功率谱密度 (dB)');
title('Jakes信道多普勒功率谱密度（理论 vs 仿真）');
legend('理论', '仿真', 'Location', 'best');
axis([-1.5 1.5 -40 5]);
hold off;

%===============================================================================
% 【三】Jakes信道统计特性验证
%===============================================================================

%------------------------------------------------------------------------------
% Figure 8: Jakes信道包络分布（瑞利分布对比）
%------------------------------------------------------------------------------
h_envelope = abs(h_jakes);
[counts, centers] = hist(h_envelope, 100);  % 增加bin数量到100以提高精度
pdf_hist = counts / (sum(counts) * (centers(2) - centers(1)));

% 理论瑞利分布：p(r) = (r/σ²) * exp(-r²/(2σ²))
% 其中 σ² = E[|h|²]/2
sigma_squared = mean(h_envelope.^2) / 2;
sigma_rayleigh = sqrt(sigma_squared);
r_theory = linspace(0, max(h_envelope)*1.5, 1000);
pdf_rayleigh = (r_theory / sigma_squared) .* exp(-r_theory.^2 / (2*sigma_squared));

figure(8);
bar(centers, pdf_hist, 'FaceColor', [0.7 0.7 0.9], 'EdgeColor', 'none');
hold on;
plot(r_theory, pdf_rayleigh, 'r-', 'LineWidth', 2);
grid on;
xlabel('包络幅度 |h(t)|');
ylabel('概率密度');
title('Jakes信道包络分布（瑞利分布对比）');
legend('仿真直方图', sprintf('理论瑞利分布 (σ=%.4f)', sigma_rayleigh), 'Location', 'best');
hold off;

%------------------------------------------------------------------------------
% Figure 9: Jakes信道相位分布（均匀分布对比）
%------------------------------------------------------------------------------
h_phase = angle(h_jakes);
[counts_phase, centers_phase] = hist(h_phase, 100);  % 增加bin数量到100以提高精度
pdf_hist_phase = counts_phase / (sum(counts_phase) * (centers_phase(2) - centers_phase(1)));

% 理论均匀分布：p(θ) = 1/(2π) for -π ≤ θ ≤ π
theta_theory = linspace(-pi, pi, 1000);
pdf_uniform = ones(size(theta_theory)) / (2*pi);

figure(9);
bar(centers_phase, pdf_hist_phase, 'FaceColor', [0.9 0.7 0.7], 'EdgeColor', 'none');
hold on;
plot(theta_theory, pdf_uniform, 'r-', 'LineWidth', 2);
grid on;
xlabel('相位 (弧度)');
ylabel('概率密度');
title('Jakes信道相位分布（均匀分布对比）');
legend('仿真直方图', '理论均匀分布', 'Location', 'best');
xlim([-pi pi]);
hold off;

%------------------------------------------------------------------------------
% Figure 10: Jakes信道自相关函数（J_0函数对比）
%------------------------------------------------------------------------------
% 计算自相关函数
max_lag = min(20000, floor(Ns/4));  % 增加max_lag到20000以提高自相关函数精度
R_h = xcorr(h_jakes, max_lag, 'normalized');
lags = (-max_lag:max_lag) * Ts * 1000;  % 转换为毫秒

% 理论自相关函数：R(τ) = J_0(2πfd*τ)，其中J_0是零阶贝塞尔函数
tau_theory = linspace(0, max(lags), 1000);
R_theory = besselj(0, 2*pi*fd*tau_theory/1000);  % 注意单位转换（ms转s）

figure(10);
subplot(2,1,1);
plot(lags, abs(R_h), 'b-', 'LineWidth', 1.5);
hold on;
plot(tau_theory, abs(R_theory), 'r--', 'LineWidth', 2);
grid on;
xlabel('时延 τ (ms)');
ylabel('|R(τ)|');
title(sprintf('Jakes信道自相关函数（fd=%.0f Hz）', fd));
legend('仿真', '理论 J_0(2πf_dτ)', 'Location', 'best');
hold off;

subplot(2,1,2);
plot(lags, angle(R_h)*180/pi, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('时延 τ (ms)');
ylabel('相位 (度)');
title('Jakes信道自相关函数（相位）');
legend('∠R(τ)', 'Location', 'best');

%===============================================================================
% 【四】Jakes信道参数影响分析
%===============================================================================

%------------------------------------------------------------------------------
% Figure 11: 不同fd下的多普勒功率谱对比
%------------------------------------------------------------------------------
fd_values = [100, 500, 926];  % 不同的最大多普勒频移
colors = {'b', 'g', 'r'};

figure(11);
hold on;
for idx = 1:length(fd_values)
    fd_test = fd_values(idx);
    
    % 重新生成Jakes信道（使用相同的N0和随机种子）
    rng(42 + idx, 'twister');  % 不同fd使用不同随机种子，但可重复，使用Mersenne Twister
    beta_n_test = pi * (1:N0) / N0;
    % 为I和Q分量分别生成独立的随机相位
    phi_I_n_test = 2*pi*rand(1, N0);  % 高质量均匀分布
    phi_Q_n_test = 2*pi*rand(1, N0);  % 高质量均匀分布
    
    h_I_test = zeros(1, Ns);
    h_Q_test = zeros(1, Ns);
    
    for n = 1:N0
        f_n_test = fd_test * cos(beta_n_test(n));
        h_I_test = h_I_test + cos(2*pi*f_n_test*t + phi_I_n_test(n));
        h_Q_test = h_Q_test + sin(2*pi*f_n_test*t + phi_Q_n_test(n));
    end
    
    h_I_test = sqrt(2/N0) * h_I_test;
    h_Q_test = sqrt(2/N0) * h_Q_test;
    h_test = (h_I_test + 1j*h_Q_test) / sqrt(mean(abs(h_I_test + 1j*h_Q_test).^2));
    
    % 计算功率谱（使用更大的FFT以提高频率分辨率）
    Nfft_test = 2^18;  % 使用2^18=262144点FFT进一步提高频率分辨率
    H_test = fft(h_test, Nfft_test);
    P_test = abs(fftshift(H_test)).^2;
    
    % 计算频率轴（基于fd_test）
    f_axis_test = (-Nfft_test/2:(Nfft_test/2-1)) * fs / Nfft_test;
    f_test_normalized = f_axis_test / fd_test;
    
    % 只保留归一化频率在[-1.5, 1.5]范围内的点（U型谱的有效范围）
    idx_valid = abs(f_test_normalized) <= 1.5;
    f_plot = f_test_normalized(idx_valid);
    P_plot = P_test(idx_valid);
    
    % 归一化功率谱
    P_plot_normalized = P_plot / max(P_plot);
    
    % 平滑处理（使用移动平均减少毛刺，提高U型谱的可见性）
    window_size = max(1, floor(length(P_plot)/200));  % 自适应窗口大小
    if window_size > 1
        % 使用conv实现移动平均（兼容旧版MATLAB）
        window = ones(1, window_size) / window_size;
        P_plot_smooth = conv(P_plot_normalized, window, 'same');
    else
        P_plot_smooth = P_plot_normalized;
    end
    
    plot(f_plot, 10*log10(P_plot_smooth + eps), ...
         [colors{idx} '-'], 'LineWidth', 1.5, 'DisplayName', sprintf('fd=%.0f Hz', fd_test));
end

% 添加理论U型谱作为参考（使用fd=926作为示例）
f_theory_ref = linspace(-1.5, 1.5, 1000);
S_theory_ref = zeros(size(f_theory_ref));
idx_valid_ref = abs(f_theory_ref) < 1;
S_theory_ref(idx_valid_ref) = 1./(pi*sqrt(1 - f_theory_ref(idx_valid_ref).^2));
S_theory_ref = S_theory_ref / max(S_theory_ref);
plot(f_theory_ref, 10*log10(S_theory_ref + eps), 'k--', 'LineWidth', 1.5, ...
     'DisplayName', '理论U型谱');

grid on;
xlabel('归一化频率 f/f_d');
ylabel('功率谱密度 (dB)');
title('不同fd下的Jakes信道多普勒功率谱对比');
legend('Location', 'best');
axis([-1.5 1.5 -40 5]);
hold off;

%------------------------------------------------------------------------------
% Figure 12: 不同N0下的包络分布对比
%------------------------------------------------------------------------------
N0_values = [10, 34, 50, 100];  % 不同的正弦波数量 - 增加100以展示更高精度
colors = {'b', 'g', 'r', 'm'};   % 增加颜色以匹配4个N0值

figure(12);
hold on;
for idx = 1:length(N0_values)
    N0_test = N0_values(idx);
    
    % 重新生成Jakes信道（使用相同的fd和随机种子）
    rng(123 + idx, 'twister');  % 不同N0使用不同随机种子，但可重复，使用Mersenne Twister
    beta_n_test = pi * (1:N0_test) / N0_test;
    % 为I和Q分量分别生成独立的随机相位
    phi_I_n_test = 2*pi*rand(1, N0_test);  % 高质量均匀分布
    phi_Q_n_test = 2*pi*rand(1, N0_test);  % 高质量均匀分布
    
    h_I_test = zeros(1, Ns);
    h_Q_test = zeros(1, Ns);
    
    for n = 1:N0_test
        f_n_test = fd * cos(beta_n_test(n));
        h_I_test = h_I_test + cos(2*pi*f_n_test*t + phi_I_n_test(n));
        h_Q_test = h_Q_test + sin(2*pi*f_n_test*t + phi_Q_n_test(n));
    end
    
    h_I_test = sqrt(2/N0_test) * h_I_test;
    h_Q_test = sqrt(2/N0_test) * h_Q_test;
    h_test = (h_I_test + 1j*h_Q_test) / sqrt(mean(abs(h_I_test + 1j*h_Q_test).^2));
    
    % 计算包络分布
    h_envelope_test = abs(h_test);
    [counts_test, centers_test] = hist(h_envelope_test, 80);  % 增加bin数量到80以提高精度
    pdf_test = counts_test / (sum(counts_test) * (centers_test(2) - centers_test(1)));
    
    plot(centers_test, pdf_test, [colors{idx} '-'], 'LineWidth', 1.5, ...
         'DisplayName', sprintf('N0=%d', N0_test));
end

% 理论瑞利分布
sigma_squared_theory = 0.5;  % 归一化后E[|h|²]=1，所以σ²=0.5
r_theory = linspace(0, 3, 1000);
pdf_rayleigh_theory = (r_theory / sigma_squared_theory) .* exp(-r_theory.^2 / (2*sigma_squared_theory));
plot(r_theory, pdf_rayleigh_theory, 'k--', 'LineWidth', 2, 'DisplayName', '理论瑞利分布');

grid on;
xlabel('包络幅度 |h(t)|');
ylabel('概率密度');
title('不同N0下的Jakes信道包络分布对比');
legend('Location', 'best');
hold off;

%===============================================================================
% 【统计特性验证与输出】
%===============================================================================

fprintf('\n==== Jakes信道统计特性验证 ====\n');

% 包络统计特性
h_envelope = abs(h_jakes);
fprintf('包络统计特性：\n');
fprintf('  均值 E[|h|]        : %.6f\n', mean(h_envelope));
fprintf('  理论值 (σ√(π/2))   : %.6f\n', sigma_rayleigh*sqrt(pi/2));
fprintf('  标准差 std(|h|)    : %.6f\n', std(h_envelope));
fprintf('  理论值 (σ√(2-π/2)) : %.6f\n', sigma_rayleigh*sqrt(2-pi/2));
fprintf('  方差 var(|h|)      : %.6f\n', var(h_envelope));
fprintf('  理论值 (σ²(2-π/2)) : %.6f\n', sigma_squared*(2-pi/2));

% 相位统计特性
h_phase = angle(h_jakes);
fprintf('\n相位统计特性：\n');
fprintf('  均值 E[∠h]         : %.6f rad (%.2f 度)\n', mean(h_phase), mean(h_phase)*180/pi);
fprintf('  理论值（均匀分布）  : 0 rad (0 度)\n');
fprintf('  标准差 std(∠h)     : %.6f rad (%.2f 度)\n', std(h_phase), std(h_phase)*180/pi);
fprintf('  理论值（均匀分布）  : %.6f rad (%.2f 度)\n', pi/sqrt(3), pi/sqrt(3)*180/pi);

% 功率谱特性
fprintf('\n功率谱特性：\n');
fprintf('  最大多普勒频移 fd  : %.2f Hz\n', fd);
fprintf('  功率谱范围         : [-fd, fd] = [%.2f, %.2f] Hz\n', -fd, fd);
fprintf('  功率谱形状         : U型谱（在f=0处功率最大）\n');

% 自相关函数特性
fprintf('\n自相关函数特性：\n');
fprintf('  理论形式           : R(τ) = J_0(2πf_d*τ)\n');
fprintf('  其中 J_0 为零阶贝塞尔函数\n');
fprintf('  相关时间（近似）   : τ_c ≈ 1/(2π*fd) = %.6f s\n', 1/(2*pi*fd));

fprintf('=====================================\n\n');

toc
%===============================================================================
% 文件结束
%===============================================================================

