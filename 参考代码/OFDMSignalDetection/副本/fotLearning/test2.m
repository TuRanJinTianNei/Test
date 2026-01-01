%% 平稳信号 vs 循环平稳信号 - 循环谱方法对比仿真
% 本仿真对比展示：
% 1. 平稳信号（Stationary Signal）的循环谱特性
% 2. 循环平稳信号（Cyclostationary Signal）的循环谱特性
% 3. 解释为什么循环谱方法只适用于循环平稳信号

clear; close all; clc;

%% 1. 参数设置
fs = 40e6;          % 采样频率 (Hz)
fc = 10e6;          % 载波频率 (Hz)
N = 4096;           % 信号长度
M = 16;             % 平滑长度（用于循环谱计算）

% 时间轴
t = (0:N-1) / fs;
dt = 1/fs;

% 循环频率轴 α
d_alpha = fs / N;    % 循环频率分辨率
alpha_max = fc;      % 最大循环频率
alpha = -alpha_max:d_alpha:alpha_max;
alpha_len = length(alpha);

% 延迟轴 τ
tau_max_samples = min(round(N/fs/4*fs), N/4);
tau_samples = -tau_max_samples:tau_max_samples;
tau_len = length(tau_samples);
tau = tau_samples / fs;

% 频率轴 f（用于FFT）
N_fft = 2^nextpow2(tau_len);
f = (-N_fft/2:N_fft/2-1) * fs / N_fft;
f_len = length(f);

fprintf('=== 平稳信号 vs 循环平稳信号对比仿真 ===\n\n');

%% 2. 生成平稳信号（Stationary Signal）
% 平稳信号：统计特性不随时间变化
% 示例1：白噪声（完全平稳）
fprintf('生成平稳信号...\n');
rng(42);  % 设置随机种子以便复现
x_stationary1 = randn(1, N);  % 高斯白噪声
x_stationary1 = x_stationary1 / std(x_stationary1);

% 示例2：正弦波（确定性信号，可视为平稳）
f_sin = 5e6;  % 正弦波频率
x_stationary2 = cos(2*pi*f_sin*t);
x_stationary2 = x_stationary2 / std(x_stationary2);

% 示例3：带噪声的正弦波
x_stationary3 = cos(2*pi*f_sin*t) + 0.1*randn(1, N);
x_stationary3 = x_stationary3 / std(x_stationary3);

%% 3. 生成循环平稳信号（Cyclostationary Signal）
% 循环平稳信号：统计特性随时间呈周期性变化
fprintf('生成循环平稳信号...\n');

% 示例1：OFDM信号（典型的循环平稳信号）
if exist('ofdm.m', 'file')
    N_symbols = 20;
    para = 128;
    ratio = 1/4;
    sig_ofdm = ofdm(N_symbols, para, ratio);
    
    % 上变频到载波频率
    trst_rate = 5e5;
    s_n = ceil(fs/trst_rate);
    sign = sig_ofdm(ones(s_n,1),:);
    sign = reshape(sign, 1, s_n*length(sign));
    
    if length(sign) > N
        sign = sign(1:N);
    else
        sign = [sign, zeros(1, N-length(sign))];
    end
    
    x_cyclo1 = sign .* exp(1i*2*pi*fc/fs*(0:N-1));
    x_cyclo1 = real(x_cyclo1);
    x_cyclo1 = x_cyclo1 / std(x_cyclo1);
else
    % 如果没有ofdm函数，使用周期调制的正弦波
    T_sym = 1e-6;  % 符号周期
    modulation_freq = 1/T_sym;
    x_cyclo1 = cos(2*pi*fc*t + 0.5*cos(2*pi*modulation_freq*t));
    x_cyclo1 = x_cyclo1 / std(x_cyclo1);
end

% 示例2：周期调制的信号
T_mod = 1e-6;  % 调制周期
modulation_freq = 1/T_mod;
x_cyclo2 = cos(2*pi*fc*t) .* (1 + 0.5*cos(2*pi*modulation_freq*t));
x_cyclo2 = x_cyclo2 / std(x_cyclo2);

%% 4. 计算循环自相关函数和循环谱的辅助函数
% 定义计算函数（将在文件末尾定义）

%% 5. 计算各信号的循环谱
fprintf('\n计算循环谱...\n');

% 平稳信号1：白噪声
fprintf('计算平稳信号1（白噪声）的循环谱...\n');
[R_stat1, S_stat1] = compute_cyclic_spectrum(x_stationary1, N, fs, alpha, tau_samples, N_fft, f);

% 平稳信号2：正弦波
fprintf('计算平稳信号2（正弦波）的循环谱...\n');
[R_stat2, S_stat2] = compute_cyclic_spectrum(x_stationary2, N, fs, alpha, tau_samples, N_fft, f);

% 循环平稳信号1：OFDM
fprintf('计算循环平稳信号1（OFDM）的循环谱...\n');
[R_cyclo1, S_cyclo1] = compute_cyclic_spectrum(x_cyclo1, N, fs, alpha, tau_samples, N_fft, f);

% 循环平稳信号2：周期调制信号
fprintf('计算循环平稳信号2（周期调制）的循环谱...\n');
[R_cyclo2, S_cyclo2] = compute_cyclic_spectrum(x_cyclo2, N, fs, alpha, tau_samples, N_fft, f);

%% 6. 可视化对比结果

%% 6.1 信号时域波形对比
figure('Position', [50, 50, 1400, 900], 'Name', '平稳信号 vs 循环平稳信号对比');

% 平稳信号时域波形
subplot(3, 4, 1);
plot(t(1:min(500, N))*1e6, x_stationary1(1:min(500, N)));
xlabel('时间 t (μs)');
ylabel('幅度');
title('平稳信号1: 白噪声');
grid on;

subplot(3, 4, 2);
plot(t(1:min(500, N))*1e6, x_stationary2(1:min(500, N)));
xlabel('时间 t (μs)');
ylabel('幅度');
title('平稳信号2: 正弦波');
grid on;

% 循环平稳信号时域波形
subplot(3, 4, 3);
plot(t(1:min(500, N))*1e6, x_cyclo1(1:min(500, N)));
xlabel('时间 t (μs)');
ylabel('幅度');
title('循环平稳信号1: OFDM');
grid on;

subplot(3, 4, 4);
plot(t(1:min(500, N))*1e6, x_cyclo2(1:min(500, N)));
xlabel('时间 t (μs)');
ylabel('幅度');
title('循环平稳信号2: 周期调制');
grid on;

%% 6.2 循环谱在α=0处的切片（传统功率谱）
subplot(3, 4, 5);
alpha_zero_idx = find(abs(alpha) < d_alpha, 1);
if isempty(alpha_zero_idx)
    alpha_zero_idx = floor(alpha_len/2) + 1;
end
S_alpha0_stat1 = abs(S_stat1(alpha_zero_idx, :));
plot(f/1e6, 20*log10(S_alpha0_stat1 + eps));
xlabel('频率 f (MHz)');
ylabel('功率谱 (dB)');
title('平稳信号1: 功率谱 (α=0)');
grid on;
xlim([-fs/2/1e6, fs/2/1e6]);

subplot(3, 4, 6);
S_alpha0_stat2 = abs(S_stat2(alpha_zero_idx, :));
plot(f/1e6, 20*log10(S_alpha0_stat2 + eps));
xlabel('频率 f (MHz)');
ylabel('功率谱 (dB)');
title('平稳信号2: 功率谱 (α=0)');
grid on;
xlim([-fs/2/1e6, fs/2/1e6]);

subplot(3, 4, 7);
S_alpha0_cyclo1 = abs(S_cyclo1(alpha_zero_idx, :));
plot(f/1e6, 20*log10(S_alpha0_cyclo1 + eps));
xlabel('频率 f (MHz)');
ylabel('功率谱 (dB)');
title('循环平稳信号1: 功率谱 (α=0)');
grid on;
xlim([-fs/2/1e6, fs/2/1e6]);

subplot(3, 4, 8);
S_alpha0_cyclo2 = abs(S_cyclo2(alpha_zero_idx, :));
plot(f/1e6, 20*log10(S_alpha0_cyclo2 + eps));
xlabel('频率 f (MHz)');
ylabel('功率谱 (dB)');
title('循环平稳信号2: 功率谱 (α=0)');
grid on;
xlim([-fs/2/1e6, fs/2/1e6]);

%% 6.3 循环谱在特定频率f处的切片（关键对比）
% 对于平稳信号，除了α=0外，其他α处应该接近零
% 对于循环平稳信号，在α=fc处应该有峰值

f_plot_idx = find(abs(f - fc) < fs/N*10, 1);
if isempty(f_plot_idx)
    f_plot_idx = find(abs(f) < fs/4, 1);
end

subplot(3, 4, 9);
S_f_stat1 = abs(S_stat1(:, f_plot_idx));
plot(alpha/1e6, 20*log10(S_f_stat1 + eps));
xlabel('循环频率 α (MHz)');
ylabel('循环谱 (dB)');
title(sprintf('平稳信号1: S_α(f), f=%.2f MHz', f(f_plot_idx)/1e6));
grid on;
xlim([-alpha_max/1e6, alpha_max/1e6]);
hold on;
plot([0, 0], ylim, 'r--', 'LineWidth', 1);
plot([fc/1e6, fc/1e6], ylim, 'g--', 'LineWidth', 1);
hold off;
legend('循环谱', 'α=0', 'α=fc', 'Location', 'best');

subplot(3, 4, 10);
S_f_stat2 = abs(S_stat2(:, f_plot_idx));
plot(alpha/1e6, 20*log10(S_f_stat2 + eps));
xlabel('循环频率 α (MHz)');
ylabel('循环谱 (dB)');
title(sprintf('平稳信号2: S_α(f), f=%.2f MHz', f(f_plot_idx)/1e6));
grid on;
xlim([-alpha_max/1e6, alpha_max/1e6]);
hold on;
plot([0, 0], ylim, 'r--', 'LineWidth', 1);
plot([fc/1e6, fc/1e6], ylim, 'g--', 'LineWidth', 1);
hold off;
legend('循环谱', 'α=0', 'α=fc', 'Location', 'best');

subplot(3, 4, 11);
S_f_cyclo1 = abs(S_cyclo1(:, f_plot_idx));
plot(alpha/1e6, 20*log10(S_f_cyclo1 + eps));
xlabel('循环频率 α (MHz)');
ylabel('循环谱 (dB)');
title(sprintf('循环平稳信号1: S_α(f), f=%.2f MHz', f(f_plot_idx)/1e6));
grid on;
xlim([-alpha_max/1e6, alpha_max/1e6]);
hold on;
plot([0, 0], ylim, 'r--', 'LineWidth', 1);
plot([fc/1e6, fc/1e6], ylim, 'g--', 'LineWidth', 1);
hold off;
legend('循环谱', 'α=0', 'α=fc', 'Location', 'best');

subplot(3, 4, 12);
S_f_cyclo2 = abs(S_cyclo2(:, f_plot_idx));
plot(alpha/1e6, 20*log10(S_f_cyclo2 + eps));
xlabel('循环频率 α (MHz)');
ylabel('循环谱 (dB)');
title(sprintf('循环平稳信号2: S_α(f), f=%.2f MHz', f(f_plot_idx)/1e6));
grid on;
xlim([-alpha_max/1e6, alpha_max/1e6]);
hold on;
plot([0, 0], ylim, 'r--', 'LineWidth', 1);
plot([fc/1e6, fc/1e6], ylim, 'g--', 'LineWidth', 1);
hold off;
legend('循环谱', 'α=0', 'α=fc', 'Location', 'best');

%% 7. 循环谱三维图对比
figure('Position', [100, 100, 1400, 700], 'Name', '循环谱三维对比');

% 平稳信号1
subplot(2, 2, 1);
alpha_display_idx = 1:10:alpha_len;
f_display_idx = 1:5:f_len;
[Alpha_mesh, F_mesh] = meshgrid(alpha(alpha_display_idx)/1e6, f(f_display_idx)/1e6);
S_display = abs(S_stat1(alpha_display_idx, f_display_idx))';
S_display = 20*log10(S_display + eps);
mesh(Alpha_mesh, F_mesh, S_display);
xlabel('循环频率 α (MHz)');
ylabel('频率 f (MHz)');
zlabel('循环谱 (dB)');
title('平稳信号1: 循环谱三维图');
colorbar;
view(45, 30);

% 平稳信号2
subplot(2, 2, 2);
S_display = abs(S_stat2(alpha_display_idx, f_display_idx))';
S_display = 20*log10(S_display + eps);
mesh(Alpha_mesh, F_mesh, S_display);
xlabel('循环频率 α (MHz)');
ylabel('频率 f (MHz)');
zlabel('循环谱 (dB)');
title('平稳信号2: 循环谱三维图');
colorbar;
view(45, 30);

% 循环平稳信号1
subplot(2, 2, 3);
S_display = abs(S_cyclo1(alpha_display_idx, f_display_idx))';
S_display = 20*log10(S_display + eps);
mesh(Alpha_mesh, F_mesh, S_display);
xlabel('循环频率 α (MHz)');
ylabel('频率 f (MHz)');
zlabel('循环谱 (dB)');
title('循环平稳信号1: 循环谱三维图');
colorbar;
view(45, 30);

% 循环平稳信号2
subplot(2, 2, 4);
S_display = abs(S_cyclo2(alpha_display_idx, f_display_idx))';
S_display = 20*log10(S_display + eps);
mesh(Alpha_mesh, F_mesh, S_display);
xlabel('循环频率 α (MHz)');
ylabel('频率 f (MHz)');
zlabel('循环谱 (dB)');
title('循环平稳信号2: 循环谱三维图');
colorbar;
view(45, 30);

%% 8. 定量分析对比
fprintf('\n=== 定量分析对比 ===\n');

% 计算α=0和α=fc处的能量
alpha_zero_idx = find(abs(alpha) < d_alpha, 1);
if isempty(alpha_zero_idx)
    alpha_zero_idx = floor(alpha_len/2) + 1;
end

alpha_fc_idx = find(abs(alpha - fc) < d_alpha*5);
if isempty(alpha_fc_idx)
    alpha_fc_idx = find(abs(alpha) < alpha_max/4, 1);
end

fprintf('\n【平稳信号1：白噪声】\n');
energy_alpha0_stat1 = sum(abs(S_stat1(alpha_zero_idx, :)).^2);
if ~isempty(alpha_fc_idx)
    energy_alpha_fc_stat1 = sum(abs(S_stat1(alpha_fc_idx, :)).^2);
    ratio_stat1 = energy_alpha_fc_stat1 / energy_alpha0_stat1;
    fprintf('  α=0处能量: %.2e\n', energy_alpha0_stat1);
    fprintf('  α≈fc处能量: %.2e\n', energy_alpha_fc_stat1);
    fprintf('  能量比 (α=fc/α=0): %.6f\n', ratio_stat1);
    if ratio_stat1 < 0.01
        fprintf('  ✓ 符合平稳信号特征：α≠0处能量很小\n');
    else
        fprintf('  ⚠ 不符合平稳信号特征\n');
    end
end

fprintf('\n【平稳信号2：正弦波】\n');
energy_alpha0_stat2 = sum(abs(S_stat2(alpha_zero_idx, :)).^2);
if ~isempty(alpha_fc_idx)
    energy_alpha_fc_stat2 = sum(abs(S_stat2(alpha_fc_idx, :)).^2);
    ratio_stat2 = energy_alpha_fc_stat2 / energy_alpha0_stat2;
    fprintf('  α=0处能量: %.2e\n', energy_alpha0_stat2);
    fprintf('  α≈fc处能量: %.2e\n', energy_alpha_fc_stat2);
    fprintf('  能量比 (α=fc/α=0): %.6f\n', ratio_stat2);
    if ratio_stat2 < 0.01
        fprintf('  ✓ 符合平稳信号特征：α≠0处能量很小\n');
    else
        fprintf('  ⚠ 不符合平稳信号特征\n');
    end
end

fprintf('\n【循环平稳信号1：OFDM】\n');
energy_alpha0_cyclo1 = sum(abs(S_cyclo1(alpha_zero_idx, :)).^2);
if ~isempty(alpha_fc_idx)
    energy_alpha_fc_cyclo1 = sum(abs(S_cyclo1(alpha_fc_idx, :)).^2);
    ratio_cyclo1 = energy_alpha_fc_cyclo1 / energy_alpha0_cyclo1;
    fprintf('  α=0处能量: %.2e\n', energy_alpha0_cyclo1);
    fprintf('  α≈fc处能量: %.2e\n', energy_alpha_fc_cyclo1);
    fprintf('  能量比 (α=fc/α=0): %.6f\n', ratio_cyclo1);
    if ratio_cyclo1 > 0.01
        fprintf('  ✓ 符合循环平稳信号特征：α=fc处有明显能量\n');
    else
        fprintf('  ⚠ 循环平稳特性较弱\n');
    end
end

fprintf('\n【循环平稳信号2：周期调制】\n');
energy_alpha0_cyclo2 = sum(abs(S_cyclo2(alpha_zero_idx, :)).^2);
if ~isempty(alpha_fc_idx)
    energy_alpha_fc_cyclo2 = sum(abs(S_cyclo2(alpha_fc_idx, :)).^2);
    ratio_cyclo2 = energy_alpha_fc_cyclo2 / energy_alpha0_cyclo2;
    fprintf('  α=0处能量: %.2e\n', energy_alpha0_cyclo2);
    fprintf('  α≈fc处能量: %.2e\n', energy_alpha_fc_cyclo2);
    fprintf('  能量比 (α=fc/α=0): %.6f\n', ratio_cyclo2);
    if ratio_cyclo2 > 0.01
        fprintf('  ✓ 符合循环平稳信号特征：α=fc处有明显能量\n');
    else
        fprintf('  ⚠ 循环平稳特性较弱\n');
    end
end

%% 9. 总结说明
fprintf('\n=== 理论总结 ===\n');
fprintf('\n【平稳信号的特征】\n');
fprintf('1. 统计特性不随时间变化：E[X(t)] = μ (常数)\n');
fprintf('2. 自相关函数只与时间差有关：R_x(t,τ) = R_x(τ)\n');
fprintf('3. 循环自相关函数Rα(τ)只在α=0处有非零值\n');
fprintf('4. 循环谱Sα(f)在α≠0处接近零\n');
fprintf('5. 无法利用循环频率特征进行载波频率估计\n');

fprintf('\n【循环平稳信号的特征】\n');
fprintf('1. 统计特性随时间呈周期性变化：E[X(t)] = E[X(t+T)]\n');
fprintf('2. 自相关函数是周期函数：R_x(t,τ) = R_x(t+T,τ)\n');
fprintf('3. 循环自相关函数Rα(τ)在特定α处有非零值\n');
fprintf('4. 循环谱Sα(f)在α=0和α=k/T处有峰值（k为整数）\n');
fprintf('5. 可以利用循环频率特征进行载波频率估计\n');

fprintf('\n【循环谱方法的应用】\n');
fprintf('✓ 适用于循环平稳信号（如OFDM、QAM等调制信号）\n');
fprintf('✗ 不适用于平稳信号（如白噪声、纯正弦波等）\n');
fprintf('✗ 对于平稳信号，循环谱方法无法提取载波频率信息\n');

fprintf('\n仿真完成！\n');

%% 辅助函数：计算循环自相关函数和循环谱
function [R_alpha_tau, S_alpha_f] = compute_cyclic_spectrum(x, N, fs, alpha, tau_samples, N_fft, f)
    % 计算循环自相关函数 Rα(τ)
    alpha_len = length(alpha);
    tau_len = length(tau_samples);
    R_alpha_tau = zeros(alpha_len, tau_len);
    
    t_vec = (0:N-1) / fs;
    dt = 1/fs;
    
    % 预计算指数项
    exp_alpha_t = zeros(alpha_len, N);
    for i = 1:alpha_len
        exp_alpha_t(i, :) = exp(-1i*2*pi*alpha(i)*t_vec);
    end
    
    % 向量化计算
    for j = 1:tau_len
        tau_j_samples = tau_samples(j);
        
        idx1_vec = (1:N) + round(tau_j_samples/2);
        idx2_vec = (1:N) - round(tau_j_samples/2);
        
        valid_idx = (idx1_vec >= 1 & idx1_vec <= N & idx2_vec >= 1 & idx2_vec <= N);
        
        if sum(valid_idx) > 0
            x1 = zeros(1, N);
            x2 = zeros(1, N);
            x1(valid_idx) = x(idx1_vec(valid_idx));
            x2(valid_idx) = x(idx2_vec(valid_idx));
            
            x_product = x1 .* conj(x2);
            
            for i = 1:alpha_len
                integrand = x_product .* exp_alpha_t(i, :);
                R_alpha_tau(i, j) = mean(integrand(valid_idx));
            end
        end
        
        if mod(j, 50) == 0
            fprintf('  进度: %d/%d\n', j, tau_len);
        end
    end
    
    % 计算循环谱密度 Sα(f)
    f_len = length(f);
    S_alpha_f = zeros(alpha_len, f_len);
    
    for i = 1:alpha_len
        R_tau = R_alpha_tau(i, :);
        
        pad_left = floor((N_fft - tau_len) / 2);
        pad_right = N_fft - tau_len - pad_left;
        R_tau_padded = [zeros(1, pad_left), R_tau, zeros(1, pad_right)];
        
        S_f = fftshift(fft(R_tau_padded, N_fft));
        S_alpha_f(i, :) = S_f * dt;
    end
end

