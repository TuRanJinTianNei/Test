%% 循环谱（Cyclostationary Spectrum）仿真
% 根据循环自相关函数计算循环谱密度
% 实现图片中描述的理论
%
% ========================================================================
% 【如何理解信号的循环平稳特性】
% ========================================================================
%
% 1. 什么是循环平稳信号？
%    --------------------------------------------
%    循环平稳信号是指其统计特性随时间呈周期性变化的信号。
%    
%    数学定义：
%    - 均值周期性：E[x(t)] = E[x(t+T)]，对于所有t和某个周期T
%    - 方差周期性：Var[x(t)] = Var[x(t+T)]
%    - 自相关周期性：R_x(t,τ) = R_x(t+T,τ)
%    
%    通俗理解：
%    - 信号的"统计规律"（平均值、相关性等）每隔一段时间T就重复一次
%    - 虽然信号本身不是周期性的，但它的统计特性是周期性的
%
% 2. 为什么OFDM信号是循环平稳的？
%    --------------------------------------------
%    OFDM信号具有以下周期性特征：
%    
%    (a) 符号结构周期性：
%        - 每个OFDM符号都有相同的结构：IFFT部分 + 循环前缀
%        - 符号周期 T_sym = IFFT点数 + 循环前缀长度
%        - 每个符号的统计特性相同（因为数据是随机但结构相同）
%    
%    (b) 循环前缀的周期性：
%        - 每个符号的循环前缀都是符号尾部的复制
%        - 这创造了明显的周期性结构
%        - 导致信号在时域上呈现周期性模式
%    
%    (c) 统计特性的周期性：
%        - 由于符号结构相同，信号的均值、方差、自相关函数
%          在每个符号周期内重复
%        - 即：R_x(t,τ) = R_x(t+T_sym,τ)
%
% 3. 如何观察循环平稳特性？
%    --------------------------------------------
%    方法1：观察时域波形
%    - 查看多个符号的波形，应该看到相似的结构
%    - 循环前缀部分应该与符号尾部相同
%    
%    方法2：分析统计特性
%    - 计算滑动窗口的均值和方差，应该呈现周期性
%    - 验证自相关函数：R_x(t,τ) ≈ R_x(t+T_sym,τ)
%    
%    方法3：循环谱分析
%    - 在循环频率α=k/T_sym处（k为整数）应该出现峰值
%    - 在α=fc（载波频率）处也应该出现峰值
%
% 4. 循环平稳特性与循环谱的关系
%    --------------------------------------------
%    循环自相关函数：Rα(τ) = E[x(t+τ/2)x*(t-τ/2)e^(-j2παt)]
%    
%    对于循环平稳信号：
%    - Rα(τ)只在特定的循环频率α处有非零值
%    - α = 0: 对应传统功率谱（所有信号都有）
%    - α = k/T_sym: 对应符号周期的谐波（循环平稳信号特有）
%    - α = fc: 对应载波频率（上变频后的信号）
%    
%    循环谱密度：Sα(f) = FFT{Rα(τ)}
%    - 在(α, f)平面上，循环平稳信号会在特定位置出现峰值
%    - 这些峰值位置揭示了信号的周期性调制特征
%
% 5. 为什么循环谱方法需要循环平稳信号？
%    --------------------------------------------
%    - 对于平稳信号：Rα(τ)在所有α≠0处都接近零
%    - 对于循环平稳信号：Rα(τ)在特定α处有非零值
%    - 只有循环平稳信号才能在循环频率域提取载波频率等信息
%
% ========================================================================

clear; close all; clc;

%% 1. 参数设置
fs = 40e6;          % 采样频率 (Hz)
fc = 10e6;          % 载波频率 (Hz)
T0 = 1/fc;          % 信号统计特性的周期 (s)
alpha_max = fc;     % 最大循环频率 (Hz)
N = 2048;           % 信号长度
M = 16;             % 平滑长度

% 时间轴
t = (0:N-1) / fs;
dt = 1/fs;

% 延迟轴 τ
tau_max = N/fs / 4;  % 最大延迟
tau_max_samples = min(round(tau_max*fs), N/4);  % 最大延迟（采样点数）
tau_samples = -tau_max_samples:tau_max_samples;
tau_len = length(tau_samples);
tau = tau_samples / fs;

% 循环频率轴 α
d_alpha = fs / N;    % 循环频率分辨率
alpha = -alpha_max:d_alpha:alpha_max;
alpha_len = length(alpha);

% 频率轴 f（用于FFT）
N_fft = 2^nextpow2(tau_len);  % FFT长度
f = (-N_fft/2:N_fft/2-1) * fs / N_fft;
f_len = length(f);

%% 2. 生成循环平稳信号（OFDM信号示例）
% 使用简单的周期调制信号作为示例
% x(t) = A*cos(2*pi*fc*t + phi(t))
% 其中 phi(t) 是周期性的相位调制

% 生成OFDM信号（使用项目中的函数）
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
    
    % 截取到指定长度
    if length(sign) > N
        sign = sign(1:N);
    else
        sign = [sign, zeros(1, N-length(sign))];
    end
    
    % 上变频
    x = sign .* exp(1i*2*pi*fc/fs*(0:N-1));
    x = real(x);  % 取实部
else
    % 如果没有ofdm函数，使用简单的周期调制信号
    A = 1;
    modulation_freq = fc / 10;  % 调制频率
    x = A * cos(2*pi*fc*t + 0.5*cos(2*pi*modulation_freq*t));
end

% 归一化
x = x / std(x);

%% 2.1 分析信号的循环平稳特性
% 理解循环平稳特性的关键：
% 1. 信号的统计特性（均值、方差、自相关）是周期性的
% 2. 对于OFDM信号，周期T = 符号周期T_sym
% 3. 这意味着：E[x(t)] = E[x(t+T)], R_x(t,τ) = R_x(t+T,τ)

fprintf('\n=== 循环平稳特性分析 ===\n');

% 计算符号周期（如果使用OFDM信号）
if exist('ofdm.m', 'file')
    % OFDM符号周期 = IFFT点数 + 循环前缀长度
    IFFT_length = para;  % 子载波数 = IFFT点数
    CP_length = round(IFFT_length * ratio);  % 循环前缀长度
    T_sym_samples = IFFT_length + CP_length;  % 符号周期（采样点数）
    T_sym = T_sym_samples / fs;  % 符号周期（秒）
    
    fprintf('OFDM信号参数：\n');
    fprintf('  - 子载波数（IFFT点数）: %d\n', IFFT_length);
    fprintf('  - 循环前缀长度: %d 采样点\n', CP_length);
    fprintf('  - 符号周期 T_sym: %.2f μs (%d 采样点)\n', T_sym*1e6, T_sym_samples);
    fprintf('  - 符号速率: %.2f kbps\n', 1/T_sym/1e3);
    
    % 可视化符号结构
    figure('Position', [50, 50, 1400, 600], 'Name', 'OFDM信号循环平稳特性分析');
    
    % 显示前几个符号的时域波形
    num_symbols_show = min(5, floor(N/T_sym_samples));
    subplot(2, 3, 1);
    for i = 1:num_symbols_show
        idx_start = (i-1)*T_sym_samples + 1;
        idx_end = min(i*T_sym_samples, N);
        if idx_end > idx_start
            t_symbol = (idx_start:idx_end) / fs * 1e6;
            plot(t_symbol, x(idx_start:idx_end), 'LineWidth', 1.5);
            hold on;
        end
    end
    xlabel('时间 t (μs)');
    ylabel('幅度');
    title(sprintf('前%d个OFDM符号的时域波形', num_symbols_show));
    grid on;
    legend(arrayfun(@(i) sprintf('符号%d', i), 1:num_symbols_show, 'UniformOutput', false), ...
        'Location', 'best');
    
    % 分析统计特性的周期性
    % 计算滑动窗口的均值和方差
    window_size = round(T_sym_samples / 4);  % 窗口大小为符号周期的1/4
    num_windows = floor(N / window_size);
    
    mean_values = zeros(1, num_windows);
    var_values = zeros(1, num_windows);
    window_centers = zeros(1, num_windows);
    
    for i = 1:num_windows
        idx_start = (i-1)*window_size + 1;
        idx_end = min(i*window_size, N);
        window_centers(i) = (idx_start + idx_end) / 2 / fs * 1e6;
        mean_values(i) = mean(x(idx_start:idx_end));
        var_values(i) = var(x(idx_start:idx_end));
    end
    
    % 绘制均值随时间的变化
    subplot(2, 3, 2);
    plot(window_centers, mean_values, 'b-', 'LineWidth', 1.5);
    hold on;
    % 标记符号周期边界
    for i = 1:floor(N/T_sym_samples)
        xline(i*T_sym*1e6, 'r--', 'LineWidth', 1);
    end
    hold off;
    xlabel('时间 t (μs)');
    ylabel('局部均值');
    title('信号均值随时间的变化（应呈现周期性）');
    grid on;
    legend('局部均值', '符号周期边界', 'Location', 'best');
    
    % 绘制方差随时间的变化
    subplot(2, 3, 3);
    plot(window_centers, var_values, 'g-', 'LineWidth', 1.5);
    hold on;
    % 标记符号周期边界
    for i = 1:floor(N/T_sym_samples)
        xline(i*T_sym*1e6, 'r--', 'LineWidth', 1);
    end
    hold off;
    xlabel('时间 t (μs)');
    ylabel('局部方差');
    title('信号方差随时间的变化（应呈现周期性）');
    grid on;
    legend('局部方差', '符号周期边界', 'Location', 'best');
    
    % 计算自相关函数在不同时间点的值
    % 验证 R_x(t,τ) = R_x(t+T,τ)
    tau_test = round(T_sym_samples / 8);  % 测试延迟
    num_test_points = floor(N / T_sym_samples) - 1;
    
    autocorr_t1 = zeros(1, num_test_points);
    autocorr_t2 = zeros(1, num_test_points);
    t_test = zeros(1, num_test_points);
    
    for i = 1:num_test_points
        t1 = (i-1)*T_sym_samples + 1;
        t2 = i*T_sym_samples + 1;  % t2 = t1 + T_sym
        
        if t1 + tau_test <= N && t2 + tau_test <= N
            % 计算 R_x(t1, tau_test) 和 R_x(t2, tau_test)
            autocorr_t1(i) = mean(x(t1:t1+tau_test-1) .* x(t1+tau_test:t1+2*tau_test-1));
            autocorr_t2(i) = mean(x(t2:t2+tau_test-1) .* x(t2+tau_test:t2+2*tau_test-1));
            t_test(i) = t1 / fs * 1e6;
        end
    end
    
    % 绘制自相关函数的周期性
    subplot(2, 3, 4);
    plot(t_test, autocorr_t1, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4);
    hold on;
    plot(t_test, autocorr_t2, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 4);
    hold off;
    xlabel('时间 t (μs)');
    ylabel('自相关值');
    title(sprintf('自相关函数周期性验证 R_x(t,τ) vs R_x(t+T,τ), τ=%.2f μs', tau_test/fs*1e6));
    grid on;
    legend('R_x(t,τ)', 'R_x(t+T,τ)', 'Location', 'best');
    
    % 计算相关性（应该接近1）
    if length(autocorr_t1) > 1 && length(autocorr_t2) > 1
        corr_coef = corrcoef(autocorr_t1, autocorr_t2);
        if size(corr_coef, 1) == 2
            fprintf('  自相关函数周期性验证：\n');
            fprintf('    R_x(t,τ) 与 R_x(t+T,τ) 的相关系数: %.4f\n', corr_coef(1,2));
            if corr_coef(1,2) > 0.7
                fprintf('    ✓ 信号具有明显的循环平稳特性\n');
            else
                fprintf('    ⚠ 循环平稳特性较弱\n');
            end
        end
    end
    
    % 绘制循环前缀的周期性
    subplot(2, 3, 5);
    % 显示第一个符号的循环前缀和主体部分
    symbol1_start = 1;
    symbol1_cp_end = CP_length;
    symbol1_end = T_sym_samples;
    
    if symbol1_end <= N
        t_symbol1 = (symbol1_start:symbol1_end) / fs * 1e6;
        plot(t_symbol1, x(symbol1_start:symbol1_end), 'b-', 'LineWidth', 1.5);
        hold on;
        % 标记循环前缀部分
        t_cp = (symbol1_start:symbol1_cp_end) / fs * 1e6;
        plot(t_cp, x(symbol1_start:symbol1_cp_end), 'r-', 'LineWidth', 2);
        % 标记主体部分（应该与循环前缀的尾部相同）
        t_body_tail = ((symbol1_end-CP_length+1):symbol1_end) / fs * 1e6;
        plot(t_body_tail, x((symbol1_end-CP_length+1):symbol1_end), 'g--', 'LineWidth', 2);
        hold off;
        xlabel('时间 t (μs)');
        ylabel('幅度');
        title('OFDM符号结构：循环前缀的周期性');
        grid on;
        legend('完整符号', '循环前缀', '符号尾部（应与CP相同）', 'Location', 'best');
    end
    
    % 功率谱密度（展示频域特性）
    subplot(2, 3, 6);
    X_fft = fftshift(fft(x, N));
    f_fft = (-N/2:N/2-1) * fs / N;
    plot(f_fft/1e6, 20*log10(abs(X_fft) + eps));
    xlabel('频率 f (MHz)');
    ylabel('功率谱 (dB)');
    title('信号功率谱密度');
    grid on;
    xlim([-fs/2/1e6, fs/2/1e6]);
    
    fprintf('\n循环平稳特性的关键特征：\n');
    fprintf('  1. 符号结构周期性：每个符号都有相同的结构（IFFT部分+循环前缀）\n');
    fprintf('  2. 统计特性周期性：均值、方差、自相关函数在符号周期上重复\n');
    fprintf('  3. 循环前缀周期性：每个符号的循环前缀都是符号尾部的复制\n');
    fprintf('  4. 这些周期性导致信号在循环频率α=k/T_sym处出现峰值\n');
    fprintf('    其中k为整数，T_sym为符号周期\n');
end

%% 3. 计算循环自相关函数 Rα(τ)
% Rα(τ) = lim(T→∞) [1/T * ∫ x(t+τ/2) x*(t-τ/2) e^(-j2παt) dt]
% 使用向量化方法提高计算效率

fprintf('正在计算循环自相关函数 Rα(τ)...\n');

R_alpha_tau = zeros(alpha_len, tau_len);

% 预计算指数项
exp_alpha_t = zeros(alpha_len, N);
for i = 1:alpha_len
    exp_alpha_t(i, :) = exp(-1i*2*pi*alpha(i)*t);
end

% 向量化计算
for j = 1:tau_len
    tau_j_samples = tau_samples(j);
    tau_j = tau(j);
    
    % 计算 x(t+τ/2) 和 x(t-τ/2)
    idx1_vec = (1:N) + round(tau_j_samples/2);
    idx2_vec = (1:N) - round(tau_j_samples/2);
    
    % 边界处理
    valid_idx = (idx1_vec >= 1 & idx1_vec <= N & idx2_vec >= 1 & idx2_vec <= N);
    
    if sum(valid_idx) > 0
        x1 = zeros(1, N);
        x2 = zeros(1, N);
        x1(valid_idx) = x(idx1_vec(valid_idx));
        x2(valid_idx) = x(idx2_vec(valid_idx));
        
        % 计算 x(t+τ/2) * x*(t-τ/2)
        x_product = x1 .* conj(x2);
        
        % 对每个循环频率计算积分
        for i = 1:alpha_len
            integrand = x_product .* exp_alpha_t(i, :);
            R_alpha_tau(i, j) = mean(integrand(valid_idx));
        end
    end
    
    if mod(j, 10) == 0
        fprintf('进度: %d/%d\n', j, tau_len);
    end
end

%% 4. 计算循环谱密度 Sα(f)
% Sα(f) = ∫ Rα(τ) e^(-j2πfτ) dτ
% 对每个循环频率α，对τ进行傅里叶变换

fprintf('正在计算循环谱密度 Sα(f)...\n');
S_alpha_f = zeros(alpha_len, f_len);

for i = 1:alpha_len
    % 对Rα(τ)进行FFT
    R_tau = R_alpha_tau(i, :);
    
    % 补零到N_fft长度（居中补零）
    pad_left = floor((N_fft - tau_len) / 2);
    pad_right = N_fft - tau_len - pad_left;
    R_tau_padded = [zeros(1, pad_left), R_tau, zeros(1, pad_right)];
    
    % FFT变换
    S_f = fftshift(fft(R_tau_padded, N_fft));
    
    % 归一化
    S_alpha_f(i, :) = S_f * dt;
end

%% 5. 可视化结果

% 5.1 原始信号时域波形
figure('Position', [100, 100, 1200, 800]);

subplot(2, 3, 1);
plot(t(1:min(1000, N))*1e6, x(1:min(1000, N)));
xlabel('时间 t (μs)');
ylabel('幅度');
title('原始信号时域波形');
grid on;

% 5.2 信号的功率谱密度（α=0的循环谱切片）
subplot(2, 3, 2);
alpha_zero_idx = find(abs(alpha) < d_alpha, 1);
if isempty(alpha_zero_idx)
    alpha_zero_idx = floor(alpha_len/2) + 1;
end
S_alpha0 = abs(S_alpha_f(alpha_zero_idx, :));
plot(f/1e6, 20*log10(S_alpha0 + eps));
xlabel('频率 f (MHz)');
ylabel('功率谱密度 (dB)');
title('功率谱密度 S_0(f) (α=0的切片)');
grid on;
xlim([-fs/2/1e6, fs/2/1e6]);

% 5.3 循环自相关函数 Rα(τ) 在特定α下的切片
subplot(2, 3, 3);
alpha_plot_idx = find(abs(alpha - fc) < d_alpha*5, 1);
if isempty(alpha_plot_idx)
    alpha_plot_idx = find(abs(alpha) < alpha_max/4, 1);
end
R_alpha_plot = abs(R_alpha_tau(alpha_plot_idx, :));
plot(tau*1e6, R_alpha_plot);
xlabel('延迟 τ (μs)');
ylabel('|R_α(τ)|');
title(sprintf('循环自相关函数 R_α(τ), α=%.2f MHz', alpha(alpha_plot_idx)/1e6));
grid on;

% 5.4 循环谱密度 Sα(f) 的三维图
subplot(2, 3, [4, 5, 6]);
% 选择部分数据进行显示以提高速度
alpha_display_idx = 1:5:alpha_len;
f_display_idx = 1:2:f_len;
[Alpha_mesh, F_mesh] = meshgrid(alpha(alpha_display_idx)/1e6, f(f_display_idx)/1e6);
S_display = abs(S_alpha_f(alpha_display_idx, f_display_idx))';
S_display = 20*log10(S_display + eps);

mesh(Alpha_mesh, F_mesh, S_display);
xlabel('循环频率 α (MHz)');
ylabel('频率 f (MHz)');
zlabel('循环谱密度 (dB)');
title('循环谱密度 S_α(f) 三维图');
colorbar;
view(45, 30);

% 5.5 循环谱密度在特定频率f=fc处的切片
figure('Position', [200, 200, 1000, 600]);

subplot(1, 2, 1);
f_plot_idx = find(abs(f - fc) < fs/N*10, 1);
if isempty(f_plot_idx)
    f_plot_idx = find(abs(f) < fs/4, 1);
end
S_f_plot = abs(S_alpha_f(:, f_plot_idx));
plot(alpha/1e6, 20*log10(S_f_plot + eps));
xlabel('循环频率 α (MHz)');
ylabel('循环谱密度 (dB)');
title(sprintf('循环谱密度 S_α(f), f=%.2f MHz', f(f_plot_idx)/1e6));
grid on;
xlim([-alpha_max/1e6, alpha_max/1e6]);

% 5.6 循环谱密度的二维等高线图
subplot(1, 2, 2);
alpha_contour_idx = 1:10:alpha_len;
f_contour_idx = 1:5:f_len;
[Alpha_contour, F_contour] = meshgrid(alpha(alpha_contour_idx)/1e6, f(f_contour_idx)/1e6);
S_contour = abs(S_alpha_f(alpha_contour_idx, f_contour_idx))';
S_contour = 20*log10(S_contour + eps);

contour(Alpha_contour, F_contour, S_contour, 20);
xlabel('循环频率 α (MHz)');
ylabel('频率 f (MHz)');
title('循环谱密度 S_α(f) 等高线图');
colorbar;
grid on;

%% 6. 理论验证与循环平稳性分析
fprintf('\n=== 理论验证 ===\n');
fprintf('载波频率 fc = %.2f MHz\n', fc/1e6);
fprintf('循环频率分辨率 d_alpha = %.2f kHz\n', d_alpha/1e3);
fprintf('频率分辨率 df = %.2f kHz\n', fs/N/1e3);

% 在α=fc处查找峰值
alpha_fc_idx = find(abs(alpha - fc) < d_alpha*5);
if ~isempty(alpha_fc_idx)
    [max_val, max_idx] = max(abs(S_alpha_f(alpha_fc_idx, :)));
    fprintf('在α≈fc处找到峰值，幅度 = %.2f dB\n', 20*log10(max_val + eps));
end

fprintf('\n=== 循环平稳性验证 ===\n');
fprintf('【循环平稳信号的特征】\n');
fprintf('1. 循环自相关函数Rα(τ)只在特定循环频率α处有非零值\n');
fprintf('2. 对于OFDM信号，主要峰值出现在：\n');
fprintf('   - α = 0: 传统功率谱（所有信号都有）\n');
fprintf('   - α = ±fc: 载波频率（循环平稳信号特有）\n');
fprintf('   - α = k/T_sym: 符号周期相关的循环频率\n');

% 分析循环频率α=0和α=fc处的能量
alpha_zero_idx = find(abs(alpha) < d_alpha, 1);
if isempty(alpha_zero_idx)
    alpha_zero_idx = floor(alpha_len/2) + 1;
end
energy_alpha0 = sum(abs(S_alpha_f(alpha_zero_idx, :)).^2);
fprintf('\n3. 能量分析：\n');
fprintf('   - α=0处的总能量: %.2e\n', energy_alpha0);

if ~isempty(alpha_fc_idx)
    energy_alpha_fc = sum(abs(S_alpha_f(alpha_fc_idx, :)).^2);
    fprintf('   - α≈fc处的总能量: %.2e\n', energy_alpha_fc);
    fprintf('   - 能量比 (α=fc/α=0): %.4f\n', energy_alpha_fc/energy_alpha0);
    
    if energy_alpha_fc/energy_alpha0 > 0.01
        fprintf('   ✓ 信号在α=fc处有明显的循环平稳特性\n');
    else
        fprintf('   ⚠ 信号在α=fc处的循环平稳特性较弱\n');
        fprintf('   ⚠ 可能原因：信号不是循环平稳的，或数据长度不足\n');
    end
end

fprintf('\n【重要提示】\n');
fprintf('循环谱方法要求信号必须是循环平稳信号。如果信号不是循环平稳的：\n');
fprintf('- 循环自相关函数Rα(τ)在所有α处都接近零（除了α=0）\n');
fprintf('- 无法利用循环频率α=fc处的峰值来估计载波频率\n');
fprintf('- 估计结果可能不准确或失败\n');
fprintf('\nOFDM信号是典型的循环平稳信号，因为：\n');
fprintf('- 具有周期性的符号结构\n');
fprintf('- 循环前缀的周期性重复\n');
fprintf('- 统计特性在符号周期上重复\n');

fprintf('\n仿真完成！\n');

