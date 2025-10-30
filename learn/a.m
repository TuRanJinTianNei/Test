% 清除环境变量，设置字体（确保中文显示）
clear; clc; close all;
set(groot,'DefaultAxesFontName','SimHei');
set(groot,'DefaultTextFontName','SimHei');

%% 1. 参数设置（核心参数与信号生成）
fs = 200;          % 采样频率：200Hz
T_total = 1;       % 信号总时长：1秒
t = 0:1/fs:T_total-1/fs;  % 时间向量
M = length(t);     % 原始信号长度：M = fs*T_total = 200点

% 生成含两个相近频率的信号（便于观察分辨率差异）
f1 = 50;   % 第一个频率：50Hz
f2 = 52;   % 第二个频率：52Hz（与f1相差2Hz）
x = sin(2*pi*f1*t) + sin(2*pi*f2*t);  % 原始信号


%% 2. 定义不同DFT点数（覆盖截断、补零、2的幂次等场景）
N_small = 64;      % 小点数（<M，需截断原始信号）
N_large = 200;     % 中等点数（=M，无截断/补零）
N_pad = 1024;      % 大点数（>M，需补零至1024点，2的幂次）
N_nonpow2 = 999;   % 非2的幂次点数（用于对比FFT效率）


%% 3. 对不同点数进行信号处理（截断或补零）
% 截断：原始信号取前N_small点
x_trunc = x(1:N_small);

% 补零：原始信号补零至N_pad点
x_pad = [x, zeros(1, N_pad - M)];

% 非2的幂次补零：原始信号补零至N_nonpow2点
x_nonpow2 = [x, zeros(1, N_nonpow2 - M)];


%% 4. 计算各情况下的FFT（DFT的快速实现）
X_small = fft(x_trunc, N_small);    % 64点FFT（截断信号）
X_large = fft(x, N_large);          % 200点FFT（原始信号）
X_pad = fft(x_pad, N_pad);          % 1024点FFT（补零信号）


%% 5. 计算频率轴（转换为实际频率值）
f_small = (0:N_small-1)*fs/N_small;    % 64点频率轴
f_large = (0:N_large-1)*fs/N_large;    % 200点频率轴
f_pad = (0:N_pad-1)*fs/N_pad;          % 1024点频率轴


%% 6. 绘制频谱图（取单边谱，实信号频谱对称）
figure(1); set(gcf,'Name','不同点数DFT的频谱对比','NumberTitle','off');

% 子图1：64点FFT（截断信号，低分辨率）
subplot(3,1,1);
plot(f_small, 2*abs(X_small)/N_small);  % 幅度归一化（单边谱乘2）
xlim([45, 60]);  % 聚焦感兴趣的频率范围
xlabel('频率（Hz）'); ylabel('幅度');
title(['64点DFT（截断信号）：分辨率Δf = ', num2str(fs/N_small),'Hz']);
grid on;

% 子图2：200点FFT（原始信号，中等分辨率）
subplot(3,1,2);
plot(f_large, 2*abs(X_large)/N_large);
xlim([45, 60]);
xlabel('频率（Hz）'); ylabel('幅度');
title(['200点DFT（原始信号）：分辨率Δf = ', num2str(fs/N_large),'Hz']);
grid on;

% 子图3：1024点FFT（补零信号，高计算分辨率）
subplot(3,1,3);
plot(f_pad, 2*abs(X_pad)/N_pad);
xlim([45, 60]);
xlabel('频率（Hz）'); ylabel('幅度');
title(['1024点DFT（补零信号）：分辨率Δf = ', num2str(fs/N_pad),'Hz']);
grid on;

sgtitle('点数对频率分辨率的影响');


%% 7. 验证补零不增加实际信息（对比补零前后的频谱细节）
figure(2); set(gcf,'Name','补零对频谱的影响（放大细节）','NumberTitle','off');
% 原始信号200点FFT的频谱（取45-60Hz）
idx_large = find(f_large >=45 & f_large <=60);
plot(f_large(idx_large), 2*abs(X_large(idx_large))/N_large, 'b', 'LineWidth',1.5);
hold on;
% 补零到1024点的频谱（同一频率范围）
idx_pad = find(f_pad >=45 & f_pad <=60);
plot(f_pad(idx_pad), 2*abs(X_pad(idx_pad))/N_pad, 'r--', 'LineWidth',1);
xlabel('频率（Hz）'); ylabel('幅度');
legend('200点原始信号','1024点补零信号');
title('补零仅使频谱更平滑，不增加实际频率信息');
grid on;


%% 7A. 更多可视化：时域对比（截断/原始/补零）
figure(3); set(gcf,'Name','时域对比：截断/原始/补零','NumberTitle','off');
subplot(2,1,1);
plot(t(1:N_small), x(1:N_small), 'b', 'LineWidth',1.0); hold on;
plot(t(1:N_small), x_trunc, 'r--', 'LineWidth',1.0);
hold off; grid on;
xlabel('时间（s）'); ylabel('幅度');
title(['前 ', num2str(N_small), ' 个样本：原始 vs 截断']);
legend('原始','截断');
subplot(2,1,2);
plot(t, x, 'b', 'LineWidth',1.0); hold on;
plot(0:1/fs:(length(x_pad)-1)/fs, x_pad, 'k--', 'LineWidth',1.0);
hold off; grid on;
xlabel('时间（s）'); ylabel('幅度');
title(['全时长：原始 vs 补零（至 ', num2str(N_pad), ' 点）']);
legend('原始','补零');


%% 7B. 更多可视化：窗函数对泄漏的影响（以64点为例）
figure(4); set(gcf,'Name','窗函数对泄漏的影响（64点DFT）','NumberTitle','off');
w_hann = hann(N_small)';
X_rect = fft(x_trunc, N_small);
X_hann = fft(x_trunc .* w_hann, N_small);
plot(f_small, 2*abs(X_rect)/N_small, 'r-', 'LineWidth',1.0); hold on;
plot(f_small, 2*abs(X_hann)/N_small, 'g--', 'LineWidth',1.0);
hold off; grid on;
xlim([45, 60]);
xlabel('频率（Hz）'); ylabel('幅度');
title('矩形窗 vs Hann窗：旁瓣降低但主瓣加宽');
legend('矩形窗','Hann窗');


%% 7C. 更多可视化：短时傅里叶变换（STFT）时频图
figure(5); set(gcf,'Name','STFT时频图','NumberTitle','off');
spectrogram(x, hann(64), 32, 256, fs, 'yaxis');
title('STFT：双频稳态信号的时频分布');
ylim([0 120]); colorbar;


%% 7D. 更多可视化：相位谱（200点DFT）
figure(6); set(gcf,'Name','相位谱（200点DFT）','NumberTitle','off');
phi_large = unwrap(angle(X_large));
plot(f_large, phi_large, 'm'); grid on;
xlim([45, 60]);
xlabel('频率（Hz）'); ylabel('相位（rad）');
title('相位谱（聚焦45–60 Hz）');


%% 8. 对比FFT计算效率（2的幂次 vs 非2的幂次）
% 重复多次计算，取平均时间（减少误差）
n_rep = 1000;  % 重复次数

% 2的幂次点数（1024）的计算时间
tic;
for i=1:n_rep
    fft(x_pad, N_pad);
end
t_pow2 = toc/n_rep;  % 平均时间

% 非2的幂次点数（999）的计算时间
tic;
for i=1:n_rep
    fft(x_nonpow2, N_nonpow2);
end
t_nonpow2 = toc/n_rep;  % 平均时间

% 输出效率对比
fprintf('FFT计算效率对比：\n');
fprintf('1024点（2的幂次）平均时间：%.6f 秒\n', t_pow2);
fprintf('999点（非2的幂次）平均时间：%.6f 秒\n', t_nonpow2);
fprintf('2的幂次点数计算速度提升：%.2f 倍\n', t_nonpow2/t_pow2);