% OFDM 加窗影响的详细仿真（时域与频域对比）
% 目标：展示加窗后主瓣略微展宽、旁瓣显著抑制，信息不丢失

tic
clc;
clear;
close all;

%================= 参数设置 =================
carrier_count = 200;            % 有效子载波数
symbols_per_carrier = 16;       % 一帧内OFDM符号数（用于形成串流做频谱对比）
bits_per_symbol = 4;            % 16QAM
IFFT_bin_length = 512;          % IFFT 点数
PrefixRatio = 1/4;              % CP 比例
GI = PrefixRatio*IFFT_bin_length;   % 循环前缀长度
beta = 1/16;                    % 滚降系数（窗的过渡比例）
GIP = floor(beta*(IFFT_bin_length+GI)); % 后缀长度，用于窗的尾部过渡

%================= 频域映射（构造实信号） =================
baseband_out_length = carrier_count*symbols_per_carrier*bits_per_symbol;
rand('twister', 0);
tx_bits = randi([0 1], 1, baseband_out_length);

% 16QAM 调制
tx_symbols = qam16(tx_bits); % 依赖同目录下 qam16.m
tx_symbols_matrix = reshape(tx_symbols', carrier_count, symbols_per_carrier)';

% 进行子载波共轭映射，使 IFFT 输出为实信号
carriers = (1:carrier_count) + (floor(IFFT_bin_length/4) - floor(carrier_count/2));
conj_carriers = IFFT_bin_length - carriers + 2;

X = zeros(symbols_per_carrier, IFFT_bin_length);
X(:, carriers) = tx_symbols_matrix;
X(:, conj_carriers) = conj(tx_symbols_matrix);

%================= IFFT（单符号时域） =================
time_symbols = ifft(X, IFFT_bin_length, 2);   % 大小: [S, N]

%================= 添加CP和后缀（逐符号） =================
S = symbols_per_carrier; N = IFFT_bin_length; 
XX = zeros(S, N + GI + GIP);
for k = 1:S
    % 符号体写入（中间）
    XX(k, GI + (1:N)) = time_symbols(k, 1:N);
    % CP = 符号尾部复制到开头
    XX(k, 1:GI) = time_symbols(k, N-GI+1:N);
    % 后缀 = 符号头部复制到末尾（用于窗的右侧过渡）
    if GIP > 0
        XX(k, GI + N + (1:GIP)) = time_symbols(k, 1:GIP);
    end
end

%================= 加窗（仅边缘过渡位于CP/后缀内） =================
% 采用余弦滚降窗：长度为 N+GI，过渡段主要落在CP内；右侧后缀承接右边过渡
window_full = rcoswindow(beta, N + GI)';          % 依赖 rcoswindow.m，返回长度为 (1+beta)*(N+GI)+1
window_vec = window_full(1:(N+GI));               % 截取前 N+GI 个元素匹配符号体长度
windowed_symbols = XX;                             % 拷贝一份
for k = 1:S
    % 仅对 [1 : N+GI] 区段乘窗，尾部后缀不乘（已作为过渡延拓）
    windowed_symbols(k, 1:(N+GI)) = real(XX(k, 1:(N+GI))) .* window_vec;
    if GIP > 0
        windowed_symbols(k, (N+GI+1):(N+GI+GIP)) = real(XX(k, (N+GI+1):(N+GI+GIP)));
    end
end

%================= 形成串行发送信号（对比：未加窗 vs 加窗） =================
% 未加窗（逐符号含 CP+后缀，全保留串接）
tx_no_window_serial = reshape(XX', 1, []);

% 加窗：按步长 (N+GI) 串接，仅在整帧末尾保留一个后缀
tx_window_serial = zeros(1, S*(N+GI) + GIP);
tx_window_serial(1:N+GI+GIP) = windowed_symbols(1, :);
for i = 1:(S-1)
    head_idx = (N+GI)*i + 1;
    tail_idx = (N+GI)*(i+1) + GIP;
    tx_window_serial(head_idx:tail_idx) = windowed_symbols(i+1, :);
end

%================= 可视化 1：第二个符号（时域 + 单符号频域） =================
sym_idx = min(2, S);        % 取第2个符号（若S=1则取1）
t_axis = 0:(N+GI+GIP-1);

% 未加窗：时域
figure(12);
subplot(2,1,1);
plot(t_axis, real(XX(sym_idx, :)), 'b-'); grid on;
axis([0, N+GI+GIP, -0.25, 0.25]);
ylabel('Amplitude'); xlabel('Time (samples)');
title('Second OFDM Symbol (Before Windowing) - Time');

% 未加窗：单符号频域（幅度谱，0~fs/2）
subplot(2,1,2);
seq_before = real(XX(sym_idx, :));
Nb = length(seq_before);
spec_before = abs(fft(seq_before));
f_before = (0:(Nb-1))/Nb;
plot(f_before(1:floor(Nb/2)), spec_before(1:floor(Nb/2)), 'b-'); grid on;
xlabel('Normalized Frequency (0.5 = fs/2)'); ylabel('Magnitude');
title('Second OFDM Symbol (Before Windowing) - Spectrum');

% 加窗：时域
figure(13);
subplot(2,1,1);
plot(t_axis, windowed_symbols(sym_idx, :), 'r-'); grid on;
axis([0, N+GI+GIP, -0.25, 0.25]);
ylabel('Amplitude'); xlabel('Time (samples)');
title('Second OFDM Symbol (After Windowing) - Time');

% 加窗：单符号频域（幅度谱，0~fs/2）
subplot(2,1,2);
seq_after = windowed_symbols(sym_idx, :);
Na = length(seq_after);
spec_after = abs(fft(seq_after));
f_after = (0:(Na-1))/Na;
plot(f_after(1:floor(Na/2)), spec_after(1:floor(Na/2)), 'r-'); grid on;
xlabel('Normalized Frequency (0.5 = fs/2)'); ylabel('Magnitude');
title('Second OFDM Symbol (After Windowing) - Spectrum');

%================= 可视化 2：串行发送信号的平均功率谱（对比） =================
% 为降低频谱估计方差，进行分段平均（简易 Welch 思想）
avg_seg_symbols = max(4, ceil(S/4));                 % 用若干个符号长度做平均
seg_len = (N+GI+GIP)*avg_seg_symbols;                % 每段长度（未加窗串流）

% 未加窗串流：分段平均谱
L1 = min(seg_len, length(tx_no_window_serial));
num_avg1 = floor(length(tx_no_window_serial) / L1);
avg_fft1 = zeros(1, L1);
for a = 0:(num_avg1-1)
    blk = tx_no_window_serial((a*L1+1):((a+1)*L1));
    avg_fft1 = avg_fft1 + abs(fft(blk))/max(1, num_avg1);
end
avg_db1 = 20*log10(avg_fft1 + eps);

% 加窗串流：分段平均谱（段长与未加窗保持一致便于对比）
L2 = L1;
num_avg2 = floor(length(tx_window_serial) / L2);
avg_fft2 = zeros(1, L2);
for a = 0:(num_avg2-1)
    blk = tx_window_serial((a*L2+1):((a+1)*L2));
    avg_fft2 = avg_fft2 + abs(fft(blk))/max(1, num_avg2);
end
avg_db2 = 20*log10(avg_fft2 + eps);

% 频率只画 0~fs/2
f1 = (0:(L1-1))/L1; hf = 1:floor(L1/2);

figure(14);
plot(f1(hf), avg_db1(hf), 'b-', 'DisplayName', 'Before Window'); hold on;
plot(f1(hf), avg_db2(hf), 'r-', 'DisplayName', 'After Window'); grid on;
xlabel('Normalized Frequency (0.5 = fs/2)'); ylabel('Magnitude (dB)');
title('Averaged Spectrum of Serial OFDM Signal');
legend('show', 'Location', 'southwest');

%================= 提示 =================
% 观察：
% 1) Figure 12 vs 13（单符号）：加窗后主瓣略宽、旁瓣显著降低；
% 2) Figure 14（串流平均谱）：红线（加窗）旁瓣抑制更好，带外辐射更低；
% 3) 有效信息仍在各子载波幅相中，接收端在去CP后做FFT可恢复。

toc

