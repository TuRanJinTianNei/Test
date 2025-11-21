% [脚本标识] test5.m
% - IFFT=1024，GI=72（固定CP长度），使用升余弦加窗（LTE风格）
% - 频域采用"共轭镜像"映射法，支持完整 BER-SNR 曲线计算
% - 集成Jakes信道模型，详细展示信道特征
% ============================================================================
% 文件名: test5.m
% 功能: OFDM系统完整仿真程序（主程序）- 集成Jakes信道
% 描述: 
%   1. 实现完整的OFDM系统链路：发送端（16QAM调制、IFFT、CP、加窗）→ 
%      信道（Jakes衰落+AWGN）→ 接收端（去CP、FFT、16QAM解调）
%   2. 生成多个Figure详细展示Jakes信道特征，包括：
%      - Jakes信道时域特征（Figure 14-17）
%      - Jakes信道频域特征（Figure 18-19）
%      - Jakes信道统计特性（Figure 20-23）
%      - OFDM系统性能（Figure 1-13，继承test4.m的图形）
%   3. 计算BER-SNR曲线（10-30dB，步进2dB）
%   4. 在15dB SNR下进行详细分析
% ============================================================================

tic
clc;
clear all;
close all;

%===================== Figure 索引说明 =====================
%
% 【一】单符号分析与可视化（继承test4.m）
%   Figure 1 : IFFT 各频点幅度
%   Figure 2 : IFFT 各频点相位
%   Figure 3 : 单符号时域（未加 CP/后缀）
%   Figure 4 : 单符号时域（加 CP 与后缀）
%   Figure 5 : 单符号时域（加窗后）
%
% 【二】帧级、串行信号与频谱分析（继承test4.m）
%   Figure 6 : 串行发送时域对比
%   Figure 7 : 加窗前后局部窗口时域对比
%   Figure 8 : 加窗前后局部窗口频域对比
%
% 【三】信道影响与接收端评估（继承test4.m）
%   Figure 12: 局部窗口时域对比（含SNR/MSE命令行打印）
%   Figure 13: 局部窗口双边谱对比
%
% 【四】解调判决与性能（继承test4.m）
%   Figure 9 : 接收端16QAM星座图
%   Figure 10: 比特流发/收前100位对比
%   Figure 11: BER-SNR曲线（10-30dB，步进2dB）
%
% 【五】Jakes信道特征分析（新增）
%   Figure 14: Jakes信道冲激响应（时域）
%   Figure 15: Jakes信道包络（幅度）随时间变化
%   Figure 16: Jakes信道相位随时间变化
%   Figure 17: Jakes信道实部和虚部随时间变化
%   Figure 18: Jakes信道频率响应（幅度）
%   Figure 19: Jakes信道频率响应（相位）
%   Figure 20: Jakes信道多普勒功率谱密度
%   Figure 21: Jakes信道包络分布（瑞利分布对比）
%   Figure 22: Jakes信道相位分布（均匀分布对比）
%   Figure 23: Jakes信道自相关函数
%
%==========================================================

%===============================================================================
% 【系统参数配置】
%===============================================================================
% 基本时间/频率参数
subcarrier_spacing = 15e3;         % 子载波间隔：15 kHz
fs = 15.36e6;                      % 采样频率：15.36 MHz

carrier_count=300;                 % 有效数据子载波数
symbols_per_frame=50;              % 每帧 OFDM 符号数
total_symbols=100;                 % 总共要传输的 OFDM 符号数
symbols_per_carrier=total_symbols;
bits_per_symbol=4;                 % 每个子载波承载的比特数（4=16QAM）

% IFFT/FFT参数
IFFT_bin_length=1024;              % IFFT 点数

% 保护间隔参数
PrefixRatio = 72/ IFFT_bin_length;
GI = 72;                            % 保护前缀长度

% 加窗参数（LTE风格）
beta=1/16;                         % 加窗滚降系数
GIP=beta*(IFFT_bin_length+GI);     % 右端后缀长度
GIP=min(GIP, GI);                  % LTE要求：后缀长度不超过CP长度
GIP=floor(GIP);                    % 确保为整数

solution_method = 3;               % 使用LTE标准方法：重叠相加

% 信道参数
targetSNRdB = 15;                  % 目标信噪比（dB）

% Jakes信道参数
% ============================================================================
% Jakes信道模型参数：
% - 用于模拟移动通信中的多径瑞利衰落
% - 假设有大量散射体，产生多个多径分量
% - 多普勒频移导致信道随时间变化
% ============================================================================
fd = 100;                          % 最大多普勒频移（Hz），典型值：100Hz（车速约54km/h@2GHz）
N0 = 34;                           % Jakes模型中的散射体数量（建议值：34，确保统计特性）
channel_length = 1000;            % 用于分析的信道样本长度（用于展示信道特征）

%===============================================================================
% 【发送端处理流程】（继承test4.m的结构）
%===============================================================================

%------------------------------------------------------------------------------
% 步骤1: 随机比特流生成
%------------------------------------------------------------------------------
baseband_out_length=carrier_count*symbols_per_carrier*bits_per_symbol;
% 使用基于时间的随机种子，确保每次运行生成不同的随机数据
rng('default');  % 重置随机数生成器为默认状态
rng(sum(100*clock),'twister');  % 基于当前时间设置随机种子，确保每次运行都不同
baseband_out = randi([0 1], 1, baseband_out_length);

%------------------------------------------------------------------------------
% 步骤2: 16QAM调制
%------------------------------------------------------------------------------
complex_carrier_matrix=qam16(baseband_out);
complex_carrier_matrix=reshape(complex_carrier_matrix',carrier_count,symbols_per_carrier)';

%------------------------------------------------------------------------------
% 步骤3: 频域子载波映射（构造埃尔米特共轭对称）
%------------------------------------------------------------------------------
carriers=(1:carrier_count)+(floor(IFFT_bin_length/4)-floor(carrier_count/2));
conjugate_carriers=IFFT_bin_length-carriers+2; 

IFFT_modulation=zeros(symbols_per_carrier,IFFT_bin_length);
IFFT_modulation(:,carriers)=complex_carrier_matrix;
IFFT_modulation(:,conjugate_carriers)=conj(complex_carrier_matrix);

%------------------------------------------------------------------------------
% 步骤4: IFFT变换（频域 → 时域）
%------------------------------------------------------------------------------
figure(1);
stem(0:IFFT_bin_length-1, abs(IFFT_modulation(2,1:IFFT_bin_length)),'b*-')
grid on
axis ([0 IFFT_bin_length -0.5 4.5]);
ylabel('Magnitude');
xlabel('IFFT Bin');
title('OFDM子载波频率幅度');
 
figure(2);
plot(0:IFFT_bin_length-1, (180/pi)*angle(IFFT_modulation(2,1:IFFT_bin_length)), 'go')
hold on
stem(0:carriers-1, (180/pi)*angle(IFFT_modulation(2,1:carriers)),'b*-');
stem(0:conjugate_carriers-1, (180/pi)*angle(IFFT_modulation(2,1:conjugate_carriers)),'b*-');
axis ([0 IFFT_bin_length -200 +200])
grid on
ylabel('Phase (degrees)')
xlabel('IFFT Bin')
title('OFDM子载波相位')

%------------------------------------------------------------------------------
% 步骤5: IFFT变换，得到时域OFDM符号
%------------------------------------------------------------------------------
signal_after_IFFT=ifft(IFFT_modulation,IFFT_bin_length,2);
time_wave_matrix=signal_after_IFFT;

figure(3)
plot(0:IFFT_bin_length-1,time_wave_matrix(2,:));
axis([0 IFFT_bin_length -0.2 0.2]);
grid on
ylabel('Amplitude');
xlabel('Time');
title('OFDM时域信号，单符号周期');

%------------------------------------------------------------------------------
% 步骤6: 添加循环前缀(CP)与后缀
%------------------------------------------------------------------------------
XX=zeros(symbols_per_carrier,IFFT_bin_length+GI+GIP);
for k=1:symbols_per_carrier
    for i=1:IFFT_bin_length
        XX(k,i+GI)=signal_after_IFFT(k,i);
    end
    for i=1:GI
        XX(k,i)=signal_after_IFFT(k,i+IFFT_bin_length-GI);
    end
    for j=1:GIP
         XX(k,IFFT_bin_length+GI+j)=signal_after_IFFT(k,j);
    end
end  
time_wave_matrix_cp=XX;

figure(4);
plot(0:size(time_wave_matrix_cp,2)-1, time_wave_matrix_cp(2,:));
axis([0, size(time_wave_matrix_cp,2), -0.2, 0.2]);
grid on;
ylabel('Amplitude');
xlabel('Time');
title('OFDM时域信号（含循环前缀），单符号周期');

%------------------------------------------------------------------------------
% 步骤7: OFDM符号加窗处理（LTE风格）
%------------------------------------------------------------------------------
windowed_time_wave_matrix_cp=zeros(symbols_per_carrier,IFFT_bin_length+GI+GIP);

rcos_win_full = rcoswindow(beta, IFFT_bin_length+GI);
rcos_win = rcos_win_full(1:IFFT_bin_length+GI)';

for i = 1:symbols_per_carrier
    windowed_time_wave_matrix_cp(i, 1:IFFT_bin_length+GI) = ...
        real(time_wave_matrix_cp(i, 1:IFFT_bin_length+GI)) .* rcos_win;
    
    if GIP > 0
        windowed_time_wave_matrix_cp(i, IFFT_bin_length+GI+1:IFFT_bin_length+GI+GIP) = ...
            real(time_wave_matrix_cp(i, IFFT_bin_length+GI+1:IFFT_bin_length+GI+GIP));
    end
end

figure(5)
plot(0:IFFT_bin_length-1+GI+GIP,windowed_time_wave_matrix_cp(2,:)); 
axis([0, IFFT_bin_length-1+GI+GIP, -0.2, 0.2]);
grid on;
ylabel('Amplitude');
xlabel('Time');
title(sprintf('OFDM时域信号（LTE风格加窗，重叠区域限制在CP内，GIP=%d），单符号周期', GIP));

%------------------------------------------------------------------------------
% 步骤8: 生成发送信号，并串变换（按帧组织，LTE风格重叠相加）
%------------------------------------------------------------------------------
num_frames = total_symbols / symbols_per_frame;
frame_len_CP_suffix = symbols_per_frame*(IFFT_bin_length+GI)+GIP;

Tx_data = zeros(1, num_frames*frame_len_CP_suffix);

write_offset = 0;
for f = 1:num_frames
    sym_start = (f-1)*symbols_per_frame + 1;
    sym_end   = f*symbols_per_frame;

    frame_windowed = windowed_time_wave_matrix_cp(sym_start:sym_end, :);

    frame_serial_windowed = zeros(1, frame_len_CP_suffix);
    
    frame_serial_windowed(1:IFFT_bin_length+GI+GIP) = frame_windowed(1, :);
    
    for i = 1:(symbols_per_frame-1)
        next_symbol_start = i*(IFFT_bin_length+GI) + 1;
        next_symbol_end = (i+1)*(IFFT_bin_length+GI) + GIP;
        
        if GIP > 0 && next_symbol_end <= length(frame_serial_windowed)
            overlap_start = next_symbol_start;
            overlap_end = overlap_start + GIP - 1;
            
            current_suffix = frame_windowed(i, IFFT_bin_length+GI+1:IFFT_bin_length+GI+GIP);
            next_cp_prefix = frame_windowed(i+1, 1:GIP);
            
            frame_serial_windowed(overlap_start:overlap_end) = ...
                current_suffix + next_cp_prefix;
            
            if overlap_end < next_symbol_end
                non_overlap_start = overlap_end + 1;
                non_overlap_end = next_symbol_end;
                frame_serial_windowed(non_overlap_start:non_overlap_end) = ...
                    frame_windowed(i+1, GIP+1:IFFT_bin_length+GI+GIP);
            end
        else
            if next_symbol_end <= length(frame_serial_windowed)
                frame_serial_windowed(next_symbol_start:next_symbol_end) = frame_windowed(i+1, :);
            end
        end
    end

    Tx_data(write_offset+1:write_offset+frame_len_CP_suffix) = frame_serial_windowed;
    write_offset = write_offset + frame_len_CP_suffix;
end

Tx_data_withoutwindow = reshape(time_wave_matrix_cp', (total_symbols)*(IFFT_bin_length+GI+GIP), 1)';

temp_time_symbol = length(Tx_data_withoutwindow);
temp_time_frame  = length(Tx_data);

figure (6)
subplot(2,1,1);
plot(0:temp_time_symbol-1,Tx_data_withoutwindow);
grid on
ylabel('Amplitude (volts)')
xlabel('Time (samples)')
title('OFDM时域信号')
temp_time2 = temp_time_frame;
subplot(2,1,2);
plot(0:temp_time2-1,Tx_data);
grid on
ylabel('Amplitude (volts)')
xlabel('Time (samples)')
title('OFDM时域信号')

%===============================================================================
% 【信号分析与可视化：局部窗口频谱与时域对比】
%===============================================================================

symbol_length = IFFT_bin_length + GI + GIP;
zoom_len = 10 * symbol_length;

L = length(Tx_data);
zs = max(1, floor(L/2) - floor(zoom_len/2));
ze = min(L, zs + zoom_len - 1);

L2 = length(Tx_data_withoutwindow);
zs2 = max(1, floor(L2/2) - floor(zoom_len/2));
ze2 = min(L2, zs2 + zoom_len - 1);

subset_ofdm_no_window = Tx_data_withoutwindow(zs2:ze2);
subset_ofdm_window = Tx_data(zs:ze);

Nfft_no_window = length(subset_ofdm_no_window);
Nfft_window = length(subset_ofdm_window);
subset_ofdm_f_no_window = fftshift(fft(subset_ofdm_no_window));
subset_ofdm_f_log_no_window = 20*log10(abs(subset_ofdm_f_no_window)/max(abs(subset_ofdm_f_no_window)) + eps);
subset_ofdm_f_window = fftshift(fft(subset_ofdm_window));
subset_ofdm_f_log_window = 20*log10(abs(subset_ofdm_f_window)/max(abs(subset_ofdm_f_window)) + eps);

figure (7)
subplot(2,1,1)
plot(real(subset_ofdm_no_window), 'b-', 'LineWidth', 1.5)
grid on
xlabel('样本点')
ylabel('幅度')
title(sprintf('未加窗信号时域 (window [%d:%d])', zs2, ze2))

subplot(2,1,2)
plot(real(subset_ofdm_window), 'r-', 'LineWidth', 1.5)
grid on
xlabel('样本点')
ylabel('幅度')
title(sprintf('加窗信号时域 (window [%d:%d])', zs, ze))

f_norm_no_window = (-Nfft_no_window/2:(Nfft_no_window/2-1))/Nfft_no_window;
f_norm_window = (-Nfft_window/2:(Nfft_window/2-1))/Nfft_window;

figure(8)
subplot(2,1,1)
plot(f_norm_no_window, subset_ofdm_f_log_no_window, 'b-', 'LineWidth', 1.5)
grid on
axis([-0.5 0.5 -40 max(subset_ofdm_f_log_no_window)])
ylabel('幅度 (dB)')
xlabel('归一化频率 (-0.5..0.5 = ±fs/2)')
title(sprintf('未加窗信号幅度谱（双边谱，window [%d:%d]）', zs2, ze2))

subplot(2,1,2)
plot(f_norm_window, subset_ofdm_f_log_window, 'r-', 'LineWidth', 1.5)
grid on
axis([-0.5 0.5 -40 max(subset_ofdm_f_log_window)])
ylabel('幅度 (dB)')
xlabel('归一化频率 (-0.5..0.5 = ±fs/2)')
title(sprintf('加窗信号幅度谱（双边谱，window [%d:%d]）', zs, ze))

%===============================================================================
% 【Jakes信道生成与特征分析】
%===============================================================================

%------------------------------------------------------------------------------
% Jakes信道生成函数
%------------------------------------------------------------------------------
% 生成Jakes信道冲激响应
% 输入：fd - 最大多普勒频移（Hz）
%       N0 - 散射体数量
%       N_samples - 样本数量
%       fs - 采样频率
% 输出：h - 复信道冲激响应（时域）
%------------------------------------------------------------------------------
fprintf('\n==== Generating Jakes Channel ====\n');
fprintf('Maximum Doppler frequency: %.2f Hz\n', fd);
fprintf('Number of scatterers: %d\n', N0);
fprintf('Channel samples: %d\n', channel_length);
fprintf('Sampling frequency: %.2e Hz\n', fs);

% 生成Jakes信道
dt = 1/fs;  % 采样间隔
t_channel = (0:channel_length-1) * dt;  % 时间轴

% Jakes模型：h(t) = sqrt(2/N0) * sum(exp(j*(2*pi*fd*cos(alpha_n)*t + phi_n)))
% 其中 alpha_n 和 phi_n 是随机相位
alpha_n = (2*pi*(0:N0-1))/N0;  % 等间隔角度
phi_n = 2*pi*rand(1, N0);       % 随机初始相位

% 生成Jakes信道（简化版本，使用余弦和正弦分量）
h_I = zeros(1, channel_length);  % 同相分量（实部）
h_Q = zeros(1, channel_length);  % 正交分量（虚部）

for n = 1:N0
    beta_n = pi * n / N0;  % Jakes模型中的角度
    h_I = h_I + cos(2*pi*fd*cos(beta_n)*t_channel + phi_n(n));
    h_Q = h_Q + sin(2*pi*fd*cos(beta_n)*t_channel + phi_n(n));
end

% 归一化
h_I = sqrt(2/N0) * h_I;
h_Q = sqrt(2/N0) * h_Q;
h_jakes = h_I + 1j*h_Q;  % 复信道冲激响应

% 归一化信道功率
h_jakes = h_jakes / sqrt(mean(abs(h_jakes).^2));

fprintf('Channel power (normalized): %.6f\n', mean(abs(h_jakes).^2));
fprintf('=====================================\n\n');

%------------------------------------------------------------------------------
% Figure 14: Jakes信道冲激响应（时域）
%------------------------------------------------------------------------------
figure(14);
subplot(2,1,1);
plot(t_channel*1000, abs(h_jakes), 'b-', 'LineWidth', 1.5);
grid on;
xlabel('时间 (ms)');
ylabel('幅度 |h(t)|');
title(sprintf('Jakes信道冲激响应幅度（fd=%.0f Hz）', fd));
legend('|h(t)|', 'Location', 'best');

subplot(2,1,2);
plot(t_channel*1000, angle(h_jakes)*180/pi, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('时间 (ms)');
ylabel('相位 (度)');
title(sprintf('Jakes信道冲激响应相位（fd=%.0f Hz）', fd));
legend('∠h(t)', 'Location', 'best');

%------------------------------------------------------------------------------
% Figure 15: Jakes信道包络（幅度）随时间变化
%------------------------------------------------------------------------------
figure(15);
plot(t_channel*1000, 20*log10(abs(h_jakes) + eps), 'b-', 'LineWidth', 1.5);
grid on;
xlabel('时间 (ms)');
ylabel('幅度 (dB)');
title(sprintf('Jakes信道包络随时间变化（fd=%.0f Hz）', fd));
legend('20log_{10}(|h(t)|)', 'Location', 'best');

%------------------------------------------------------------------------------
% Figure 16: Jakes信道相位随时间变化
%------------------------------------------------------------------------------
figure(16);
plot(t_channel*1000, unwrap(angle(h_jakes))*180/pi, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('时间 (ms)');
ylabel('相位 (度)');
title(sprintf('Jakes信道相位随时间变化（fd=%.0f Hz）', fd));
legend('∠h(t) (unwrapped)', 'Location', 'best');

%------------------------------------------------------------------------------
% Figure 17: Jakes信道实部和虚部随时间变化
%------------------------------------------------------------------------------
figure(17);
subplot(2,1,1);
plot(t_channel*1000, real(h_jakes), 'b-', 'LineWidth', 1.5);
hold on;
plot(t_channel*1000, imag(h_jakes), 'r-', 'LineWidth', 1.5);
grid on;
xlabel('时间 (ms)');
ylabel('幅度');
title(sprintf('Jakes信道实部和虚部（fd=%.0f Hz）', fd));
legend('实部 Re[h(t)]', '虚部 Im[h(t)]', 'Location', 'best');
hold off;

subplot(2,1,2);
plot(real(h_jakes), imag(h_jakes), 'b.', 'MarkerSize', 2);
grid on;
xlabel('实部 Re[h(t)]');
ylabel('虚部 Im[h(t)]');
title('Jakes信道复平面轨迹');
axis equal;

%------------------------------------------------------------------------------
% Figure 18: Jakes信道频率响应（幅度）
%------------------------------------------------------------------------------
Nfft_channel = 2^nextpow2(channel_length);
H_jakes = fft(h_jakes, Nfft_channel);
f_channel = (-Nfft_channel/2:(Nfft_channel/2-1)) * fs / Nfft_channel;

figure(18);
subplot(2,1,1);
plot(f_channel/1000, 20*log10(abs(fftshift(H_jakes)) + eps), 'b-', 'LineWidth', 1.5);
grid on;
xlabel('频率 (kHz)');
ylabel('幅度 (dB)');
title('Jakes信道频率响应（幅度）');
legend('|H(f)|', 'Location', 'best');

subplot(2,1,2);
plot(f_channel/1000, abs(fftshift(H_jakes)), 'b-', 'LineWidth', 1.5);
grid on;
xlabel('频率 (kHz)');
ylabel('幅度');
title('Jakes信道频率响应（幅度，线性刻度）');
legend('|H(f)|', 'Location', 'best');

%------------------------------------------------------------------------------
% Figure 19: Jakes信道频率响应（相位）
%------------------------------------------------------------------------------
figure(19);
plot(f_channel/1000, angle(fftshift(H_jakes))*180/pi, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('频率 (kHz)');
ylabel('相位 (度)');
title('Jakes信道频率响应（相位）');
legend('∠H(f)', 'Location', 'best');

%------------------------------------------------------------------------------
% Figure 20: Jakes信道多普勒功率谱密度
%------------------------------------------------------------------------------
% 计算多普勒功率谱密度（通过FFT的功率谱）
P_doppler = abs(fftshift(H_jakes)).^2;
f_doppler_normalized = f_channel / fd;  % 归一化到最大多普勒频移

figure(20);
subplot(2,1,1);
plot(f_doppler_normalized, 10*log10(P_doppler/max(P_doppler) + eps), 'b-', 'LineWidth', 1.5);
grid on;
xlabel('归一化频率 f/f_d');
ylabel('功率谱密度 (dB)');
title(sprintf('Jakes信道多普勒功率谱密度（fd=%.0f Hz）', fd));
legend('S(f)', 'Location', 'best');
axis([-2 2 -40 5]);

subplot(2,1,2);
% 理论Jakes多普勒谱：S(f) = 1/(π*fd*sqrt(1-(f/fd)^2)) for |f| < fd
f_theory = linspace(-fd*1.1, fd*1.1, 1000);
S_theory = zeros(size(f_theory));
idx_valid = abs(f_theory) < fd;
S_theory(idx_valid) = 1./(pi*fd*sqrt(1 - (f_theory(idx_valid)/fd).^2));
S_theory = S_theory / max(S_theory);  % 归一化

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

%------------------------------------------------------------------------------
% Figure 21: Jakes信道包络分布（瑞利分布对比）
%------------------------------------------------------------------------------
h_envelope = abs(h_jakes);
[counts, centers] = hist(h_envelope, 50);
pdf_hist = counts / (sum(counts) * (centers(2) - centers(1)));

% 理论瑞利分布：p(r) = (r/σ²) * exp(-r²/(2σ²))
sigma_rayleigh = sqrt(mean(h_envelope.^2) / 2);  % 瑞利分布参数
r_theory = linspace(0, max(h_envelope)*1.5, 1000);
pdf_rayleigh = (r_theory / sigma_rayleigh^2) .* exp(-r_theory.^2 / (2*sigma_rayleigh^2));

figure(21);
bar(centers, pdf_hist, 'FaceColor', [0.7 0.7 0.9], 'EdgeColor', 'none');
hold on;
plot(r_theory, pdf_rayleigh, 'r-', 'LineWidth', 2);
grid on;
xlabel('包络幅度 |h(t)|');
ylabel('概率密度');
title('Jakes信道包络分布（瑞利分布对比）');
legend('仿真直方图', '理论瑞利分布', 'Location', 'best');
hold off;

%------------------------------------------------------------------------------
% Figure 22: Jakes信道相位分布（均匀分布对比）
%------------------------------------------------------------------------------
h_phase = angle(h_jakes);
[counts_phase, centers_phase] = hist(h_phase, 50);
pdf_hist_phase = counts_phase / (sum(counts_phase) * (centers_phase(2) - centers_phase(1)));

% 理论均匀分布：p(θ) = 1/(2π) for -π ≤ θ ≤ π
theta_theory = linspace(-pi, pi, 1000);
pdf_uniform = ones(size(theta_theory)) / (2*pi);

figure(22);
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
% Figure 23: Jakes信道自相关函数
%------------------------------------------------------------------------------
% 计算自相关函数
max_lag = min(500, floor(channel_length/4));
R_h = xcorr(h_jakes, max_lag, 'normalized');
lags = (-max_lag:max_lag) * dt * 1000;  % 转换为毫秒

% 理论自相关函数：R(τ) = J0(2πfd*τ)，其中J0是零阶贝塞尔函数
tau_theory = linspace(0, max(lags), 1000);
R_theory = besselj(0, 2*pi*fd*tau_theory/1000);  % 注意单位转换

figure(23);
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
% 【信道传输：Jakes衰落 + AWGN】
%===============================================================================

% 将Jakes信道应用到发送信号
% 注意：需要将Jakes信道扩展到发送信号的长度
Tx_signal_power = var(Tx_data);
L_tx = length(Tx_data);

% 如果信道长度小于发送信号长度，需要扩展信道
if channel_length < L_tx
    % 重复或插值扩展信道
    h_extended = interp1(1:channel_length, h_jakes, ...
                         linspace(1, channel_length, L_tx), 'linear');
else
    h_extended = h_jakes(1:L_tx);
end

% 归一化扩展后的信道功率
h_extended = h_extended / sqrt(mean(abs(h_extended).^2));

% 通过Jakes信道（时域卷积，简化为一阶近似：乘法）
Rx_data_jakes = Tx_data .* h_extended;

% 添加AWGN噪声
linear_SNR = 10^(targetSNRdB/10);
noise_sigma = Tx_signal_power/linear_SNR;
noise_scale_factor = sqrt(noise_sigma);
noise = randn(1, L_tx) * noise_scale_factor;
Rx_data = Rx_data_jakes + noise;

%------------------------------------------------------------------------------
% 步骤10: 局部窗口时域对比（加噪前后）
%------------------------------------------------------------------------------
figure(12);
subplot(3,1,1);
plot(zs:ze, Tx_data(zs:ze));
grid on;
ylabel('Amplitude (volts)');
xlabel('Sample index');
title(sprintf('发送信号（局部放大）[%d:%d]', zs, ze));

subplot(3,1,2);
plot(zs:ze, abs(h_extended(zs:ze)));
grid on;
ylabel('|h(t)|');
xlabel('Sample index');
title(sprintf('Jakes信道包络（局部放大）[%d:%d]', zs, ze));

subplot(3,1,3);
plot(zs:ze, Rx_data(zs:ze));
grid on;
ylabel('Amplitude (volts)');
xlabel('Sample index');
title(sprintf('接收信号（Jakes衰落+AWGN，局部放大）[%d:%d]', zs, ze));

%------------------------------------------------------------------------------
% 步骤11: 计算并打印局部窗口的 SNR 与 MSE
%------------------------------------------------------------------------------
tx_seg = Tx_data(zs:ze);
rx_seg = Rx_data(zs:ze);
noise_seg = rx_seg - tx_seg;

sig_power = mean(tx_seg.^2);
noise_power = mean(noise_seg.^2);
snr_zoom_db = 10*log10(sig_power / max(noise_power, eps));
mse_zoom = noise_power;

fprintf('\n==== Zoomed Window Metrics (Figure 12) ====\n');
fprintf('Window range       : [%d : %d] (length=%d)\n', zs, ze, numel(tx_seg));
fprintf('Signal power       : %.6g\n', sig_power);
fprintf('Noise power        : %.6g\n', noise_power);
fprintf('SNR (Tx vs Rx)     : %.4f dB\n', snr_zoom_db);
fprintf('MSE (Tx vs Rx)     : %.6g\n', mse_zoom);
fprintf('===========================================\n');

%------------------------------------------------------------------------------
% 步骤12: Figure 13 - 与 Figure 12 同一窗口的频谱对比（双边谱）
%------------------------------------------------------------------------------
Nfft_zoom = 2^nextpow2(length(tx_seg));
Tx_Fz = fftshift(fft(tx_seg, Nfft_zoom));
Rx_Fz = fftshift(fft(rx_seg, Nfft_zoom));

Tx_mag_db_z = 20*log10(abs(Tx_Fz)/max(abs(Tx_Fz)) + eps);
Rx_mag_db_z = 20*log10(abs(Rx_Fz)/max(abs(Rx_Fz)) + eps);

f_norm_z = (-Nfft_zoom/2:(Nfft_zoom/2-1))/Nfft_zoom;

figure(13);
subplot(2,1,1);
plot(f_norm_z, Tx_mag_db_z, 'b');
grid on;
ylabel('Magnitude (dB)');
xlabel('Normalized Frequency (-0.5..0.5 = \pm fs/2)');
title(sprintf('发送信号频谱（双边谱，局部放大 [%d:%d]）', zs, ze));

subplot(2,1,2);
plot(f_norm_z, Rx_mag_db_z, 'r');
grid on;
ylabel('Magnitude (dB)');
xlabel('Normalized Frequency (-0.5..0.5 = \pm fs/2)');
title(sprintf('接收信号频谱（双边谱，局部放大 [%d:%d]）', zs, ze));

%===============================================================================
% 【接收端处理流程】
%===============================================================================

%------------------------------------------------------------------------------
% 步骤13: 串并转换，去除循环前缀和后缀（LTE风格）
%------------------------------------------------------------------------------
symbol_len = IFFT_bin_length+GI+GIP;

Rx_data_matrix=zeros(total_symbols,symbol_len);
read_offset = 0;
for f = 1:num_frames
    frame_rx = Rx_data(read_offset+1:read_offset+frame_len_CP_suffix);
    for i = 1:symbols_per_frame
        global_sym_idx = (f-1)*symbols_per_frame + i;
        start_idx = (i-1)*(IFFT_bin_length+GI) + 1;
        end_idx = i*(IFFT_bin_length+GI) + GIP;
        Rx_data_matrix(global_sym_idx,:) = frame_rx(start_idx:end_idx);
    end
    read_offset = read_offset + frame_len_CP_suffix;
end

Rx_data_complex_matrix=Rx_data_matrix(:,GI+1:IFFT_bin_length+GI); 

%------------------------------------------------------------------------------
% 步骤14: FFT解调（时域 → 频域），提取子载波数据
%------------------------------------------------------------------------------
Y1=fft(Rx_data_complex_matrix,IFFT_bin_length,2);
Rx_carriers=Y1(:,carriers);
Rx_mag=abs(Rx_carriers);
Rx_phase=angle(Rx_carriers);

[M, N]=pol2cart(Rx_phase, Rx_mag); 
Rx_complex_carrier_matrix = complex(M, N);

figure(9);
plot(Rx_complex_carrier_matrix,'*b');
axis([-4, 4, -4, 4]);
title('接收信号16QAM星座图（Jakes信道）')
grid on

%------------------------------------------------------------------------------
% 步骤15: 16QAM解调（最小欧氏距离判决）
%------------------------------------------------------------------------------
Rx_serial_complex_symbols = reshape(Rx_complex_carrier_matrix',size(Rx_complex_carrier_matrix, 1)*size(Rx_complex_carrier_matrix,2),1)'; 
Rx_decoded_binary_symbols=demoduqam16(Rx_serial_complex_symbols);
baseband_in = Rx_decoded_binary_symbols;

figure(10);
subplot(2,1,1);
stem(baseband_out(1:100));
title('发送比特流（前100位）')
subplot(2,1,2);
stem(baseband_in(1:100));
title('接收比特流（前100位）')

%------------------------------------------------------------------------------
% 步骤16: 误码率计算（15dB，用于显示摘要）
%------------------------------------------------------------------------------
bit_errors=find(baseband_in ~=baseband_out);
bit_error_count = size(bit_errors, 2); 
ber_15dB=bit_error_count/baseband_out_length;

%------------------------------------------------------------------------------
% 步骤17: 命令行输出关键指标
%------------------------------------------------------------------------------
occupied_including_dc = 2*carrier_count + 1;
guard_subcarriers = IFFT_bin_length - occupied_including_dc;
null_subcarriers = guard_subcarriers;

fprintf('\n==== Simulation Summary (SNR = %.2f dB) ====\n', targetSNRdB);
fprintf('Total symbols     : %d\n', total_symbols);
fprintf('Total bits         : %d\n', baseband_out_length);
fprintf('BER               : %.6g\n', ber_15dB);
fprintf('Tx power (var)    : %.6g\n', Tx_signal_power);
fprintf('Noise variance    : %.6g\n', noise_sigma);
fprintf('IFFT length (N)   : %d\n', IFFT_bin_length);
fprintf('Active carriers   : %d\n', carrier_count);
fprintf('Occupied (incl DC): %d\n', occupied_including_dc);
fprintf('Guard carriers    : %d\n', guard_subcarriers);
fprintf('CP length (GI)    : %d (ratio=%.3f)\n', GI, GI/IFFT_bin_length);
fprintf('Suffix length     : %d (GIP <= GI，限制在CP内)\n', GIP);
fprintf('One OFDM symbol   : %d samples (N+GI+GIP，LTE风格)\n', IFFT_bin_length+GI+GIP);
fprintf('Window roll-off   : 1/%d\n', round(1/beta));
fprintf('Solution method   : %d (LTE风格：重叠相加，重叠区域限制在CP内)\n', solution_method);
fprintf('Overlap length    : %d samples (GIP <= GI，限制在CP内)\n', GIP);
fprintf('Symbols per frame : %d\n', symbols_per_carrier);
fprintf('Jakes fd          : %.2f Hz\n', fd);
fprintf('Jakes scatterers  : %d\n', N0);
fprintf('==============================\n\n');

%===============================================================================
% 【BER-SNR性能曲线计算】
%===============================================================================
fprintf('\n==== Computing BER-SNR Curve (10-30 dB, step 2 dB) ====\n');
fprintf('Note: Each SNR point uses independent random noise\n');
SNR_range = 10:2:30;
ber_results = zeros(size(SNR_range));

for idx = 1:length(SNR_range)
    snr_dB = SNR_range(idx);
    
    % 对每个SNR点：重新计算噪声并添加到发送信号
    linear_SNR = 10^(snr_dB/10);
    noise_sigma_loop = Tx_signal_power/linear_SNR;
    noise_scale_factor_loop = sqrt(noise_sigma_loop);
    noise_loop = randn(1, length(Tx_data)) * noise_scale_factor_loop;
    Rx_data_loop = Rx_data_jakes + noise_loop;
    
    % 接收端串并转换（按帧）去除循环前缀和后缀（LTE风格）
    symbol_len_loop = IFFT_bin_length+GI+GIP;
    
    Rx_data_matrix_loop = zeros(total_symbols, symbol_len_loop);
    read_offset = 0;
    for f = 1:num_frames
        frame_rx = Rx_data_loop(read_offset+1:read_offset+frame_len_CP_suffix);
        for i = 1:symbols_per_frame
            global_sym_idx = (f-1)*symbols_per_frame + i;
            start_idx = (i-1)*(IFFT_bin_length+GI) + 1;
            end_idx = i*(IFFT_bin_length+GI) + GIP;
            Rx_data_matrix_loop(global_sym_idx,:) = frame_rx(start_idx:end_idx);
        end
        read_offset = read_offset + frame_len_CP_suffix;
    end
    
    Rx_data_complex_matrix_loop = Rx_data_matrix_loop(:,GI+1:IFFT_bin_length+GI);
    
    % FFT解调，提取子载波数据
    Y1_loop = fft(Rx_data_complex_matrix_loop, IFFT_bin_length, 2);
    Rx_carriers_loop = Y1_loop(:,carriers);
    Rx_mag_loop = abs(Rx_carriers_loop);
    Rx_phase_loop = angle(Rx_carriers_loop);
    [M_loop, N_loop] = pol2cart(Rx_phase_loop, Rx_mag_loop);
    Rx_complex_carrier_matrix_loop = complex(M_loop, N_loop);
    
    % 16QAM解调
    Rx_serial_complex_symbols_loop = reshape(Rx_complex_carrier_matrix_loop', size(Rx_complex_carrier_matrix_loop, 1)*size(Rx_complex_carrier_matrix_loop, 2), 1)';
    Rx_decoded_binary_symbols_loop = demoduqam16(Rx_serial_complex_symbols_loop);
    baseband_in_loop = Rx_decoded_binary_symbols_loop;
    
    % 计算BER
    bit_errors_loop = find(baseband_in_loop ~= baseband_out);
    bit_error_count_loop = size(bit_errors_loop, 2);
    ber_results(idx) = bit_error_count_loop / baseband_out_length;
    
    fprintf('SNR = %2d dB: BER = %.6g\n', snr_dB, ber_results(idx));
end
fprintf('==========================================================\n');

%------------------------------------------------------------------------------
% 步骤18: Figure 11 - BER-SNR性能曲线绘制
%------------------------------------------------------------------------------
ber_plot = ber_results;
ber_plot(ber_plot == 0) = 1e-10;

figure(11);
semilogy(SNR_range, ber_plot, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('SNR (dB)');
ylabel('BER');
title('OFDM系统误码率性能（Jakes信道）');
grid on

non_zero_ber = ber_results(ber_results > 0);
if ~isempty(non_zero_ber)
    y_min = min(non_zero_ber) * 0.5;
    y_max = max([ber_results, 1]) * 1.1;
else
    y_min = 1e-6;
    y_max = 1;
end
axis([9 31 y_min y_max])

toc
%===============================================================================
% 文件结束
%===============================================================================

