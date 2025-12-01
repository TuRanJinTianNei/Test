% [脚本标识] test4.m
% - IFFT=1024，GI=72（固定CP长度），使用升余弦加窗（LTE风格）
% - 频域采用"共轭镜像"映射法，支持完整 BER-SNR 曲线计算
% ============================================================================
% 文件名: test4.m
% 功能: OFDM系统完整仿真程序（主程序）
% 描述: 
%   1. 实现完整的OFDM系统链路：发送端（16QAM调制、IFFT、CP、加窗）→ 
%      信道（AWGN）→ 接收端（去CP、FFT、16QAM解调）
%   2. 生成13个Figure进行可视化分析，包括：
%      - 单符号分析与可视化（Figure 1-5）
%      - 帧级、串行信号与频谱分析（Figure 6-8）
%      - 信道影响与接收端评估（Figure 12-13）
%      - 解调判决与性能（Figure 9-11，包含BER-SNR曲线）
%   3. 计算BER-SNR曲线（10-30dB，步进2dB）
%   4. 在15dB SNR下进行详细分析，其他SNR仅用于BER计算
%   5. 采用LTE风格升余弦加窗，符号间重叠相加（重叠区域限制在CP内）
% ============================================================================

tic
clc;
clear all;
close all;
%===================== Figure 索引说明 =====================
%
% 【一】单符号分析与可视化
%   Figure 1 : IFFT 各频点幅度
%   Figure 2 : IFFT 各频点相位
%   Figure 3 : 单符号时域（未加 CP/后缀）
%   Figure 4 : 单符号时域（加 CP 与后缀）
%   Figure 5 : 单符号时域（加窗后）
%
% 【二】帧级、串行信号与频谱分析
%   Figure 6 : 串行发送时域对比（上：逐符号；下：帧级+末尾后缀）
%   Figure 7 : 加窗前后局部窗口时域对比（与 Figure 12 相同范围）
%   Figure 8 : 加窗前后局部窗口频域对比（与 Figure 12 相同范围）
%   注：Figure 7 展示加窗前后信号在局部窗口内的时域波形，便于观察窗函数对信号边缘的平滑作用。
%       Figure 8 展示加窗前后信号在局部窗口内的幅度谱（dB值），用于对比加窗对旁瓣的抑制效果。
%       两者均使用与 Figure 12 相同的局部窗口进行分析。纵轴零点是任意参考。
%
% 【三】信道影响与接收端评估
%   Figure 12: 局部窗口时域对比（含SNR/MSE命令行打印）
%   Figure 13: 局部窗口双边谱对比
%
% 【四】解调判决与性能
%   Figure 9 : 接收端16QAM星座图
%   Figure 10: 比特流发/收前100位对比
%   Figure 11: BER-SNR曲线（10-30dB，步进2dB）
%
%==========================================================
%===============================================================================
% 【系统参数配置】
%===============================================================================
% 说明：以下参数定义了OFDM系统的核心配置，可根据实际需求调整
% 基本时间/频率参数（按要求设置）
subcarrier_spacing = 15e3;         % 子载波间隔：15 kHz
fs = 15.36e6;                      % 采样频率：15.36 MHz

carrier_count=300;                 % 有效数据子载波数（正频率数据，负频率为镜像）；2*300+1(DC)=601 占用
symbols_per_frame=50;              % 每帧 OFDM 符号数
total_symbols=100;                 % 总共要传输的 OFDM 符号数（改为100）
symbols_per_carrier=total_symbols; % 为兼容后续矩阵尺寸，沿用原变量表示"总符号数"
bits_per_symbol=4;                 % 每个子载波承载的比特数（4=16QAM）

% IFFT/FFT参数
IFFT_bin_length=1024;              % IFFT 点数（包含有效、空、镜像与DC）

% 保护间隔参数
PrefixRatio = 72/ IFFT_bin_length;  % 循环前缀比例（使 GI 固定为 72）
GI = 72;                            % 保护前缀长度（样本数，按要求修改为 72）

% 加窗参数（LTE风格）
% ============================================================================
% LTE下行链路OFDM时域符号加窗原理：
% 1. 目的：抑制符号边缘的陡峭跳变，减少带外辐射（频谱泄露）
% 2. 重叠机制：符号拼接时存在可控的重叠，重叠区域限制在循环前缀（CP）内
% 3. 不影响有效数据：重叠只在CP内，接收端去除CP时不受影响
beta=1/16;                         % 加窗滚降系数（过渡带占比，越小过渡越短）
GIP=beta*(IFFT_bin_length+GI);     % 右端后缀长度：配合窗的尾部过渡
GIP=min(GIP, GI);                  % LTE要求：后缀长度不超过CP长度，确保重叠区域在CP内
GIP=floor(GIP);                    % 确保为整数

% 符号交界处幅度问题解决方案（LTE标准方法）
% ============================================================================
% LTE方法：使用升余弦窗平滑符号边缘，符号拼接时在CP范围内重叠相加
% - 每个符号结构：[CP(GI) | 主体(N) | 后缀(GIP)]
% - 重叠区域：当前符号的后缀(GIP)与下一个符号的CP前GIP个样本重叠
% - 重叠相加：在重叠区域将两个符号的幅度相加，确保平滑过渡
% - 接收端：去除CP时只取主体部分，不受重叠影响
% ============================================================================
solution_method = 3;               % 使用LTE标准方法：重叠相加（方案3）

% 信道参数
targetSNRdB = 15;                  % 目标信噪比（dB），用于15dB下的详细分析

%===============================================================================
% 【发送端处理流程】
%===============================================================================

%------------------------------------------------------------------------------
% 步骤1: 随机比特流生成
%------------------------------------------------------------------------------
baseband_out_length=carrier_count*symbols_per_carrier*bits_per_symbol;
% 使用基于时间的随机种子，确保每次运行生成不同的随机数据
baseband_out = randi([0 1], 1, baseband_out_length);

%------------------------------------------------------------------------------
% 步骤2: 16QAM调制
%------------------------------------------------------------------------------
complex_carrier_matrix=qam16(baseband_out);
complex_carrier_matrix=reshape(complex_carrier_matrix',carrier_count,symbols_per_carrier)';

%------------------------------------------------------------------------------
% 步骤3: 频域子载波映射（构造埃尔米特共轭对称，使IFFT输出为实信号）
%------------------------------------------------------------------------------
% 子载波索引计算：进行子载波共轭映射，使得OFDM符号经过IFFT之后是实信号
carriers=(1:carrier_count)+(floor(IFFT_bin_length/4)-floor(carrier_count/2));
conjugate_carriers=IFFT_bin_length-carriers+2; 

IFFT_modulation=zeros(symbols_per_carrier,IFFT_bin_length);
IFFT_modulation(:,carriers)=complex_carrier_matrix;
IFFT_modulation(:,conjugate_carriers)=conj(complex_carrier_matrix);

%------------------------------------------------------------------------------
% 步骤4: IFFT变换（频域 → 时域）
%------------------------------------------------------------------------------
% Figure 1: IFFT各频点的幅度（展示激活的子载波与其共轭镜像的幅度分布）
% - 可直观看到载波分配与未用频点（幅度为0），理解频域输入结构
figure(1);
stem(0:IFFT_bin_length-1, abs(IFFT_modulation(2,1:IFFT_bin_length)),'b*-')
grid on
axis ([0 IFFT_bin_length -0.5 4.5]);
ylabel('Magnitude');
xlabel('IFFT Bin');
title('OFDM子载波频率幅度');
 
% Figure 2: IFFT各频点的相位（度）
% - 体现共轭映射带来的相位对称性，验证实值时域信号条件
figure(2);
plot(0:IFFT_bin_length-1, (180/pi)*angle(IFFT_modulation(2,1:IFFT_bin_length)), 'go')
hold on
stem(carriers-1, (180/pi)*angle(IFFT_modulation(2,carriers)),'b*-');
stem(conjugate_carriers-1, (180/pi)*angle(IFFT_modulation(2,conjugate_carriers)),'b*-');
axis ([0 IFFT_bin_length -200 +200])
grid on
ylabel('Phase (degrees)')
xlabel('IFFT Bin')
title('OFDM子载波相位')

%------------------------------------------------------------------------------
% 步骤5: IFFT变换，得到时域OFDM符号（未加保护间隔）
%------------------------------------------------------------------------------
signal_after_IFFT=ifft(IFFT_modulation,IFFT_bin_length,2);
time_wave_matrix=signal_after_IFFT;

% Figure 3: 单个OFDM符号（未加循环前缀/后缀）的时域波形
% - 展示IFFT输出的一个符号周期包络与幅度范围
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
% CP作用：对抗多径干扰，保持符号间正交性
% 后缀作用：用于加窗的平滑过渡，降低频谱旁瓣
XX=zeros(symbols_per_carrier,IFFT_bin_length+GI+GIP);
for k=1:symbols_per_carrier
    % 符号主体部分（中间）
    for i=1:IFFT_bin_length
        XX(k,i+GI)=signal_after_IFFT(k,i);
    end
    % 循环前缀：将符号尾部复制到开头
    for i=1:GI
        XX(k,i)=signal_after_IFFT(k,i+IFFT_bin_length-GI);
    end
    % 后缀：将符号头部复制到末尾（用于窗的右侧过渡）
    for j=1:GIP
         XX(k,IFFT_bin_length+GI+j)=signal_after_IFFT(k,j);
    end
end  
time_wave_matrix_cp=XX;

% Figure 4: 单个OFDM符号添加循环前缀(CP)与后缀后的时域波形
% - 左侧CP、右端后缀（用于窗口过渡），观察边沿延拓
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
% LTE原理：使用升余弦窗平滑符号边缘，抑制符号边缘的陡峭跳变，减少带外辐射
% 窗函数覆盖范围：[CP(GI) | 主体(N)]，总长度为 N+GI
% 窗函数在左边缘（CP开始）和右边缘（主体结束）提供平滑过渡
% 后缀部分（GIP）用于与下一个符号的CP重叠，实现平滑拼接
windowed_time_wave_matrix_cp=zeros(symbols_per_carrier,IFFT_bin_length+GI+GIP);

% 生成升余弦窗函数（覆盖CP+主体部分，长度为N+GI）
% 注意：rcoswindow返回的长度是(1+beta)*(N+GI)，我们需要只取前N+GI个元素
rcos_win_full = rcoswindow(beta, IFFT_bin_length+GI);  % 列向量，长度 = (1+beta)*(N+GI)
rcos_win = rcos_win_full(1:IFFT_bin_length+GI)';  % 只取前N+GI个元素，转置为行向量

% 对每个符号加窗
for i = 1:symbols_per_carrier
    % 符号结构：[CP(GI) | 主体(N) | 后缀(GIP)]
    % 窗函数应用于前 N+GI 个样本（CP+主体）
    windowed_time_wave_matrix_cp(i, 1:IFFT_bin_length+GI) = ...
        real(time_wave_matrix_cp(i, 1:IFFT_bin_length+GI)) .* rcos_win;
    
    % 后缀部分（GIP个样本）：保持原值，用于与下一个符号的CP重叠
    % 注意：后缀部分不加窗，因为它会与下一个符号的CP重叠相加
    if GIP > 0
        windowed_time_wave_matrix_cp(i, IFFT_bin_length+GI+1:IFFT_bin_length+GI+GIP) = ...
            real(time_wave_matrix_cp(i, IFFT_bin_length+GI+1:IFFT_bin_length+GI+GIP));
    end
end

% Figure 5: 加窗后的单个OFDM符号时域波形（含CP/后缀）
% - LTE风格：升余弦窗平滑边沿，降低频谱旁瓣
% - 窗函数应用于CP+主体部分，后缀部分用于与下一个符号的CP重叠
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
% LTE原理：符号拼接时在CP范围内重叠相加，实现平滑过渡
% 帧结构：每帧包含多个OFDM符号，仅在帧末尾保留一次后缀
% 
% 符号结构：[CP(GI) | 主体(N) | 后缀(GIP)]
% 重叠机制：
%   - 当前符号的后缀(GIP个样本)与下一个符号的CP前GIP个样本重叠
%   - 重叠区域：在串行序列中，两个符号的幅度相加
%   - 重叠区域限制在CP内（GIP <= GI），不影响有效数据（主体部分）
% 
% 串行序列结构：
%   [符号1: CP+主体+后缀] [符号2: CP+主体+后缀] ... [符号N: CP+主体+后缀]
%   其中：符号i的后缀与符号i+1的CP前GIP个样本重叠相加
num_frames = total_symbols / symbols_per_frame;
frame_len_CP_suffix = symbols_per_frame*(IFFT_bin_length+GI)+GIP; % 每帧串行长度（仅末尾一次后缀）

% 按帧构造：LTE风格重叠相加
Tx_data = zeros(1, num_frames*frame_len_CP_suffix);

write_offset = 0;
for f = 1:num_frames
    sym_start = (f-1)*symbols_per_frame + 1;
    sym_end   = f*symbols_per_frame;

    % 当前帧的加窗符号矩阵
    frame_windowed = windowed_time_wave_matrix_cp(sym_start:sym_end, :);

    % LTE风格重叠相加：构造帧级串行信号
    frame_serial_windowed = zeros(1, frame_len_CP_suffix);
    
    % 第一个符号：完整写入（包含CP+主体+后缀）
    frame_serial_windowed(1:IFFT_bin_length+GI+GIP) = frame_windowed(1, :);
    
    % 后续符号：重叠相加处理
    for i = 1:(symbols_per_frame-1)
        % 下一个符号在串行序列中的起始位置
        % 第i个符号占据：[1+(i-1)*(N+GI), i*(N+GI)+GIP]
        % 第i+1个符号应该从：i*(N+GI)+1 开始（与第i个符号的后缀重叠GIP个样本）
        next_symbol_start = i*(IFFT_bin_length+GI) + 1;
        next_symbol_end = (i+1)*(IFFT_bin_length+GI) + GIP;
        
        % LTE重叠机制：
        % 重叠区域：当前符号的后缀(GIP)与下一个符号的CP前GIP个样本重叠
        % 重叠位置：串行序列中的 [next_symbol_start, next_symbol_start+GIP-1]
        % 当前符号的后缀：frame_windowed(i, IFFT_bin_length+GI+1:IFFT_bin_length+GI+GIP)
        % 下一个符号的CP前GIP个样本：frame_windowed(i+1, 1:GIP)
        
        if GIP > 0 && next_symbol_end <= length(frame_serial_windowed)
            % 重叠区域：当前符号的后缀（已写入）+ 下一个符号的CP前GIP个样本（待写入）
            overlap_start = next_symbol_start;
            overlap_end = overlap_start + GIP - 1;
            
            % 当前符号的后缀部分（在frame_serial_windowed中已写入）
            current_suffix = frame_windowed(i, IFFT_bin_length+GI+1:IFFT_bin_length+GI+GIP);
            
            % 下一个符号的CP前GIP个样本
            next_cp_prefix = frame_windowed(i+1, 1:GIP);
            
            % LTE重叠相加：在重叠区域将两个符号的幅度相加
            frame_serial_windowed(overlap_start:overlap_end) = ...
                current_suffix + next_cp_prefix;
            
            % 非重叠部分：写入下一个符号的剩余部分（CP的剩余部分+主体+后缀）
            if overlap_end < next_symbol_end
                non_overlap_start = overlap_end + 1;
                non_overlap_end = next_symbol_end;
                % 下一个符号从GIP+1开始到结尾
                frame_serial_windowed(non_overlap_start:non_overlap_end) = ...
                    frame_windowed(i+1, GIP+1:IFFT_bin_length+GI+GIP);
            end
        else
            % 如果没有后缀（GIP=0），直接写入下一个符号
            if next_symbol_end <= length(frame_serial_windowed)
                frame_serial_windowed(next_symbol_start:next_symbol_end) = frame_windowed(i+1, :);
            end
        end
    end

    % 写入到整帧串行序列
    Tx_data(write_offset+1:write_offset+frame_len_CP_suffix) = frame_serial_windowed;
    write_offset = write_offset + frame_len_CP_suffix;
end

% 逐符号直接串接（含每符号后缀）：用于对比与"未加窗"频谱分析
Tx_data_withoutwindow = reshape(time_wave_matrix_cp', (total_symbols)*(IFFT_bin_length+GI+GIP), 1)';

% 两种串接方式对比（用于Figure 6）
temp_time_symbol = length(Tx_data_withoutwindow); % 每符号带后缀
temp_time_frame  = length(Tx_data);               % 每帧一次后缀

% Figure 6: 发送端串行时域信号（两种串接方式对比）
% - 子图(1): 将每个符号(含CP/后缀)直接串接
% - 子图(2): 仅在每帧末尾添加一次后缀的加窗串行
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

%------------------------------------------------------------------------------
% 步骤9: 局部窗口定义（用于详细分析，与 Figure 12 相同范围）
%------------------------------------------------------------------------------
% 计算十个时域符号的长度（用于扩大局部窗口范围）
symbol_length = IFFT_bin_length + GI + GIP;  % 单个符号长度
zoom_len = 10 * symbol_length;  % 局部窗口长度：约十个时域符号长度
L = length(Tx_data);
zs = max(1, floor(L/2) - floor(zoom_len/2));
ze = min(L, zs + zoom_len - 1);
% 为未加窗信号计算对应的窗口范围
L2 = length(Tx_data_withoutwindow);
zs2 = max(1, floor(L2/2) - floor(zoom_len/2));
ze2 = min(L2, zs2 + zoom_len - 1);

% 提取局部窗口时域信号
subset_ofdm_no_window = Tx_data_withoutwindow(zs2:ze2);  % 未加窗信号局部窗口
subset_ofdm_window = Tx_data(zs:ze);                      % 加窗信号局部窗口

% 计算频域幅度谱（双边谱，包含负频率）
Nfft_no_window = length(subset_ofdm_no_window);
Nfft_window = length(subset_ofdm_window);
subset_ofdm_f_no_window = fftshift(fft(subset_ofdm_no_window));
subset_ofdm_f_log_no_window = 20*log10(abs(subset_ofdm_f_no_window)/max(abs(subset_ofdm_f_no_window)) + eps);
subset_ofdm_f_window = fftshift(fft(subset_ofdm_window));
subset_ofdm_f_log_window = 20*log10(abs(subset_ofdm_f_window)/max(abs(subset_ofdm_f_window)) + eps);

%------------------------------------------------------------------------------
% Figure 7: 加窗前后局部窗口时域对比（与 Figure 12 相同窗口范围）
%------------------------------------------------------------------------------
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

%------------------------------------------------------------------------------
% Figure 8: 加窗前后局部窗口频域对比（与 Figure 12 相同窗口范围，双边谱）
%------------------------------------------------------------------------------
% 归一化频率轴：-0.5..0.5（对应 -fs/2..fs/2）
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
% 【信道传输：AWGN加性高斯白噪声信道】
%===============================================================================
% 噪声功率计算：根据目标SNR和发送信号功率计算噪声方差
Tx_signal_power=var(Tx_data); % 以加窗后的发送串行为基准计算功率
linear_SNR=10^(targetSNRdB/10);
noise_sigma=Tx_signal_power/linear_SNR;
noise_scale_factor=sqrt(noise_sigma);
noise=randn(1,length(Tx_data))*noise_scale_factor;
Rx_data=Tx_data+noise;

%------------------------------------------------------------------------------
% 步骤10: 局部窗口时域对比（加噪前后）
%------------------------------------------------------------------------------
% 注：窗口范围 (zs, ze) 已在 Figure 7/8 处定义，此处直接使用
figure(12);
subplot(2,1,1);
plot(zs:ze, Tx_data(zs:ze));
grid on;
ylabel('Amplitude (volts)');
xlabel('Sample index');
title(sprintf('发送信号（局部放大）[%d:%d]', zs, ze));

subplot(2,1,2);
plot(zs:ze, Rx_data(zs:ze));
grid on;
ylabel('Amplitude (volts)');
xlabel('Sample index');
title(sprintf('接收信号（局部放大）[%d:%d]', zs, ze));

%------------------------------------------------------------------------------
% 步骤11: 计算并打印局部窗口的 SNR 与 MSE
%------------------------------------------------------------------------------
tx_seg = Tx_data(zs:ze);
rx_seg = Rx_data(zs:ze);
noise_seg = rx_seg - tx_seg;

sig_power  = mean(tx_seg.^2);
noise_power = mean(noise_seg.^2);
snr_zoom_db = 10*log10(sig_power / max(noise_power, eps));
mse_zoom = noise_power; % 均方误差=噪声均方功率

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

% 归一化到自身最大值并转 dB
Tx_mag_db_z = 20*log10( abs(Tx_Fz)/max(abs(Tx_Fz)) + eps );
Rx_mag_db_z = 20*log10( abs(Rx_Fz)/max(abs(Rx_Fz)) + eps );

% 归一化频率轴：-0.5..0.5（对应 -fs/2..fs/2）
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
% LTE原理：接收端去除CP时只取主体部分，不受重叠影响
% - 发送端：符号结构为 [CP(GI) | 主体(N) | 后缀(GIP)]，重叠区域在CP内
% - 接收端：去除CP和后缀，只提取主体部分（N个样本）进行FFT解调
% - 重叠区域不影响有效数据，因为重叠只在CP范围内，主体部分不受影响
symbol_len = IFFT_bin_length+GI+GIP;  % 每个符号的总长度（含CP+主体+后缀）

Rx_data_matrix=zeros(total_symbols,symbol_len);
read_offset = 0;
for f = 1:num_frames
    frame_rx = Rx_data(read_offset+1:read_offset+frame_len_CP_suffix);
    for i = 1:symbols_per_frame
        global_sym_idx = (f-1)*symbols_per_frame + i;
        % 提取符号：从串行序列中提取每个符号（含CP+主体+后缀）
        start_idx = (i-1)*(IFFT_bin_length+GI) + 1;
        end_idx = i*(IFFT_bin_length+GI) + GIP;
        Rx_data_matrix(global_sym_idx,:) = frame_rx(start_idx:end_idx);
    end
    read_offset = read_offset + frame_len_CP_suffix;
end
% LTE原理：去除CP和后缀，只保留主体部分（N个样本）进行FFT解调
% 重叠区域在CP内，去除CP后不影响有效数据
Rx_data_complex_matrix=Rx_data_matrix(:,GI+1:IFFT_bin_length+GI); 

%------------------------------------------------------------------------------
% 步骤14: FFT解调（时域 → 频域），提取子载波数据
%------------------------------------------------------------------------------
Y1=fft(Rx_data_complex_matrix,IFFT_bin_length,2);
Rx_carriers=Y1(:,carriers);
Rx_mag=abs(Rx_carriers);
Rx_phase=angle(Rx_carriers);

% 极坐标转直角坐标
[M, N]=pol2cart(Rx_phase, Rx_mag); 
Rx_complex_carrier_matrix = complex(M, N);
% Figure 9: 接收端子载波星座图（IQ平面），SNR越高散点越集中于16QAM理想点
figure(9);
plot(Rx_complex_carrier_matrix,'*b');
axis([-4, 4, -4, 4]);
title('接收信号16QAM星座图')
grid on

%------------------------------------------------------------------------------
% 步骤15: 16QAM解调（最小欧氏距离判决）
%------------------------------------------------------------------------------
Rx_serial_complex_symbols = reshape(Rx_complex_carrier_matrix',size(Rx_complex_carrier_matrix, 1)*size(Rx_complex_carrier_matrix,2),1)' ; 
Rx_decoded_binary_symbols=demoduqam16(Rx_serial_complex_symbols);
baseband_in = Rx_decoded_binary_symbols;
% Figure 10: 比特流对比（前100比特）——上:发送，下:接收判决
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
% 统计：占用与保护子载波（按要求口径：占用=含DC，保护=其余）
occupied_including_dc = 2*carrier_count + 1;          % 601
guard_subcarriers = IFFT_bin_length - occupied_including_dc; % 423
null_subcarriers = guard_subcarriers; % 与“保护子载波”口径一致

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
fprintf('==============================\n\n');

%===============================================================================
% 【BER-SNR性能曲线计算】
%===============================================================================
% 说明：在10-30dB范围内，步进2dB，计算各SNR点的误码率
% 注意：每次循环都会生成新的随机噪声，确保每个SNR点的噪声样本都是独立的
fprintf('\n==== Computing BER-SNR Curve (10-30 dB, step 2 dB) ====\n');
SNR_range = 10:2:30;  % 信噪比范围
ber_results = zeros(size(SNR_range));  % 存储BER结果

for idx = 1:length(SNR_range)
    snr_dB = SNR_range(idx);
    
    %--------------------------------------------------------------------------
    % 对每个SNR点：重新计算噪声并添加到发送信号
    % 每次调用 randn() 都会生成新的随机噪声样本，确保每次运行都不同
    %--------------------------------------------------------------------------
    linear_SNR = 10^(snr_dB/10);
    noise_sigma_loop = Tx_signal_power/linear_SNR;
    noise_scale_factor_loop = sqrt(noise_sigma_loop);
    noise_loop = randn(1, length(Tx_data)) * noise_scale_factor_loop;
    Rx_data_loop = Tx_data + noise_loop;
    
    % 接收端串并转换（按帧）去除循环前缀和后缀（LTE风格）
    symbol_len_loop = IFFT_bin_length+GI+GIP;  % 每个符号的总长度（含CP+主体+后缀）
    
    Rx_data_matrix_loop = zeros(total_symbols, symbol_len_loop);
    read_offset = 0;
    for f = 1:num_frames
        frame_rx = Rx_data_loop(read_offset+1:read_offset+frame_len_CP_suffix);
        for i = 1:symbols_per_frame
            global_sym_idx = (f-1)*symbols_per_frame + i;
            % 提取符号：从串行序列中提取每个符号（含CP+主体+后缀）
            start_idx = (i-1)*(IFFT_bin_length+GI) + 1;
            end_idx = i*(IFFT_bin_length+GI) + GIP;
            Rx_data_matrix_loop(global_sym_idx,:) = frame_rx(start_idx:end_idx);
        end
        read_offset = read_offset + frame_len_CP_suffix;
    end
    % LTE原理：去除CP和后缀，只保留主体部分（N个样本）进行FFT解调
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
% 处理 BER 为 0 的情况，用很小的值代替以便在对数坐标下显示
ber_plot = ber_results;
ber_plot(ber_plot == 0) = 1e-10;  % 将 0 替换为很小的值

figure(11);
semilogy(SNR_range, ber_plot, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('SNR (dB)');
ylabel('BER');
title('OFDM系统误码率性能');
grid on

% 设置合理的轴范围
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
% 确保所有Figure窗口显示（非阻塞方式）
%===============================================================================
% 刷新所有图形窗口（不激活具体窗口，避免阻塞）
drawnow;  % 强制刷新所有图形窗口

fprintf('\n所有Figure窗口已创建。程序执行完成。\n');
%===============================================================================
% 文件结束
%===============================================================================











