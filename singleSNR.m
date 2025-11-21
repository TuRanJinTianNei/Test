% ============================================================================
% 文件名: singleSNR.m
% 功能: LTE标准OFDM系统完整仿真程序（主程序）
% 描述: 
%   1. 实现完整的LTE标准OFDM系统链路：发送端（16QAM调制、IFFT、CP、加窗）→ 
%      信道（AWGN）→ 接收端（去CP、FFT、16QAM解调）
%   2. 生成13个Figure进行可视化分析，包括：
%      - 单符号分析与可视化（Figure 1-5）
%      - 帧级、串行信号与频谱分析（Figure 6-8）
%      - 信道影响与接收端评估（Figure 12-13）
%      - 解调判决与性能（Figure 9-11，包含BER-SNR曲线）
%   3. 计算BER-SNR曲线（10-30dB，步进2dB）
%   4. 在15dB SNR下进行详细分析，其他SNR仅用于BER计算
%   5. 严格遵循LTE标准参数：10MHz带宽，15kHz子载波间隔，常规CP
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
% 【系统参数配置 - LTE标准10MHz系统】
%===============================================================================
% 说明：以下参数严格遵循LTE标准配置

% ============================================================================
% 一、基础固定参数（与子帧数无关）
% ============================================================================
subcarrier_spacing = 15e3;         % 子载波间隔：15kHz（LTE基本间隔）
fs = 30.72e6;                      % 采样率：30.72MHz（对应10MHz系统带宽）
symbol_duration = 1/subcarrier_spacing;  % 符号有效时长：≈66.67μs
modulation_order = 16;             % 调制方式：16QAM
bits_per_symbol = 4;               % 每个子载波承载的比特数（4=16QAM）

% ============================================================================
% 二、子载波与频域参数（固定配置）
% ============================================================================
IFFT_bin_length = 2048;           % IFFT点数：2048（对应10MHz系统）
% total_subcarriers = 600;         % 总子载波数：600个（标准索引k: -299~300，不含DC）
% null_subcarriers = 97;           % 空子载波数：97个（含1个DC + 48个低频保护带 + 48个高频保护带）

% LTE标准子载波索引映射（MATLAB IFFT bin从1开始）
% 标准索引k: -299 ~ 300，其中k=0为DC
% MATLAB bin: k=0→bin=1, k>0→bin=k+1, k<0→bin=N+k+1
% 有效子载波：k = -251 ~ 252（共503个，不含DC，其中正频率252个+负频率251个共轭对称）
% 注意：为使IFFT输出为实信号，需要共轭对称映射（数据仅映射到正频率，负频率自动共轭）

% 将标准索引k映射到MATLAB IFFT bin索引
map_k_to_bin = @(k) (k==0).*1 + (k>0).*(k+1) + (k<0).*(IFFT_bin_length+k+1);

% 有效子载波：k = -251 ~ 252（不含DC）
% 为了共轭对称，将数据映射到正频率子载波（k=1到252），负频率自动共轭
% 使用k=1到252作为正频率数据子载波
k_positive = 1:252;                 % 正频率子载波索引（k=1到252，共252个）
k_negative = -251:-1;              % 负频率子载波索引（k=-251到-1，共251个）

% 正频率子载波的MATLAB bin索引（用于数据映射）
carriers_positive = arrayfun(map_k_to_bin, k_positive);  % bin=2到253

% 负频率子载波的MATLAB bin索引（用于共轭映射）
carriers_negative = arrayfun(map_k_to_bin, k_negative);   % bin=1798到2048

% 实际数据子载波数（正频率部分，负频率自动共轭）
carrier_count = length(k_positive);  % 252个数据子载波
% 总有效子载波数（含共轭）：252 + 251 = 503个（不含DC）
% 空子载波：1个DC + 48个低频保护带（k=-299~-252）+ 48个高频保护带（k=253~300）= 97个

% ============================================================================
% 三、帧结构及时域累计参数（快速仿真：10个符号）
% ============================================================================
% 注意：为了加快仿真速度，将符号数减少到10个
num_subframes = 1;                 % 子帧数：1个（仅用于快速仿真）
subframe_duration = 1e-3;          % 单个子帧时长：1ms
total_duration = num_subframes * subframe_duration;  % 总时长
slots_per_subframe = 2;            % 每个子帧含时隙数：2个（0.5ms/时隙）
total_slots = num_subframes * slots_per_subframe;    % 总时隙数
symbols_per_slot = 7;              % 单个时隙OFDM符号数（常规CP）：7个符号/时隙
total_symbols = 10;                % 总符号数：10个（快速仿真）

% ============================================================================
% 四、循环前缀（CP）参数（常规CP，快速仿真：10个符号）
% ============================================================================
% 长CP符号（时隙第1个符号）：时长5.208μs，采样点160
% 短CP符号（时隙第2~7个符号）：时长4.6875μs，采样点144
CP_long_samples = 160;              % 长CP采样点数
CP_short_samples = 144;            % 短CP采样点数

% 快速仿真：10个符号，按照LTE标准模式分配CP
% 前2个符号：长CP（对应2个时隙的第1个符号）
% 后8个符号：短CP（对应时隙的第2-7个符号）
CP_lengths = zeros(1, total_symbols);
CP_lengths(1:2) = CP_long_samples;      % 前2个符号：长CP
CP_lengths(3:total_symbols) = CP_short_samples;  % 后8个符号：短CP

num_long_CP_symbols = sum(CP_lengths == CP_long_samples);      % 长CP符号数
num_short_CP_symbols = sum(CP_lengths == CP_short_samples);    % 短CP符号数

% ============================================================================
% 五、加窗设计参数（快速仿真：10个符号）
% ============================================================================
% 窗函数类型：升余弦窗
% 长CP符号过渡长度：40采样点（160/4）
% 短CP符号过渡长度：36采样点（144/4）
window_transition_long = 40;      % 长CP符号的窗过渡长度
window_transition_short = 36;      % 短CP符号的窗过渡长度

% 创建窗过渡长度数组（按符号顺序）
window_transitions = zeros(1, total_symbols);
window_transitions(CP_lengths == CP_long_samples) = window_transition_long;
window_transitions(CP_lengths == CP_short_samples) = window_transition_short;

% 单符号总采样点数
symbol_samples_long = IFFT_bin_length + CP_long_samples;      % 2208（2048+160）
symbol_samples_short = IFFT_bin_length + CP_short_samples;    % 2192（2048+144）

% 10个符号总采样点数
total_samples = num_long_CP_symbols * symbol_samples_long + num_short_CP_symbols * symbol_samples_short;

% ============================================================================
% 六、数据传输相关参数
% ============================================================================
% 每符号有效数据比特数：252个有效子载波 × 4bit/子载波 = 1008bit/符号
bits_per_symbol_total = carrier_count * bits_per_symbol;
% 10个符号总传输比特数（无编码）：10符号 × 1008bit/符号 = 10,080 bit
total_bits = total_symbols * bits_per_symbol_total;

% 符号交界处幅度问题解决方案（LTE标准方法）
% ============================================================================
% LTE方法：使用升余弦窗平滑符号边缘，符号拼接时在CP范围内重叠相加
% - 每个符号结构：[CP | 主体(N) | 后缀(过渡长度)]
% - 重叠区域：当前符号的后缀与下一个符号的CP前过渡长度个样本重叠
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
baseband_out_length = total_bits;  % 总比特数
rng('default'); rng('shuffle');
baseband_out = randi([0 1], 1, baseband_out_length);

%------------------------------------------------------------------------------
% 步骤2: 16QAM调制
%------------------------------------------------------------------------------
complex_carrier_matrix = qam16(baseband_out);
% 重塑为 [总符号数 × 每符号子载波数] 的矩阵
complex_carrier_matrix = reshape(complex_carrier_matrix', carrier_count, total_symbols)';

%------------------------------------------------------------------------------
% 步骤3: 频域子载波映射（共轭对称，使IFFT输出为实信号）
%------------------------------------------------------------------------------
% 为使IFFT输出为实信号，频域必须满足共轭对称：
% - DC子载波（bin=1，k=0）必须是实数
% - 正频率子载波（bin=2到N/2+1）承载数据
% - 负频率子载波（bin=N/2+2到N）是正频率的共轭
% - X[k] = conj(X[N-k+2])，其中X[1]是DC，X[2]到X[N/2+1]是正频率

IFFT_modulation = zeros(total_symbols, IFFT_bin_length);
for sym_idx = 1:total_symbols
    % 将数据映射到正频率子载波（bin=2到253，对应k=1到252）
    IFFT_modulation(sym_idx, carriers_positive) = complex_carrier_matrix(sym_idx, 1:carrier_count);
    
    % 负频率子载波：放置正频率的共轭（bin=1798到2048，对应k=-251到-1）
    % 注意：需要反向映射，使得X[N-k+2] = conj(X[k])
    for i = 1:length(carriers_negative)
        % 找到对应的正频率索引（对称位置）
        pos_idx = length(carriers_positive) - i + 1;
        if pos_idx > 0 && pos_idx <= length(carriers_positive)
            IFFT_modulation(sym_idx, carriers_negative(i)) = conj(complex_carrier_matrix(sym_idx, pos_idx));
        end
    end
    
    % DC子载波（bin=1，对应k=0）设为0（实数）
    IFFT_modulation(sym_idx, 1) = 0;
    % 其他空子载波位置保持为0（已在初始化中设置）
end

%------------------------------------------------------------------------------
% 步骤4: IFFT变换（频域 → 时域）
%------------------------------------------------------------------------------
% Figure 1: IFFT各频点的幅度（展示激活的子载波分布）
% - 可直观看到载波分配与未用频点（幅度为0），理解频域输入结构
figure(1);
stem(0:IFFT_bin_length-1, abs(IFFT_modulation(2,1:IFFT_bin_length)),'b*-')
grid on
axis ([0 IFFT_bin_length -0.5 4.5]);
ylabel('Magnitude');
xlabel('IFFT Bin');
title('LTE OFDM子载波频率幅度（10MHz系统，504个有效子载波）');
 
% Figure 2: IFFT各频点的相位（度）
% - 展示LTE标准复数OFDM的相位分布（共轭对称）
figure(2);
plot(0:IFFT_bin_length-1, (180/pi)*angle(IFFT_modulation(2,1:IFFT_bin_length)), 'go')
hold on
stem(carriers_positive-1, (180/pi)*angle(IFFT_modulation(2,carriers_positive)),'b*-');
stem(carriers_negative-1, (180/pi)*angle(IFFT_modulation(2,carriers_negative)),'r*-');
axis ([0 IFFT_bin_length -200 +200])
grid on
ylabel('Phase (degrees)')
xlabel('IFFT Bin')
title('LTE OFDM子载波相位（蓝色：正频率，红色：负频率共轭）')
legend('所有子载波', '正频率数据子载波', '负频率共轭子载波', 'Location', 'best')
hold off

%------------------------------------------------------------------------------
% 步骤5: IFFT变换，得到时域OFDM符号（未加保护间隔）
%------------------------------------------------------------------------------
signal_after_IFFT = ifft(IFFT_modulation, IFFT_bin_length, 2);
time_wave_matrix = signal_after_IFFT;

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
% 步骤6: 添加循环前缀(CP)与后缀（LTE标准：区分长CP和短CP）
%------------------------------------------------------------------------------
% CP作用：对抗多径干扰，保持符号间正交性
% 后缀作用：用于加窗的平滑过渡，降低频谱旁瓣
% LTE标准：不同符号的CP长度不同（长CP=160，短CP=144）
% 注意：IFFT输出是实信号（由于共轭对称），取实部

% 初始化：每个符号的长度不同（长CP符号更长）
max_symbol_len = IFFT_bin_length + CP_long_samples + window_transition_long;
time_wave_matrix_cp = cell(total_symbols, 1);  % 使用cell数组存储不同长度的符号

for k = 1:total_symbols
    CP_len = CP_lengths(k);
    window_trans = window_transitions(k);
    symbol_len = IFFT_bin_length + CP_len + window_trans;
    
    % 初始化当前符号
    XX = zeros(1, symbol_len);
    
    % IFFT输出是实信号（由于共轭对称），取实部
    signal_real = real(signal_after_IFFT(k, :));
    
    % 符号主体部分（中间）
    XX(CP_len+1:CP_len+IFFT_bin_length) = signal_real;
    
    % 循环前缀：将符号尾部复制到开头
    XX(1:CP_len) = signal_real(IFFT_bin_length-CP_len+1:IFFT_bin_length);
    
    % 后缀：将符号头部复制到末尾（用于窗的右侧过渡）
    if window_trans > 0
        XX(CP_len+IFFT_bin_length+1:CP_len+IFFT_bin_length+window_trans) = signal_real(1:window_trans);
    end
    
    time_wave_matrix_cp{k} = XX;
end  

% Figure 4: 单个OFDM符号添加循环前缀(CP)与后缀后的时域波形
% - 左侧CP、右端后缀（用于窗口过渡），观察边沿延拓
% 显示一个长CP符号（第1个符号）和一个短CP符号（第2个符号）
figure(4);
subplot(2,1,1);
plot(0:length(time_wave_matrix_cp{1})-1, time_wave_matrix_cp{1});
axis([0, length(time_wave_matrix_cp{1}), -0.2, 0.2]);
grid on;
ylabel('Amplitude');
xlabel('Time (samples)');
title(sprintf('LTE OFDM时域信号（长CP符号，CP=%d，总长度=%d，实信号）', CP_lengths(1), length(time_wave_matrix_cp{1})));
subplot(2,1,2);
plot(0:length(time_wave_matrix_cp{2})-1, time_wave_matrix_cp{2});
axis([0, length(time_wave_matrix_cp{2}), -0.2, 0.2]);
grid on;
ylabel('Amplitude');
xlabel('Time (samples)');
title(sprintf('LTE OFDM时域信号（短CP符号，CP=%d，总长度=%d，实信号）', CP_lengths(2), length(time_wave_matrix_cp{2})));

%------------------------------------------------------------------------------
% 步骤7: OFDM符号加窗处理（LTE风格，区分长CP和短CP）
%------------------------------------------------------------------------------
% LTE原理：使用升余弦窗平滑符号边缘，抑制符号边缘的陡峭跳变，减少带外辐射
% 窗函数覆盖范围：[CP | 主体(N)]，总长度为 N+CP_len
% 窗函数在左边缘（CP开始）和右边缘（主体结束）提供平滑过渡
% 后缀部分（window_trans）用于与下一个符号的CP重叠，实现平滑拼接
windowed_time_wave_matrix_cp = cell(total_symbols, 1);

% 对每个符号加窗（根据CP长度使用不同的窗函数）
for i = 1:total_symbols
    CP_len = CP_lengths(i);
    window_trans = window_transitions(i);
    symbol_data = time_wave_matrix_cp{i};
    symbol_len = length(symbol_data);
    
    % 生成升余弦窗函数（覆盖CP+主体部分，长度为N+CP_len）
    % 注意：rcoswindow返回的长度是(1+beta)*(N+CP_len)，我们需要只取前N+CP_len个元素
    beta_local = window_trans / (IFFT_bin_length + CP_len);  % 局部滚降系数
    rcos_win_full = rcoswindow(beta_local, IFFT_bin_length + CP_len);
    rcos_win = rcos_win_full(1:IFFT_bin_length + CP_len)';  % 转置为行向量
    
    % 初始化加窗后的符号
    windowed_symbol = zeros(1, symbol_len);
    
    % 符号结构：[CP | 主体(N) | 后缀(window_trans)]
    % 窗函数应用于前 N+CP_len 个样本（CP+主体）
    % 注意：symbol_data已经是实信号，直接使用
    windowed_symbol(1:IFFT_bin_length + CP_len) = ...
        symbol_data(1:IFFT_bin_length + CP_len) .* rcos_win;
    
    % 后缀部分（window_trans个样本）：保持原值，用于与下一个符号的CP重叠
    % 注意：后缀部分不加窗，因为它会与下一个符号的CP重叠相加
    if window_trans > 0
        windowed_symbol(IFFT_bin_length + CP_len + 1:end) = ...
            symbol_data(IFFT_bin_length + CP_len + 1:end);
    end
    
    windowed_time_wave_matrix_cp{i} = windowed_symbol;
end

% Figure 5: 加窗后的单个OFDM符号时域波形（含CP/后缀）
% - LTE风格：升余弦窗平滑边沿，降低频谱旁瓣
% - 窗函数应用于CP+主体部分，后缀部分用于与下一个符号的CP重叠
figure(5)
subplot(2,1,1);
plot(0:length(windowed_time_wave_matrix_cp{1})-1, windowed_time_wave_matrix_cp{1}); 
axis([0, length(windowed_time_wave_matrix_cp{1}), -0.2, 0.2]);
grid on;
ylabel('Amplitude');
xlabel('Time (samples)');
title(sprintf('LTE OFDM时域信号（长CP加窗，过渡长度=%d，实信号）', window_transitions(1)));
subplot(2,1,2);
plot(0:length(windowed_time_wave_matrix_cp{2})-1, windowed_time_wave_matrix_cp{2}); 
axis([0, length(windowed_time_wave_matrix_cp{2}), -0.2, 0.2]);
grid on;
ylabel('Amplitude');
xlabel('Time (samples)');
title(sprintf('LTE OFDM时域信号（短CP加窗，过渡长度=%d，实信号）', window_transitions(2)));

%------------------------------------------------------------------------------
% 步骤8: 生成发送信号，并串变换（LTE风格重叠相加，支持不同CP长度）
%------------------------------------------------------------------------------
% LTE原理：符号拼接时在CP范围内重叠相加，实现平滑过渡
% 符号结构：[CP | 主体(N) | 后缀(window_trans)]
% 重叠机制：
%   - 当前符号的后缀(window_trans个样本)与下一个符号的CP前window_trans个样本重叠
%   - 重叠区域：在串行序列中，两个符号的幅度相加
%   - 重叠区域限制在CP内，不影响有效数据（主体部分）

% 计算总长度（考虑重叠）
Tx_data = [];
current_pos = 1;

for sym_idx = 1:total_symbols
    symbol_data = windowed_time_wave_matrix_cp{sym_idx};
    CP_len = CP_lengths(sym_idx);
    window_trans = window_transitions(sym_idx);
    symbol_len = length(symbol_data);
    
    if sym_idx == 1
        % 第一个符号：完整写入
        % 注意：symbol_data已经是实信号，直接使用
        Tx_data = [Tx_data, symbol_data];
        current_pos = length(Tx_data) + 1;
    else
    % 后续符号：重叠相加处理
        prev_window_trans = window_transitions(sym_idx-1);
        
        if prev_window_trans > 0 && window_trans > 0
            % 重叠区域：前一个符号的后缀 + 当前符号的CP前window_trans个样本
            overlap_len = min([prev_window_trans, window_trans, CP_len]);
            overlap_start = current_pos - overlap_len;
            overlap_end = current_pos - 1;
            
            % 前一个符号的后缀部分（已在Tx_data中）
            prev_suffix = Tx_data(overlap_start:overlap_end);
            
            % 当前符号的CP前overlap_len个样本（已经是实信号）
            curr_cp_prefix = symbol_data(1:overlap_len);
            
            % 重叠相加
            Tx_data(overlap_start:overlap_end) = prev_suffix + curr_cp_prefix;
            
            % 写入当前符号的剩余部分（CP剩余+主体+后缀）
            if overlap_len < symbol_len
                Tx_data = [Tx_data, symbol_data(overlap_len+1:end)];
            end
        else
            % 无重叠：直接追加
            Tx_data = [Tx_data, symbol_data];
        end
        
        current_pos = length(Tx_data) + 1;
    end
end

% 用于对比：未加窗信号的直接串接（不含重叠）
% 注意：time_wave_matrix_cp已经是实信号，直接使用
Tx_data_withoutwindow = [];
for sym_idx = 1:total_symbols
    Tx_data_withoutwindow = [Tx_data_withoutwindow, time_wave_matrix_cp{sym_idx}];
end

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
zoom_len = 5000; % 局部窗口长度（与 Figure 12 相同）
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
plot(subset_ofdm_no_window, 'b-', 'LineWidth', 1.5)
grid on
xlabel('样本点')
ylabel('幅度')
title(sprintf('未加窗信号时域 (window [%d:%d])', zs2, ze2))

subplot(2,1,2)
plot(subset_ofdm_window, 'r-', 'LineWidth', 1.5)
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
% 步骤13: 串并转换，去除循环前缀和后缀（LTE风格，支持不同CP长度）
%------------------------------------------------------------------------------
% LTE原理：接收端去除CP时只取主体部分，不受重叠影响
% - 发送端：符号结构为 [CP | 主体(N) | 后缀(window_trans)]，重叠区域在CP内
% - 接收端：去除CP和后缀，只提取主体部分（N个样本）进行FFT解调
% - 重叠区域不影响有效数据，因为重叠只在CP范围内，主体部分不受影响
% - 注意：不同符号的CP长度不同（长CP=160，短CP=144）

% 从串行序列中提取每个符号（考虑重叠和不同CP长度）
Rx_data_matrix = zeros(total_symbols, IFFT_bin_length);  % 只存储主体部分

current_read_pos = 1;
for sym_idx = 1:total_symbols
    CP_len = CP_lengths(sym_idx);
    window_trans = window_transitions(sym_idx);
    
    % 计算当前符号在串行序列中的位置
    % 第一个符号：从位置1开始
    % 后续符号：需要考虑前一个符号的长度和重叠
    if sym_idx == 1
        symbol_start = 1;
    else
        % 前一个符号的长度
        prev_CP_len = CP_lengths(sym_idx-1);
        prev_window_trans = window_transitions(sym_idx-1);
        prev_symbol_len = IFFT_bin_length + prev_CP_len + prev_window_trans;
        
        % 当前符号的起始位置 = 前一个符号的起始位置 + 前一个符号长度 - 重叠长度
        if prev_window_trans > 0 && window_trans > 0
            overlap_len = min([prev_window_trans, window_trans, CP_len]);
            symbol_start = current_read_pos + prev_symbol_len - overlap_len;
        else
            symbol_start = current_read_pos + prev_symbol_len;
        end
    end
    
    % 提取当前符号（含CP+主体+后缀）
    symbol_end = symbol_start + IFFT_bin_length + CP_len + window_trans - 1;
    if symbol_end <= length(Rx_data)
        symbol_with_cp = Rx_data(symbol_start:symbol_end);
        
        % 去除CP和后缀，只保留主体部分（N个样本）
        Rx_data_matrix(sym_idx, :) = symbol_with_cp(CP_len+1:CP_len+IFFT_bin_length);
        
        current_read_pos = symbol_start;
    else
        error('接收数据长度不足');
    end
end

% 接收到的时域信号（实信号，主体部分）
Rx_data_complex_matrix = Rx_data_matrix; 

%------------------------------------------------------------------------------
% 步骤14: FFT解调（时域 → 频域），提取子载波数据
%------------------------------------------------------------------------------
Y1=fft(Rx_data_complex_matrix,IFFT_bin_length,2);
% 从正频率子载波提取数据（负频率是共轭，不需要单独提取）
Rx_carriers=Y1(:,carriers_positive);
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
% 计算空子载波数（包含DC/保护带）
% 空子载波 = 总子载波 - 有效子载波（含共轭）- DC = 600 - 503 - 1 = 96个
% 但根据实际配置：低频保护带48个（k=-299~-252）+ 高频保护带48个（k=253~300）+ DC 1个 = 97个
null_subcarriers_calc = 97;  % 1个DC + 48个低频保护带 + 48个高频保护带

fprintf('\n==== Simulation Summary (SNR = %.2f dB) ====\n', targetSNRdB);
fprintf('Total symbols     : %d\n', total_symbols);
fprintf('Total bits        : %d\n', baseband_out_length);
fprintf('BER               : %.6g\n', ber_15dB);
fprintf('Tx power (var)    : %.6g\n', Tx_signal_power);
fprintf('Noise variance    : %.6g\n', noise_sigma);
fprintf('IFFT length (N)   : %d\n', IFFT_bin_length);
fprintf('Active carriers   : %d (正频率数据子载波)\n', carrier_count);
fprintf('Total active      : %d (含共轭对称)\n', carrier_count + length(k_negative));
fprintf('Null carriers     : %d (含DC和保护带)\n', null_subcarriers_calc);
fprintf('CP length (long)  : %d samples (200个符号)\n', CP_long_samples);
fprintf('CP length (short) : %d samples (1200个符号)\n', CP_short_samples);
fprintf('Window transition (long)  : %d samples\n', window_transition_long);
fprintf('Window transition (short) : %d samples\n', window_transition_short);
fprintf('Solution method   : %d (LTE风格：重叠相加，重叠区域限制在CP内)\n', solution_method);
fprintf('Signal type       : 实信号（频域共轭对称）\n');
fprintf('==============================\n\n');

%===============================================================================
% 【BER-SNR性能曲线计算】（已禁用，仅仿真15dB SNR）
%===============================================================================
% 说明：为了加快仿真速度，仅仿真15dB SNR，不计算BER-SNR曲线
% 如需计算BER-SNR曲线，请取消以下注释并恢复原代码

% fprintf('\n==== Computing BER-SNR Curve (10-30 dB, step 2 dB) ====\n');
% fprintf('注意：快速仿真模式，仅计算15dB SNR点\n');
% fprintf('==========================================================\n');

% Figure 11: 单点BER显示（15dB）
figure(11);
semilogy(targetSNRdB, ber_15dB, 'b-o', 'LineWidth', 2, 'MarkerSize', 10);
xlabel('SNR (dB)');
ylabel('BER');
title(sprintf('OFDM系统误码率（快速仿真：10个符号，SNR=%.1f dB）', targetSNRdB));
grid on
axis([targetSNRdB-5 targetSNRdB+5 max(ber_15dB*0.1, 1e-6) min(ber_15dB*10, 1)])

toc
%===============================================================================
% 文件结束
%===============================================================================











