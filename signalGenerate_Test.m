%===============================================================================
% signalGenerate.m - OFDM信号生成函数（可独立运行）
% 
% 功能说明:
%   1. 生成完整的OFDM发送信号（基带实信号，通过共轭对称映射实现）
%   2. 支持16QAM调制、IFFT、循环前缀、升余弦加窗等完整流程
%   3. 采用LTE风格的符号重叠相加机制
%   4. 可独立运行，生成信号后保存在变量 Tx_data 中
%
% 依赖关系:
%   - 所有依赖函数已内置在文件末尾（qam16, rcoswindow, jakesChannel）
%   - 完全独立运行，无需外部函数文件
%   - 带宽估计部分：可选，需要调用 estimateBandwidth.m 函数
%     estimateBandwidth.m 依赖于 PSD_OFDM.m 及相关函数
%
% 使用示例:
%   1. 仅生成信号：
%      signalGenerate
%      信号保存在变量 Tx_data 中
%
%   2. 生成信号并估计带宽：
%      signalGenerate
%      [B_welch, B_ar, B_ideal, results] = estimateBandwidth(Tx_data, fs, ...
%          carrier_count, subcarrier_spacing, 'fc', 0, 'snr', 20, 'plot', true);
%
% 创建日期: 2025.12.10
% 基于: test3.m的信号生成部分
% 修改日期: 
%   2025.12.10 - 拆分带宽估计功能到 estimateBandwidth.m
%   2025.12.23 - 内置所有依赖函数（qam16, rcoswindow, jakesChannel），实现完全独立运行
%===============================================================================

clc;
clear all;

fprintf('========================================\n');
fprintf('OFDM信号生成程序（signalGenerate.m）\n');
fprintf('========================================\n\n');

%===============================================================================
% 【系统参数配置】- LTE 3MHz标准配置
%===============================================================================
% 基本时间/频率参数（LTE 3MHz标准）
subcarrier_spacing = 15e3;         % 子载波间隔：15 kHz（与20MHz系统一致）
fs = 3.84e6;                       % 采样频率：3.84 MHz（LTE 3MHz标准：15 kHz × 256）

carrier_count = 180;               % 有效数据子载波数：180个（正频率90个 + 负频率90个）
                                    % 注：LTE 3MHz标准为180个数据子载波 + 1个DC = 181个有效子载波
symbols_per_frame = 50;            % 每帧 OFDM 符号数
total_symbols = 100;                % 总共要传输的 OFDM 符号数
symbols_per_carrier = total_symbols; % 为兼容后续矩阵尺寸，沿用原变量表示"总符号数"
bits_per_symbol = 4;               % 每个子载波承载的比特数（4=16QAM）

% IFFT/FFT参数（LTE 3MHz标准）
IFFT_bin_length = 256;             % IFFT 点数：256（LTE 3MHz标准）

% 保护间隔参数（LTE 3MHz标准）
% LTE 3MHz正常CP：第一个符号20个采样点，后续符号18个采样点
% 这里使用后续符号的CP长度（18）作为默认值
GI_first_symbol = 20;              % 第一个符号的CP长度：20个采样点（≈5.2µs）
GI = 18;                            % 后续符号的CP长度：18个采样点（≈4.7µs）
PrefixRatio = GI / IFFT_bin_length; % 循环前缀比例：18/256 ≈ 7.03%

% 加窗参数（LTE风格）
beta = 1/16;                        % 加窗滚降系数（过渡带占比，越小过渡越短）
GIP = beta*(IFFT_bin_length+GI);    % 右端后缀长度：配合窗的尾部过渡
GIP = min(GIP, GI);                 % LTE要求：后缀长度不超过CP长度，确保重叠区域在CP内
GIP = floor(GIP);                   % 确保为整数

% 导频参数（可选，默认不使用）
use_pilot_equalization = false;     % 是否使用导频进行信道估计和均衡
pilot_spacing = 6;                  % 导频间隔（如果使用导频）
pilot_symbol = 1 + 1i*0;            % 导频符号（已知的复数值，归一化功率）

% 信道参数（高斯+衰落）
use_channel = true;                 % 是否通过信道传输（true=使用信道，false=仅生成发送信号）
use_rayleigh_fading = true;         % 是否使用瑞利衰落信道（true=使用，false=仅AWGN）
targetSNRdB = 15;                   % 目标信噪比（dB）
fd = 100;                           % 最大多普勒频移（Hz），典型值：100Hz（对应约30km/h@2.4GHz）
N0 = 34;                            % Jakes模型散射体数量（基础振荡器数量），典型值：34

% 输出参数
fprintf('========== 系统参数（LTE 3MHz标准）==========\n');
fprintf('信道带宽: 3.0 MHz\n');
fprintf('占用带宽: 2.7 MHz（90%%规则）\n');
fprintf('子载波间隔: %.1f kHz（与20MHz系统一致）\n', subcarrier_spacing/1e3);
fprintf('采样频率: %.2f MHz（15 kHz × 256）\n', fs/1e6);
fprintf('IFFT点数: %d\n', IFFT_bin_length);
fprintf('资源块数量: 15 RB（15 × 180 kHz = 2.7 MHz）\n');
fprintf('有效数据子载波数: %d（正频率90 + 负频率90 = 180）\n', carrier_count);
fprintf('有效子载波总数: 181（180数据 + 1个DC）\n');
fprintf('OFDM符号总数: %d\n', total_symbols);
fprintf('每帧符号数: %d\n', symbols_per_frame);
fprintf('循环前缀长度: %d（第一个符号）/ %d（后续符号）\n', GI_first_symbol, GI);
fprintf('后缀长度: %d\n', GIP);
fprintf('加窗滚降系数: 1/%d\n', round(1/beta));
if use_channel
    fprintf('信道传输: 启用\n');
    if use_rayleigh_fading
        fprintf('信道类型: Jakes瑞利衰落 + AWGN\n');
        fprintf('最大多普勒频移: %.1f Hz\n', fd);
    else
        fprintf('信道类型: AWGN only\n');
    end
    fprintf('目标SNR: %.1f dB\n', targetSNRdB);
else
    fprintf('信道传输: 禁用（仅生成发送信号）\n');
end
fprintf('================================\n\n');

%===============================================================================
% 【发送端处理流程】
%===============================================================================

%------------------------------------------------------------------------------
% 步骤1: 随机比特流生成
%------------------------------------------------------------------------------
% LTE 3MHz标准：180个数据子载波（正频率90个 + 负频率90个共轭对称）
% 实际需要生成90个子载波的数据（正频率部分），负频率通过共轭对称获得
positive_carrier_count = 90;  % 正频率子载波数（LTE 3MHz标准）

% 计算实际数据子载波数量（如果使用导频，需要排除导频占用的子载波）
if use_pilot_equalization
    % 计算导频位置（在正频率90个子载波中）
    pilot_positions_temp = 1:pilot_spacing:positive_carrier_count;
    data_carrier_count_actual = positive_carrier_count - length(pilot_positions_temp);
else
    data_carrier_count_actual = positive_carrier_count;  % 90个子载波
end

baseband_out_length = data_carrier_count_actual * symbols_per_carrier * bits_per_symbol;
% 使用基于时间的随机种子，确保每次运行生成不同的随机数据
rng('shuffle');  % 基于当前时间设置随机种子，确保每次运行都不同
baseband_out = randi([0 1], 1, baseband_out_length);

fprintf('步骤1: 生成随机比特流...\n');
fprintf('  - 比特流长度: %d bits\n', baseband_out_length);
fprintf('  - 数据子载波数: %d\n', data_carrier_count_actual);
fprintf('  完成！\n\n');

%------------------------------------------------------------------------------
% 步骤2: 16QAM调制（仅对数据子载波）
%------------------------------------------------------------------------------
fprintf('步骤2: 16QAM调制...\n');
% 调用内置的qam16函数（文件末尾定义）
complex_carrier_matrix = qam16(baseband_out);
complex_carrier_matrix = reshape(complex_carrier_matrix', data_carrier_count_actual, symbols_per_carrier)';
fprintf('  - 调制符号数: %d\n', length(complex_carrier_matrix(:)));
fprintf('  完成！\n\n');

%------------------------------------------------------------------------------
% 步骤2.5: 导频插入准备（计算导频位置和数据子载波位置）
%------------------------------------------------------------------------------
% LTE 3MHz标准子载波分配（256点FFT）：
%   [0]: DC子载波（置0）
%   [1~90]: 正频率部分（90个子载波）
%   [91~165]: 保护带（置零，75个子载波）
%   [166~255]: 负频率部分（90个子载波，共轭对称）
% 总有效子载波：180个数据子载波 + 1个DC = 181个
% 子载波索引计算：LTE 3MHz标准分配
% 修正：MATLAB FFT索引中，1是DC，2是第一个正频率
positive_carriers = 2:91;               % 正频率子载波：[2~91] (90个)
negative_carriers = 167:256;            % 负频率子载波：[167~256] (90个，共轭对称)
dc_carrier = 0;                         % DC子载波：[0] (MATLAB索引1)
guard_band = 92:166;                    % 保护带：[92~166] (75个)

% 为了兼容现有代码结构，将正频率和负频率合并
carriers = positive_carriers;       % 正频率子载波索引
conjugate_carriers = negative_carriers;  % 负频率子载波索引（共轭对称）

% 计算导频位置和数据子载波位置（LTE 3MHz标准）
% 注意：data_carrier_count_actual已经是90（正频率子载波数）
% complex_carrier_matrix的列数应该等于data_carrier_count_actual（90列）

if use_pilot_equalization
    % 在正频率子载波范围内插入导频
    pilot_positions_in_data = 1:pilot_spacing:data_carrier_count_actual;  % 导频在数据子载波中的相对位置
    pilot_carriers = positive_carriers(pilot_positions_in_data);  % 导频的绝对bin位置（正频率）
    pilot_conjugate_carriers = negative_carriers(pilot_positions_in_data);  % 导频的共轭位置（负频率）
    
    % 数据子载波位置（排除导频位置）
    data_positions_in_data = setdiff(1:data_carrier_count_actual, pilot_positions_in_data);
    data_carriers = positive_carriers(data_positions_in_data);  % 数据子载波的绝对bin位置（正频率）
    data_conjugate_carriers = negative_carriers(data_positions_in_data);  % 数据子载波的共轭位置（负频率）
    
    % 数据矩阵：complex_carrier_matrix已经是data_carrier_count_actual列
    data_matrix = complex_carrier_matrix(:, data_positions_in_data);  % 排除导频位置的数据
    pilot_matrix = repmat(pilot_symbol, symbols_per_carrier, length(pilot_positions_in_data));  % 导频矩阵
else
    % 不使用导频，所有子载波都是数据子载波
    % LTE 3MHz：正频率90个，负频率90个（共轭对称）
    data_carriers = positive_carriers;  % [2~91]
    data_conjugate_carriers = negative_carriers;  % [167~256]
    pilot_carriers = [];
    pilot_conjugate_carriers = [];
    
    % 数据矩阵：complex_carrier_matrix应该正好是90列（data_carrier_count_actual）
    data_matrix = complex_carrier_matrix;  % 直接使用，应该正好是90列
    pilot_matrix = [];
end

%------------------------------------------------------------------------------
% 步骤3: 频域子载波映射（LTE 3MHz标准分配）
%------------------------------------------------------------------------------
fprintf('步骤3: 频域子载波映射（LTE 3MHz标准）...\n');
IFFT_modulation = zeros(symbols_per_carrier, IFFT_bin_length);

% LTE 3MHz标准分配：
%   [1]: DC子载波（置零）
%   [2~91]: 正频率数据子载波
%   [92~166]: 保护带（置零）
%   [167~256]: 负频率数据子载波（共轭对称）

% 1. DC子载波置零
IFFT_modulation(:, dc_carrier + 1) = 0;  % MATLAB索引从1开始，所以是索引1

% 2. 映射正频率数据子载波 [2~91]
IFFT_modulation(:, data_carriers) = data_matrix;

% 3. 保护带置零 [92~166]
IFFT_modulation(:, guard_band) = 0;

% 4. 映射负频率数据子载波 [167~256]（共轭对称）
IFFT_modulation(:, data_conjugate_carriers) = conj(data_matrix);

% 5. 映射导频子载波（如果使用导频）
if use_pilot_equalization && ~isempty(pilot_matrix)
    IFFT_modulation(:, pilot_carriers) = pilot_matrix;
    IFFT_modulation(:, pilot_conjugate_carriers) = conj(pilot_matrix);
end

% 计算实际频率范围
freq_bin_spacing = fs / IFFT_bin_length;  % 频率bin间隔（Hz）
carrier_start_bin = min(data_carriers);  % 起始bin位置
carrier_end_bin = max(data_carriers);    % 结束bin位置
% 转换为相对于DC（0频率）的频率（Hz）
freq_start_hz = (carrier_start_bin - 1) * freq_bin_spacing;
freq_end_hz = (carrier_end_bin - 1) * freq_bin_spacing;
actual_bandwidth_hz = freq_end_hz - freq_start_hz;  % 实际占用带宽（Hz）

fprintf('  - LTE 3MHz标准子载波分配（256点FFT）：\n');
fprintf('    [1]: DC子载波（置0）\n');
fprintf('    [2~91]: 正频率数据子载波（%d个）\n', length(data_carriers));
fprintf('    [92~166]: 保护带（置零，75个子载波）\n');
fprintf('    [167~256]: 负频率数据子载波（%d个，共轭对称）\n', length(data_conjugate_carriers));
fprintf('  - 子载波bin范围: %d 到 %d（正频率）\n', carrier_start_bin, carrier_end_bin);
fprintf('  - 频率范围: %.3f MHz 到 %.3f MHz（相对于DC）\n', ...
    freq_start_hz/1e6, freq_end_hz/1e6);
fprintf('  - 占用带宽: %.3f MHz (%.0f kHz) - 15 RB × 180 kHz\n', ...
    actual_bandwidth_hz/1e6, actual_bandwidth_hz/1e3);
if use_pilot_equalization
    fprintf('  - 导频子载波位置: %d 个\n', length(pilot_carriers));
end
fprintf('  完成！\n\n');

%------------------------------------------------------------------------------
% 步骤4: IFFT变换（频域 → 时域）
%------------------------------------------------------------------------------
fprintf('步骤4: IFFT变换（频域 → 时域）...\n');
signal_after_IFFT = ifft(IFFT_modulation, IFFT_bin_length, 2);
time_wave_matrix = signal_after_IFFT;
fprintf('  - IFFT点数: %d\n', IFFT_bin_length);
fprintf('  完成！\n\n');

%------------------------------------------------------------------------------
% 步骤5: 添加循环前缀(CP)与后缀（LTE 3MHz标准）
%------------------------------------------------------------------------------
fprintf('步骤5: 添加循环前缀和后缀（LTE 3MHz标准）...\n');
% LTE 3MHz标准：第一个符号CP=20，后续符号CP=18
% 计算每个符号的CP长度
GI_per_symbol = zeros(1, symbols_per_carrier);
GI_per_symbol(1) = GI_first_symbol;  % 第一个符号：20个采样点
GI_per_symbol(2:end) = GI;            % 后续符号：18个采样点

% 计算最大符号长度（用于矩阵预分配）
max_symbol_length = IFFT_bin_length + GI_first_symbol + GIP;
XX = zeros(symbols_per_carrier, max_symbol_length);

for k = 1:symbols_per_carrier
    GI_current = GI_per_symbol(k);  % 当前符号的CP长度
    
    % 符号主体部分（中间）
    for i = 1:IFFT_bin_length
        XX(k, i+GI_current) = signal_after_IFFT(k, i);
    end
    % 循环前缀：将符号尾部复制到开头
    for i = 1:GI_current
        XX(k, i) = signal_after_IFFT(k, i+IFFT_bin_length-GI_current);
    end
    % 后缀：将符号头部复制到末尾（用于窗的右侧过渡）
    for j = 1:GIP
        XX(k, IFFT_bin_length+GI_current+j) = signal_after_IFFT(k, j);
    end
end
time_wave_matrix_cp = XX;
fprintf('  - 第一个符号CP长度: %d 个采样点（≈5.2µs）\n', GI_first_symbol);
fprintf('  - 后续符号CP长度: %d 个采样点（≈4.7µs）\n', GI);
fprintf('  - 后缀长度: %d\n', GIP);
fprintf('  - 符号总长度: %d（第一个符号）/%d（后续符号）\n', ...
    IFFT_bin_length+GI_first_symbol+GIP, IFFT_bin_length+GI+GIP);
fprintf('  完成！\n\n');

%------------------------------------------------------------------------------
% 步骤6: OFDM符号加窗处理（LTE风格）
%------------------------------------------------------------------------------
fprintf('步骤6: 升余弦加窗处理（LTE风格）...\n');
% 由于不同符号的CP长度不同，需要为每个符号单独处理
windowed_time_wave_matrix_cp = zeros(size(time_wave_matrix_cp));

% 对每个符号加窗
% 注意：由于共轭对称映射，IFFT后理论上应为实信号，但使用real()确保数值精度
for i = 1:symbols_per_carrier
    GI_current = GI_per_symbol(i);  % 当前符号的CP长度
    symbol_length = IFFT_bin_length + GI_current + GIP;  % 当前符号总长度
    
    % 生成升余弦窗函数（覆盖CP+主体+后缀）
    % rcoswindow生成长度为 round((1+beta)*Ts)+1，我们需要截取到 symbol_length
    Ts_current = IFFT_bin_length + GI_current;
    rcos_win_full = rcoswindow(beta, Ts_current);  % 列向量
    
    % 截取窗函数以匹配符号长度
    % 注意：rcoswindow生成的长度可能比symbol_length多1个点或正好
    if length(rcos_win_full) >= symbol_length
        rcos_win = rcos_win_full(1:symbol_length)';  % 转置为行向量
    else
        % 理论上不应发生，但作为防御
        rcos_win = [rcos_win_full; zeros(symbol_length-length(rcos_win_full), 1)]';
    end
    
    % 应用窗函数（覆盖整个符号：CP + 主体 + 后缀）
    windowed_time_wave_matrix_cp(i, 1:symbol_length) = ...
        real(time_wave_matrix_cp(i, 1:symbol_length)) .* rcos_win;
    
    % 之前的代码错误地截断了窗并手动复制了后缀，这里已修正为全长加窗
end
fprintf('  - 窗函数滚降系数: 1/%d\n', round(1/beta));
fprintf('  完成！\n\n');

%------------------------------------------------------------------------------
% 步骤7: 生成发送信号，并串变换（按帧组织，LTE风格重叠相加）
%------------------------------------------------------------------------------
fprintf('步骤7: 生成发送信号（LTE风格重叠相加）...\n');
% LTE原理：符号拼接时在CP范围内重叠相加，实现平滑过渡
% LTE 3MHz：第一个符号CP=20，后续符号CP=18，需要分别处理
% 检查total_symbols是否能被symbols_per_frame整除
if mod(total_symbols, symbols_per_frame) ~= 0
    error('错误：total_symbols (%d) 必须能被 symbols_per_frame (%d) 整除', ...
        total_symbols, symbols_per_frame);
end
num_frames = total_symbols / symbols_per_frame;

% 计算每帧长度：需要考虑不同符号的CP长度
% LTE 3MHz：第一个符号CP=20，后续符号CP=18
% 每帧长度 = 第一个符号长度 + (后续符号数-1)×后续符号长度 - (符号数-1)×重叠长度(GIP)
frame_len_CP_suffix = (IFFT_bin_length + GI_first_symbol + GIP) + ...
    (symbols_per_frame - 1) * (IFFT_bin_length + GI + GIP) - ...
    (symbols_per_frame - 1) * GIP;

% 按帧构造：LTE风格重叠相加
Tx_data = [];

for f = 1:num_frames
    sym_start = (f-1)*symbols_per_frame + 1;
    sym_end   = f*symbols_per_frame;

    % 当前帧的加窗符号矩阵
    frame_windowed = windowed_time_wave_matrix_cp(sym_start:sym_end, :);
    
    % 获取每个符号的实际长度
    symbol_lengths = zeros(1, symbols_per_frame);
    for i = 1:symbols_per_frame
        GI_frame = GI_per_symbol(sym_start + i - 1);
        symbol_lengths(i) = IFFT_bin_length + GI_frame + GIP;
    end

    % LTE风格重叠相加：构造帧级串行信号
    frame_serial_windowed = [];
    
    % 第一个符号：完整写入（包含CP+主体+后缀）
    frame_serial_windowed = frame_windowed(1, 1:symbol_lengths(1));
    
    % 后续符号：重叠相加处理
    for i = 1:(symbols_per_frame-1)
        GI_prev = GI_per_symbol(sym_start + i - 1);  % 前一个符号的CP长度
        GI_next = GI_per_symbol(sym_start + i);      % 当前符号的CP长度
        
        % 前一个符号的后缀部分
        prev_suffix_start = IFFT_bin_length + GI_prev + 1;
        prev_suffix_end = symbol_lengths(i);
        prev_suffix = frame_windowed(i, prev_suffix_start:prev_suffix_end);
        
        % 当前符号的CP前GIP个样本
        next_cp_prefix = frame_windowed(i+1, 1:GIP);
        
        % LTE重叠相加：在重叠区域将两个符号的幅度相加
        if GIP > 0
            overlap_region = prev_suffix + next_cp_prefix;
            % 更新重叠区域
            frame_serial_windowed(end-GIP+1:end) = overlap_region;
        end
        
        % 添加当前符号的剩余部分（CP的剩余部分+主体+后缀）
        if GI_next > GIP
            next_remaining_start = GIP + 1;
            next_remaining_end = symbol_lengths(i+1);
            frame_serial_windowed = [frame_serial_windowed, ...
                frame_windowed(i+1, next_remaining_start:next_remaining_end)];
        else
            % 如果CP长度小于等于GIP，直接添加主体+后缀
            frame_serial_windowed = [frame_serial_windowed, ...
                frame_windowed(i+1, GI_next+1:symbol_lengths(i+1))];
        end
    end

    % 写入到整帧串行序列
    Tx_data = [Tx_data, frame_serial_windowed];
end

fprintf('  - 发送信号长度: %d 样本\n', length(Tx_data));
fprintf('  - 信号时长: %.4f ms\n', length(Tx_data)/fs*1e3);
fprintf('  完成！\n\n');

%===============================================================================
% 【信号参数总结】
%===============================================================================
fprintf('========================================\n');
fprintf('信号生成完成！\n');
fprintf('========================================\n');
fprintf('发送信号参数：\n');
fprintf('  - 信号长度: %d 样本\n', length(Tx_data));
fprintf('  - 采样频率: %.2f MHz\n', fs/1e6);
fprintf('  - 信号时长: %.4f ms\n', length(Tx_data)/fs*1e3);
B_ideal = carrier_count * subcarrier_spacing;  % 占用带宽（Hz）：180 × 15 kHz = 2.7 MHz
B_channel = 3.0e6;                             % 信道带宽：3.0 MHz
fprintf('  - 信道带宽: %.1f MHz\n', B_channel/1e6);
fprintf('  - 占用带宽: %.3f MHz (%.0f kHz) - 180个子载波 × 15 kHz\n', B_ideal/1e6, B_ideal/1e3);
fprintf('  - 实际占用频率范围: %.3f MHz 到 %.3f MHz（正频率部分）\n', ...
    freq_start_hz/1e6, freq_end_hz/1e6);
fprintf('  - 实际占用带宽: %.3f MHz (%.0f kHz)\n', ...
    actual_bandwidth_hz/1e6, actual_bandwidth_hz/1e3);
fprintf('  - OFDM符号数: %d\n', total_symbols);
fprintf('  - 符号长度: %d样本（第一个，CP=20）/%d样本（后续，CP=18）\n', ...
    IFFT_bin_length+GI_first_symbol+GIP, IFFT_bin_length+GI+GIP);
fprintf('  - 有用符号时长(Tu): %.2f µs（1/15 kHz）\n', 1/subcarrier_spacing*1e6);
fprintf('  - 采样间隔(Ts): %.2f ns（1/%.2f MHz）\n', 1/fs*1e9, fs/1e6);
fprintf('========================================\n\n');

%===============================================================================
% 【可视化（可选）】
%===============================================================================
plot_signal = true;  % 是否绘制信号

if plot_signal
    fprintf('正在绘制信号波形和频谱...\n');
    
    % 计算带宽（用于显示）
    B_ideal = carrier_count * subcarrier_spacing;  % 占用带宽（Hz）：2.7 MHz
    B_channel = 3.0e6;                             % 信道带宽：3.0 MHz
    
    % 绘制时域信号（前1000个样本）
    figure('Name', 'OFDM基带发送信号（LTE 3MHz标准）', 'Position', [100, 100, 1200, 600]);
    subplot(2, 1, 1);
    plot_samples = min(1000, length(Tx_data));
    plot(1:plot_samples, Tx_data(1:plot_samples), 'b-', 'LineWidth', 1);
    grid on;
    xlabel('样本索引');
    ylabel('幅度');
    title(sprintf('LTE 3MHz基带发送信号时域波形（前%d个样本）\n信号类型：实信号（共轭对称映射）', plot_samples));
    
    % 绘制频谱
    subplot(2, 1, 2);
    Nfft = 2^nextpow2(length(Tx_data));
    Tx_Fz = fftshift(fft(Tx_data, Nfft));
    f_axis = (-Nfft/2:(Nfft/2-1)) * fs / Nfft;
    
    plot(f_axis/1e6, 20*log10(abs(Tx_Fz) / max(abs(Tx_Fz)) + eps));
    grid on;
    xlabel('频率 (MHz)');
    ylabel('幅度 (dB)');
    title(sprintf('LTE 3MHz基带发送信号频谱（归一化，双边频谱）\n信道带宽：%.1f MHz | 占用带宽：%.3f MHz', ...
        B_channel/1e6, B_ideal/1e6));
    xlim([-fs/2/1e6, fs/2/1e6]);
    hold on;
    % 标记实际频率范围
    plot([freq_start_hz/1e6, freq_start_hz/1e6], ylim, 'r--', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('起始频率 %.3f MHz', freq_start_hz/1e6));
    plot([freq_end_hz/1e6, freq_end_hz/1e6], ylim, 'r--', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('结束频率 %.3f MHz', freq_end_hz/1e6));
    plot([-freq_end_hz/1e6, -freq_end_hz/1e6], ylim, 'r--', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('负频率边界 -%.3f MHz', freq_end_hz/1e6));
    plot([-freq_start_hz/1e6, -freq_start_hz/1e6], ylim, 'r--', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('负频率边界 -%.3f MHz', freq_start_hz/1e6));
    legend('Location', 'best');
    hold off;
    
    fprintf('图形绘制完成！\n');
    fprintf('  Figure 1 说明（LTE 3MHz标准）：\n');
    fprintf('    - 信号类型：OFDM基带发送信号（Tx_data）\n');
    fprintf('    - 信号特性：实信号（通过共轭对称映射实现）\n');
    fprintf('    - 频谱类型：双边频谱（关于0频率对称）\n');
    fprintf('    - 信道带宽：%.1f MHz\n', B_channel/1e6);
    fprintf('    - 占用带宽：%.3f MHz (%.0f kHz) - 180个子载波 × 15 kHz\n', B_ideal/1e6, B_ideal/1e3);
    fprintf('    - 实际占用频率范围：%.3f MHz 到 %.3f MHz（正频率部分）\n', ...
        freq_start_hz/1e6, freq_end_hz/1e6);
    fprintf('    - 子载波分配：[2~91]正频率，[92~166]保护带，[167~256]负频率\n');
    fprintf('    - IFFT点数：256，采样频率：%.2f MHz\n', fs/1e6);
    fprintf('\n');

    % 绘制子载波映射图
    figure('Name', 'OFDM子载波映射（LTE 3MHz标准）', 'Position', [100, 100, 1000, 400]);
    
    % 准备绘图数据
    carriers_plot = zeros(1, IFFT_bin_length);
    
    % 标记数据子载波（正频率和负频率）
    carriers_plot(data_carriers) = 1;
    carriers_plot(data_conjugate_carriers) = 1;
    
    % 标记导频子载波（如果有）
    if use_pilot_equalization
        carriers_plot(pilot_carriers) = 1.5;
        carriers_plot(pilot_conjugate_carriers) = 1.5;
    end
    
    % 绘制
    stem(1:IFFT_bin_length, carriers_plot, 'b.', 'MarkerSize', 4);
    hold on;
    
    % 标记DC子载波
    plot(dc_carrier+1, 0, 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'DC子载波');
    
    % 标记保护带区域
    % 保护带是连续的区域 [92~166]
    fill([min(guard_band) max(guard_band) max(guard_band) min(guard_band)], ...
         [-0.1 -0.1 1.6 1.6], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', '保护带（空子载波）');
         
    % 设置坐标轴
    xlim([1 IFFT_bin_length]);
    ylim([-0.2 1.8]);
    xlabel('子载波索引 (1-256)');
    ylabel('类型');
    title(sprintf('LTE 3MHz 子载波映射 (256点FFT, %d个有效子载波)', carrier_count+1));
    
    % 自定义Y轴刻度
    set(gca, 'YTick', [0 1 1.5]);
    if use_pilot_equalization
        set(gca, 'YTickLabel', {'Null/DC', 'Data', 'Pilot'});
    else
        set(gca, 'YTickLabel', {'Null/DC', 'Data', ''});
    end
    
    % 添加图例
    % 创建用于图例的虚拟对象
    h_data = plot(nan, nan, 'b.', 'DisplayName', '数据子载波');
    h_dc = plot(nan, nan, 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'DC子载波');
    h_guard = fill(nan, nan, 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', '保护带');
    
    if use_pilot_equalization
        h_pilot = plot(nan, nan, 'b.', 'MarkerSize', 8, 'DisplayName', '导频子载波');
        legend([h_data, h_pilot, h_dc, h_guard], 'Location', 'best');
    else
        legend([h_data, h_dc, h_guard], 'Location', 'best');
    end
    
    grid on;
    fprintf('子载波映射图绘制完成！\n');
end

fprintf('所有步骤完成！\n');
fprintf('发送信号已保存在变量 Tx_data 中\n');

%===============================================================================
% 【信道传输：瑞利衰落信道 + AWGN加性高斯白噪声信道】（可选）
%===============================================================================
if use_channel
    fprintf('\n========================================\n');
    fprintf('开始信道传输...\n');
    fprintf('========================================\n\n');
    
    %------------------------------------------------------------------------------
    % 步骤1: 生成Jakes瑞利衰落信道（如果启用）
    %------------------------------------------------------------------------------
    if use_rayleigh_fading
        fprintf('步骤1: 应用Jakes瑞利衰落信道...\n');
        % 调用内置的jakesChannel函数（文件末尾定义）
        % 注意：jakesChannel返回comm.RayleighChannel对象，参数已内置
        % 参数：fs=3.84e6 Hz（与主程序一致），fd=100 Hz, 六径多径配置（频率选择性衰落）
        rchan = jakesChannel();
        
        % 应用瑞利衰落信道：使用comm.RayleighChannel对象
        % 注意：comm.RayleighChannel需要列向量输入，返回列向量
        Tx_data_col = Tx_data(:);  % 转换为列向量
        Rx_faded_col = rchan(Tx_data_col);  % 通过信道传输
        Rx_faded = Rx_faded_col(:)';  % 转换回行向量
        
        % 计算衰落后的信号功率（用于SNR计算）
        Tx_signal_power = var(Rx_faded);  % 使用衰落后的信号功率计算噪声
        
        % 重置信道状态（为下次使用准备）
        rchan.reset();
        
        fprintf('  - 信道类型: Jakes瑞利衰落 + AWGN\n');
        fprintf('  - 最大多普勒频移: %.1f Hz\n', fd);
        fprintf('  - 散射体数量: %d\n', N0);
        fprintf('  完成！\n\n');
    else
        % 不使用衰落，直接使用原始信号功率
        fprintf('步骤1: 跳过瑞利衰落（仅使用AWGN）...\n');
        Rx_faded = Tx_data;
        Tx_signal_power = var(Tx_data);
        fprintf('  完成！\n\n');
    end
    
    %------------------------------------------------------------------------------
    % 步骤2: 添加AWGN噪声
    %------------------------------------------------------------------------------
    fprintf('步骤2: 添加AWGN噪声...\n');
    % 噪声功率计算：根据目标SNR和信号功率计算噪声方差
    linear_SNR = 10^(targetSNRdB/10);
    noise_sigma = Tx_signal_power / linear_SNR;
    noise_scale_factor = sqrt(noise_sigma);
    noise = randn(1, length(Rx_faded)) * noise_scale_factor;
    Rx_data = Rx_faded + noise;
    
    fprintf('  - 目标SNR: %.1f dB\n', targetSNRdB);
    fprintf('  - 信号功率: %.6g\n', Tx_signal_power);
    fprintf('  - 噪声方差: %.6g\n', noise_sigma);
    fprintf('  完成！\n\n');
    
    %------------------------------------------------------------------------------
    % 信道传输结果总结
    %------------------------------------------------------------------------------
    fprintf('========================================\n');
    fprintf('信道传输完成！\n');
    fprintf('========================================\n');
    fprintf('接收信号参数：\n');
    fprintf('  - 接收信号长度: %d 样本\n', length(Rx_data));
    fprintf('  - 发送信号功率: %.6g\n', Tx_signal_power);
    fprintf('  - 噪声方差: %.6g\n', noise_sigma);
    if use_rayleigh_fading
        fprintf('  - 信道类型: Jakes瑞利衰落 + AWGN\n');
    else
        fprintf('  - 信道类型: AWGN only\n');
    end
    fprintf('  - 目标SNR: %.1f dB\n', targetSNRdB);
    fprintf('========================================\n\n');
    
    fprintf('发送信号已保存在变量 Tx_data 中\n');
    fprintf('接收信号已保存在变量 Rx_data 中\n');
    
    %------------------------------------------------------------------------------
    % 绘制接收信号的时域和频域图
    %------------------------------------------------------------------------------
    if plot_signal
        fprintf('\n正在绘制接收信号波形和频谱...\n');
        
        % 绘制接收信号时域和频域图
        figure('Name', 'OFDM接收信号时域和频域', 'Position', [150, 150, 1200, 600]);
        
        % 时域信号（前1000个样本）
        subplot(2, 1, 1);
        plot_samples_rx = min(1000, length(Rx_data));
        plot(1:plot_samples_rx, Rx_data(1:plot_samples_rx), 'r-', 'LineWidth', 1);
        grid on;
        xlabel('样本索引');
        ylabel('幅度');
        title(sprintf('接收信号时域波形（前%d个样本，SNR=%.1f dB）', plot_samples_rx, targetSNRdB));
        
        % 频域信号
        subplot(2, 1, 2);
        Nfft_rx = 2^nextpow2(length(Rx_data));
        Rx_Fz = fftshift(fft(Rx_data, Nfft_rx));
        f_axis_rx = (-Nfft_rx/2:(Nfft_rx/2-1)) * fs / Nfft_rx;
        
        plot(f_axis_rx/1e6, 20*log10(abs(Rx_Fz) / max(abs(Rx_Fz)) + eps), 'r-', 'LineWidth', 1);
        grid on;
        xlabel('频率 (MHz)');
        ylabel('幅度 (dB)');
        title(sprintf('接收信号频谱（归一化，SNR=%.1f dB）', targetSNRdB));
        xlim([-fs/2/1e6, fs/2/1e6]);
        
        fprintf('接收信号图形绘制完成！\n\n');
        
        % 可选：绘制发送和接收信号对比图
        fprintf('正在绘制发送与接收信号对比图...\n');
        figure('Name', '发送与接收信号对比', 'Position', [200, 200, 1400, 800]);
        
        % 时域对比（前1000个样本）
        subplot(2, 2, 1);
        plot_samples_comp = min(1000, min(length(Tx_data), length(Rx_data)));
        plot(1:plot_samples_comp, Tx_data(1:plot_samples_comp), 'b-', 'LineWidth', 1);
        grid on;
        xlabel('样本索引');
        ylabel('幅度');
        title(sprintf('发送信号时域（前%d个样本）', plot_samples_comp));
        legend('发送信号', 'Location', 'best');
        
        subplot(2, 2, 2);
        plot(1:plot_samples_comp, Rx_data(1:plot_samples_comp), 'r-', 'LineWidth', 1);
        grid on;
        xlabel('样本索引');
        ylabel('幅度');
        title(sprintf('接收信号时域（前%d个样本）', plot_samples_comp));
        legend('接收信号', 'Location', 'best');
        
        % 频域对比
        subplot(2, 2, 3);
        Nfft_tx = 2^nextpow2(length(Tx_data));
        Tx_Fz_comp = fftshift(fft(Tx_data, Nfft_tx));
        f_axis_tx = (-Nfft_tx/2:(Nfft_tx/2-1)) * fs / Nfft_tx;
        plot(f_axis_tx/1e6, 20*log10(abs(Tx_Fz_comp) / max(abs(Tx_Fz_comp)) + eps), 'b-', 'LineWidth', 1);
        grid on;
        xlabel('频率 (MHz)');
        ylabel('幅度 (dB)');
        title('发送信号频谱（归一化）');
        xlim([-fs/2/1e6, fs/2/1e6]);
        legend('发送信号', 'Location', 'best');
        
        subplot(2, 2, 4);
        plot(f_axis_rx/1e6, 20*log10(abs(Rx_Fz) / max(abs(Rx_Fz)) + eps), 'r-', 'LineWidth', 1);
        grid on;
        xlabel('频率 (MHz)');
        ylabel('幅度 (dB)');
        title(sprintf('接收信号频谱（归一化，SNR=%.1f dB）', targetSNRdB));
        xlim([-fs/2/1e6, fs/2/1e6]);
        legend('接收信号', 'Location', 'best');
        
        fprintf('对比图绘制完成！\n\n');
    end
else
    fprintf('\n注意：未启用信道传输，仅生成发送信号\n');
    fprintf('如需通过信道传输，请设置 use_channel = true\n\n');
    
    %------------------------------------------------------------------------------
    % 保存发送信号到mat文件（未启用信道传输时）
    %------------------------------------------------------------------------------
    fprintf('正在保存发送信号到mat文件...\n');
    
    % 生成文件名（包含参数信息）
    mat_filename = sprintf('OFDM_Tx_signal_%dsubcarriers_%s.mat', ...
        carrier_count, datestr(now, 'yyyymmdd_HHMMSS'));
    
    % 保存发送信号和相关参数
    save(mat_filename, 'Tx_data', 'fs', 'carrier_count', 'subcarrier_spacing', ...
        'IFFT_bin_length', 'GI', 'GIP', 'total_symbols', 'symbols_per_frame', ...
        '-v7.3');
    
    fprintf('发送信号已保存到文件: %s\n', mat_filename);
    fprintf('  包含变量: Tx_data, fs, carrier_count, ...\n');
    fprintf('  文件格式: MATLAB v7.3 (支持大文件)\n\n');
end

%===============================================================================
% 【带宽估计】（可选）
%===============================================================================
% 如果需要估计带宽，可以调用estimateBandwidth函数
% 示例：
%   [B_welch, B_ar, B_ideal, results] = estimateBandwidth(Tx_data, fs, ...
%       carrier_count, subcarrier_spacing, 'fc', 0, 'snr', 20, 'plot', true);
%
% 参数说明：
%   Tx_data          - 发送信号（已生成）
%   fs               - 采样频率
%   carrier_count    - 子载波数
%   subcarrier_spacing - 子载波间隔
%   'fc'             - 载波频率（可选，默认0）
%   'snr'            - 信噪比（可选，默认20dB）
%   'plot'           - 是否绘图（可选，默认true）
%
% 输出：
%   B_welch          - Welch算法估计的带宽
%   B_ar             - AR模型法估计的带宽
%   B_ideal          - 理论带宽
%   results          - 详细的估计结果结构体
%
% 注意：estimateBandwidth函数依赖于PSD_OFDM.m及相关函数
%===============================================================================

%===============================================================================
% 【内置依赖函数】
%===============================================================================
% 以下函数已内置，确保signalGenerate.m可以完全独立运行
%===============================================================================

function [complex_qam_data] = qam16(bitdata)
% ============================================================================
% 内置函数: qam16 - 16QAM调制函数
% 功能: 将输入的比特流调制为16QAM复信号
% 输入参数:
%   bitdata - 输入比特流（行向量，长度为4的倍数）
% 输出参数:
%   complex_qam_data - 16QAM调制后的复数符号序列
% ============================================================================

% 步骤1: 将比特流按4位一组进行分组（16QAM需要4比特表示一个符号）
x1 = reshape(bitdata, 4, length(bitdata)/4)';
d = 1;  % 最小距离（符号间最小欧氏距离）

% 步骤2: 将4位二进制数转换为十进制索引（0-15）
for i = 1:length(bitdata)/4
    for j = 1:4
        x1(i,j) = x1(i,j) * (2^(4-j)); % 16进制数对位转化成10进制
    end
    source(i,1) = 1 + sum(x1(i,:)); % 将16进制数转化为十进制后，再加一对位（MATLAB索引从1开始）
end

% 步骤3: 16QAM星座映射表（4×4格点，I/Q坐标）
mapping = [-3*d  3*d;
           -d    3*d;
            d    3*d;
           3*d   3*d;
           -3*d   d;
            -d    d;
             d    d;
            3*d   d;
           -3*d  -d; 
            -d   -d; 
             d   -d;
            3*d  -d;
           -3*d -3*d;
            -d  -3*d;
             d  -3*d;
            3*d -3*d];

% 步骤4: 根据索引映射到对应的星座点
for i = 1:length(bitdata)/4
    qam_data(i,:) = mapping(source(i),:); % 根据索引映射到对应的星座点
end

% 步骤5: 转换为复数表示（I+jQ）
complex_qam_data = complex(qam_data(:,1), qam_data(:,2));

end

function [rcosw] = rcoswindow(beta, Ts)
% ============================================================================
% 内置函数: rcoswindow - 升余弦窗函数生成器
% 功能: 生成升余弦（Raised Cosine）窗函数，用于OFDM符号的平滑加窗
% 输入参数:
%   beta - 滚降系数（过渡带占比，通常为1/16或1/32等）
%   Ts   - 窗函数覆盖长度（通常为包含循环前缀的OFDM符号长度 N+GI）
% 输出参数:
%   rcosw - 升余弦窗向量（列向量），长度为 round((1+beta)*Ts)+1
% ============================================================================

% 窗函数总长度：round((1+beta)*Ts)+1，分为三段

% 步骤1: 初始化参数（确保所有长度都是整数）
Ts_int = round(Ts);  % 确保 Ts 为整数
beta_Ts_int = round(beta*Ts_int);  % 确保 beta*Ts 为整数
win_length = round((1+beta)*Ts_int);  % 确保总长度为整数

% 步骤2: 初始化时间轴和窗函数数组
t = 0:win_length;  % 时间轴，长度为 win_length+1（索引 1 到 win_length+1）
rcosw = zeros(1, win_length+1);  % 窗函数数组，长度为 win_length+1（索引 1 到 win_length+1）

% 步骤3: 左过渡段 [0, beta*Ts]：从0平滑上升到1
% 公式：w(t) = 0.5 + 0.5*cos(π + t*π/(beta*Ts))
for i = 1:beta_Ts_int
    if i <= length(t) && beta_Ts_int > 0
        rcosw(i) = 0.5 + 0.5*cos(pi + t(i)*pi/(beta*Ts_int));
    end
end

% 步骤4: 平坦段 [beta*Ts, Ts]：保持为1（不衰减）
for i = (beta_Ts_int+1):Ts_int
    if i <= length(rcosw)
        rcosw(i) = 1;
    end
end

% 步骤5: 右过渡段 [Ts, (1+beta)*Ts]：从1平滑下降到0
% 公式：w(t) = 0.5 + 0.5*cos((t-Ts)*π/(beta*Ts))
for i = (Ts_int+1):(win_length+1)
    if i <= length(t) && beta_Ts_int > 0
        rcosw(i) = 0.5 + 0.5*cos((t(i)-Ts_int)*pi/(beta*Ts_int));
    end
end

% 步骤6: 转换为列向量（便于矩阵运算）
rcosw = rcosw';

end

function rchan = jakesChannel()
% ============================================================================
% 内置函数: jakesChannel - Jakes模型瑞利衰落信道创建函数
% 功能：创建并返回一个配置好的瑞利衰落信道对象
% 输出参数：
%   rchan - comm.RayleighChannel对象，配置好的瑞利衰落信道
% 信道参数（六径多径配置）：
%   - 采样率：3.84 MHz（与主程序一致，LTE 3MHz标准）
%   - 路径延迟：[0, 0.1e-6, 0.3e-6, 0.6e-6, 1.2e-6, 2.5e-6] 秒（六径）
%     对应时间：[0, 0.1, 0.3, 0.6, 1.2, 2.5] 微秒
%   - 路径功率：[0, -8, -12, -16, -20, -25] dB（主径功率最高，多径功率递减）
%     功率比：主径100%，多径依次为15.8%, 6.3%, 2.5%, 1.0%, 0.32%
%   - 最大多普勒频移：100 Hz
% ============================================================================

% 信道参数设置（六径多径配置，中等严重程度）
fs = 3.84e6;             % Sample rate (Hz) - 与主程序一致，LTE 3MHz标准
pathDelays = [0, 0.1e-6, 0.3e-6, 0.6e-6, 1.2e-6, 2.5e-6];  % Path delays (s) - 六径：0秒（主径），0.1, 0.3, 0.6, 1.2, 2.5微秒（多径）
pathPower = [0, -8, -12, -16, -20, -25];  % Path power (dB) - 六径功率：0 dB（主径），-8, -12, -16, -20, -25 dB（多径）
fD = 100;                % Maximum Doppler shift (Hz)

% 创建瑞利信道对象
rchan = comm.RayleighChannel('SampleRate', fs, ...
    'PathDelays', pathDelays, ...
    'AveragePathGains', pathPower, ...
    'MaximumDopplerShift', fD);

end
