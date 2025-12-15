function [B_welch, B_ar, B_ideal, results] = estimateBandwidth(Tx_data, fs, carrier_count, subcarrier_spacing, varargin)
%===============================================================================
% estimateBandwidth.m - OFDM信号带宽估计函数（使用signalGenerate.m的基带信号）
% 
% 功能说明:
%   使用Welch算法和AR模型法估计OFDM信号的带宽
%   - 输入信号：signalGenerate.m生成的基带信号（Tx_data）
%   - 信号处理：对基带信号进行上变频和加噪声处理
%   - 带宽估计：使用处理后的信号进行功率谱估计和带宽计算
%   - 可独立运行，自动调用signalGenerate.m生成信号
%
% 输入参数:
%   Tx_data          - OFDM发送信号（基带复数信号，行向量，来自signalGenerate.m）
%                      如果未提供，将从工作空间获取
%                      如果工作空间也没有，将自动调用signalGenerate.m生成信号
%   fs               - 采样频率（Hz）
%                      如果未提供，将从工作空间获取（默认15.36e6）
%   carrier_count    - 有效数据子载波数
%                      如果未提供，将从工作空间获取（默认300）
%   subcarrier_spacing - 子载波间隔（Hz）
%                      如果未提供，将从工作空间获取（默认15e3）
%   varargin         - 可选参数（名称-值对）:
%                      'fc'        - 载波频率（Hz），默认0（基带信号）
%                      'snr'       - 信噪比（dB），默认20
%                      'plot'      - 是否绘制PSD图，默认true
%
% 输出参数:
%   B_welch          - Welch算法估计的带宽值（Hz）
%   B_ar             - AR模型法估计的带宽值（Hz）
%   B_ideal          - 理论带宽值（Hz）
%   results          - 结构体，包含详细的估计结果和误差信息
%
% 信号处理流程:
%   1. 获取Tx_data（基带复数信号，来自signalGenerate.m）
%   2. 上变频：real(Tx_data * exp(j*2*pi*fc/fs*t))
%   3. 加噪声：awgn(sig_processed, snr, 'measured')
%   4. 使用处理后的信号进行功率谱估计和带宽计算
%
% 使用示例:
%   1. 无参数调用（自动调用signalGenerate.m生成Tx_data）:
%      [B_welch, B_ar, B_ideal, results] = estimateBandwidth;
%
%   2. 完整参数调用:
%      [B_welch, B_ar, B_ideal, results] = estimateBandwidth(Tx_data, fs, ...
%          carrier_count, subcarrier_spacing, 'fc', 0, 'snr', 20, 'plot', true);
%
% 创建日期: 2017.03.01 (原始版本)
% 修改日期: 2025.12.10 - 改造成函数，接受signalGenerate.m的基带信号
%===============================================================================

% 处理无参数调用情况：从工作空间获取参数或自动调用signalGenerate
if nargin == 0
    fprintf('\n========================================\n');
    fprintf('检测到无参数调用\n');
    fprintf('========================================\n');
    fprintf('正在从工作空间获取参数...\n');
    
    % 尝试从基础工作空间获取Tx_data
    tx_data_found = false;
    try
        Tx_data = evalin('base', 'Tx_data');
        if ~isempty(Tx_data)
            tx_data_found = true;
            fprintf('  ✓ 从基础工作空间获取 Tx_data\n');
        end
    catch
        % Tx_data不存在
    end
    
    % 如果基础工作空间没有，尝试调用者工作空间
    if ~tx_data_found
        try
            Tx_data = evalin('caller', 'Tx_data');
            if ~isempty(Tx_data)
                tx_data_found = true;
                fprintf('  ✓ 从调用者工作空间获取 Tx_data\n');
            end
        catch
            % Tx_data不存在
        end
    end
    
    % 如果仍然没有找到Tx_data，自动调用signalGenerate.m生成信号
    if ~tx_data_found
        fprintf('  ✗ 工作空间中未找到 Tx_data\n');
        fprintf('  正在自动调用 signalGenerate.m 生成信号...\n');
        
        % 检查signalGenerate.m是否存在
        if exist('signalGenerate.m', 'file') == 2
            fprintf('  ✓ 找到 signalGenerate.m 文件\n');
            fprintf('----------------------------------------\n');
            
            % 运行signalGenerate.m（在基础工作空间中）
            try
                evalin('base', 'signalGenerate');
                pause(0.1);
                
                % 从基础工作空间获取生成的发送信号（Tx_data）
                Tx_data = evalin('base', 'Tx_data');
                if isempty(Tx_data)
                    error('Tx_data为空');
                end
                
                % 同时获取其他参数
                try
                    fs = evalin('base', 'fs');
                catch
                    fprintf('  ⚠ 未找到fs变量，使用默认值15.36e6 Hz\n');
                    fs = 15.36e6;
                end
                
                try
                    carrier_count = evalin('base', 'carrier_count');
                catch
                    fprintf('  ⚠ 未找到carrier_count变量，使用默认值300\n');
                    carrier_count = 300;
                end
                
                try
                    subcarrier_spacing = evalin('base', 'subcarrier_spacing');
                catch
                    fprintf('  ⚠ 未找到subcarrier_spacing变量，使用默认值15e3 Hz\n');
                    subcarrier_spacing = 15e3;
                end
                
                fprintf('----------------------------------------\n');
                fprintf('  ✓ 信号生成完成\n');
                fprintf('    Tx_data长度: %d 样本\n', length(Tx_data));
                fprintf('    fs: %.2e Hz\n', fs);
                fprintf('    carrier_count: %d\n', carrier_count);
                fprintf('    subcarrier_spacing: %.1f Hz\n', subcarrier_spacing);
                fprintf('========================================\n\n');
                
                skip_param_fetch = true;
            catch ME
                error('错误：自动调用signalGenerate.m失败。\n错误信息：%s\n\n请手动运行signalGenerate.m生成信号，或提供Tx_data参数。', ME.message);
            end
        else
            error('错误：工作空间中未找到Tx_data变量，且signalGenerate.m文件不存在。\n请先运行signalGenerate.m生成信号，或提供Tx_data参数。');
        end
    end
    
    % 如果自动调用了signalGenerate，参数已经获取，跳过后续步骤
    if ~exist('skip_param_fetch', 'var') || ~skip_param_fetch
        % 尝试从工作空间获取fs
        try
            fs = evalin('base', 'fs');
        catch
            try
                fs = evalin('caller', 'fs');
            catch
                fprintf('  ⚠ 未找到fs变量，使用默认值15.36e6 Hz\n');
                fs = 15.36e6;
            end
        end
        
        % 尝试从工作空间获取carrier_count
        try
            carrier_count = evalin('base', 'carrier_count');
        catch
            try
                carrier_count = evalin('caller', 'carrier_count');
            catch
                fprintf('  ⚠ 未找到carrier_count变量，使用默认值300\n');
                carrier_count = 300;
            end
        end
        
        % 尝试从工作空间获取subcarrier_spacing
        try
            subcarrier_spacing = evalin('base', 'subcarrier_spacing');
        catch
            try
                subcarrier_spacing = evalin('caller', 'subcarrier_spacing');
            catch
                fprintf('  ⚠ 未找到subcarrier_spacing变量，使用默认值15e3 Hz\n');
                subcarrier_spacing = 15e3;
            end
        end
        
        fprintf('  ✓ 参数获取完成\n');
        fprintf('    Tx_data长度: %d 样本\n', length(Tx_data));
        fprintf('    fs: %.2e Hz\n', fs);
        fprintf('    carrier_count: %d\n', carrier_count);
        fprintf('    subcarrier_spacing: %.1f Hz\n', subcarrier_spacing);
        fprintf('========================================\n\n');
    end
    
    % 解析可选参数（无参数调用时，varargin为空）
    varargin = {};
end

% 记录仿真开始时间
sim_start_time = tic;
start_time_str = datestr(now, 'yyyy-mm-dd HH:MM:SS');

% 解析可选参数
p = inputParser;
addParameter(p, 'fc', 0, @isnumeric);           % 载波频率，默认0
addParameter(p, 'snr', 20, @isnumeric);         % 信噪比，默认20dB
addParameter(p, 'plot', true, @islogical);      % 是否绘图，默认true
parse(p, varargin{:});

fc = p.Results.fc;
snr_estimate = p.Results.snr;
plot_flag = p.Results.plot;

% 保存原始采样频率（Tx_data的采样频率，来自signalGenerate.m）
fs_original = fs;  % 保存Tx_data的原始采样频率

% 计算理论带宽
B_ideal = carrier_count * subcarrier_spacing;

% 输出仿真参数
fprintf('\n========================================\n');
fprintf('OFDM信号带宽估计仿真\n');
fprintf('========================================\n');
fprintf('开始时间: %s\n', start_time_str);
fprintf('\n【输入参数】\n');
fprintf('  基带信号来源: signalGenerate.m 生成的 Tx_data\n');
fprintf('  信号长度: %d 样本（Tx_data，基带信号）\n', length(Tx_data));
fprintf('  采样频率: %.2e Hz (%.2f MHz) [输入参数，Tx_data的采样频率]\n', fs_original, fs_original/1e6);
fprintf('  有效数据子载波数: %d\n', carrier_count);
fprintf('  子载波间隔: %.1f Hz (%.3f kHz)\n', subcarrier_spacing, subcarrier_spacing/1e3);
fprintf('  载波频率: %.6e Hz (%.3f MHz) [输入参数，将被覆盖]\n', fc, fc/1e6);
fprintf('  信噪比: %.1f dB [输入参数，用于单次估计]\n', snr_estimate);
fprintf('  是否绘图: %s\n', mat2str(plot_flag));
fprintf('\n【理论参数】\n');
fprintf('  理论带宽: %.6e Hz (%.3f MHz)\n', B_ideal, B_ideal/1e6);
fprintf('========================================\n\n');

%**************************************************************************
% 信号预处理：检查Tx_data类型并转换为统一格式
%**************************************************************************
fprintf('\n【步骤1】信号预处理\n');
fprintf('----------------------------------------\n');
fprintf('使用从signalGenerate.m引用的基带信号 Tx_data...\n');
fprintf('  信号来源: signalGenerate.m 生成的 Tx_data（基带OFDM发送信号）\n');

% 检查Tx_data是否为复数信号
if isreal(Tx_data)
    fprintf('  ✓ 检测到Tx_data为实数信号（来自signalGenerate.m的加窗后实信号）\n');
    fprintf('    说明: signalGenerate.m使用共轭对称映射，IFFT输出为实信号\n');
    sig_baseband = Tx_data + 0*1j;  % 转换为复数形式
else
    fprintf('  ✓ 检测到Tx_data为复数信号\n');
    sig_baseband = Tx_data;
end

% 计算信号长度信息（使用原始采样频率）
sig_length = length(sig_baseband);
sig_duration = sig_length / fs_original;  % 使用原始采样频率
sig_duration_us = sig_duration * 1e6;
sig_duration_ms = sig_duration * 1e3;

fprintf('  信号长度: %d 个采样点\n', sig_length);
fprintf('  采样频率: %.2e Hz (%.2f MHz) [Tx_data的原始采样频率]\n', fs_original, fs_original/1e6);
fprintf('  信号时长: %.4f 秒 = %.2f 毫秒 = %.2f 微秒\n', ...
    sig_duration, sig_duration_ms, sig_duration_us);
fprintf('----------------------------------------\n\n');

%**************************************************************************
% 设置仿真参数（用于Rayleigh信道仿真和带宽检测率计算）
%**************************************************************************
fprintf('【步骤2】设置仿真参数\n');
fprintf('----------------------------------------\n');
N = 20;          % 构成一个帧结构的OFDM信号的个数（用于PSD_generate，已弃用）
snr = -4:2:10;   % 信噪比范围（用于带宽检测率计算）
fc = 10e6;       % 载波频率（用于上变频和Rayleigh信道仿真）
fs_sim = 40e6;   % 仿真采样频率（用于Rayleigh信道仿真，覆盖输入的fs）
itau = floor([0,1e-8,2e-8,5e-8,2e-7,5e-7].*fs_sim);  % 多径延时（样本数）
power = [0,-1.0,-7.0,-10.0,-12.0,-17.0];  % 多径信道的每径功率（dB）
fmax = 20;       % 最大多普勒频率（Hz）
itn = [10000,20000,30000,40000,50000,60000];  % 瑞利信道的记录次数

% 注意：fs_sim会覆盖输入的fs参数，用于Rayleigh信道仿真
fprintf('  载波频率 fc = %.2e Hz (%.2f MHz)\n', fc, fc/1e6);
fprintf('  仿真采样频率 fs = %.2e Hz (%.2f MHz) [用于Rayleigh信道]\n', fs_sim, fs_sim/1e6);
fprintf('  信噪比范围 SNR = [%d:%d:%d] dB\n', snr(1), snr(2)-snr(1), snr(end));
fprintf('  最大多普勒频率 fmax = %d Hz\n', fmax);
fprintf('  多径数 = %d 径\n', length(itau));
fprintf('  多径延时 itau = [%s] 样本\n', num2str(itau));
fprintf('  多径功率 power = [%s] dB\n', num2str(power));
fprintf('----------------------------------------\n\n');

% 更新fs为仿真采样频率（用于后续处理）
fs = fs_sim;
% 
% %**************************************************************************
% [c42_ofdm,c42_qpsk,c42_qam16,c42_qam64,c42_fsk8] = cumulant(snr1,N,para,ratio,K);
% [VAR_FSK8,VAR_QAM16,VAR_QAM64,VAR_QPSK,VAR_OFDM] = sc_ofdm_wavelet(N,snr1,para,ratio,K);
% [VAR_FSK8_Ray,VAR_QAM16_Ray,VAR_QAM64_Ray,VAR_QPSK_Ray,VAR_OFDM_Ray] = sc_ofdm_wavelet_Ray(N,snr1,para,ratio,K,itau,power,itn,fmax,fs);
% %**************************************************************************
% %�߽������ķ�ʽ��ͬ��������źŵ�ʶ����
% %rate_ofdm��ʾ�źŵ���ȷʶ����
% %**************************************************************************
% [rate_ofdm] = OFDM_rate(snr,N,para,ratio);
% [rate_ofdm1] = OFDM_rate_xiaobo(snr,N,para,ratio);
% [rate_ofdm_Ray] = OFDM_rate_xiaobo_Ray(snr,N,para,ratio,itau,power,itn,fmax,fs);
% %**************************************************************************
% %OFDM不同子载波数目对高阶累积量估计值的不确定性影响
% %c42_ofdm64,c42_ofdm128,c42_ofdm256分别表示
% %子载波数目为128,256,512时的累积量
% %**************************************************************************
% [V_ofdm128,V_ofdm256,V_ofdm512] = ofdm_dif_para(snr,N,ratio,K);
%  
% %**************************************************************************
% %ofdm�ز�Ƶ�ʲ�������
% %ofdmѭ���׹���ofdm�ź��ز�Ƶ��
% %**************************************************************************
% trst_rate = 5e5;  % �źŷ�������,����Ƭ����
% M = 16; % ƽ������ 
% sig = ofdm(N,para,ratio);
% Ns = length(sig);       % ѭ���׼��������ȣ�����С�ڵ����ź����г���
%  
% %**************************************************************************
% %��������OFDM�źŲ���ѭ��ƽ����
% %**************************************************************************
% s_n = ceil(fs/trst_rate); % �������ʽ���ΪOFDM�ź����ʵ������� 
% sign = sig(ones(s_n,1),:); % ��Ԫ����
% sign = reshape(sign, 1, s_n*length(sign));
% data = sign.*exp(1i*2*pi*fc/fs*(0:length(sign)-1));
% %�߰��ŵ���ѭ���׹�����Ƶ
% data_awgn = awgn(data,25,'measured');
% [f_awgn] = cyclic_spectrum(real(data_awgn), Ns, fs, M,K);  % ѭ������Ƶ��f
% 
% %�����ŵ���ѭ���׹���
% %data_rayleigh = RayFade(itau,power,fmax,fs,data,25);
% data_rayleigh = (MUL_RAYLEIGH(data,itau,power,itn,length(itau),length(data),1/fs,fmax,0));
% data_rayleigh = data_rayleigh/std(data_rayleigh);
% data_rayleigh = awgn(data_rayleigh,25,'measured');
% [f_rayleigh] = cyclic_spectrum(real(data_rayleigh), Ns, fs, M,K);  % ѭ������Ƶ��f
%**************************************************************************
% 绘制输入基带信号 Tx_data（时域和频域）
%**************************************************************************
if plot_flag
    fprintf('【步骤3】绘制输入基带信号 Tx_data\n');
    fprintf('----------------------------------------\n');
    fprintf('  注意: Tx_data来自signalGenerate.m，采样频率为 %.2f MHz\n', fs_original/1e6);
    fprintf('        绘制Tx_data的频域图时使用原始采样频率，而不是后续仿真用的 %.2f MHz\n', fs/1e6);
    fprintf('        这样可以正确显示Tx_data的实际频率范围和带宽\n');
    
    % 绘制基带信号的时域和频域图
    figure('Name', '从signalGenerate.m引用的基带信号 Tx_data', 'Position', [100, 100, 1400, 800]);
    
    % 时域信号 - 实部（或实数信号）
    % 注意：使用Tx_data的原始采样频率fs_original计算时间轴
    subplot(2, 2, 1);
    time_axis = (0:sig_length-1) / fs_original * 1e6;  % 使用原始采样频率，转换为微秒
    sig_duration_original = sig_length / fs_original;  % 使用原始采样频率计算时长
    sig_duration_us_original = sig_duration_original * 1e6;
    % 如果信号太长，只显示一部分以便观察细节
    display_samples = min(5000, sig_length);  % 最多显示5000个采样点
    plot(time_axis(1:display_samples), real(sig_baseband(1:display_samples)), 'b-', 'LineWidth', 0.5);
    xlabel('时间 (μs)');
    ylabel('幅度');
    title(sprintf('Tx_data 时域 - 实部\n总长度: %d 个采样点, 时长: %.2f μs (fs=%.2f MHz)\n显示: 前 %d 个采样点', ...
        sig_length, sig_duration_us_original, fs_original/1e6, display_samples));
    grid on;
    
    % 时域信号 - 虚部（如果Tx_data是复数）
    subplot(2, 2, 2);
    if ~isreal(sig_baseband) && any(imag(sig_baseband) ~= 0)
        plot(time_axis(1:display_samples), imag(sig_baseband(1:display_samples)), 'r-', 'LineWidth', 0.5);
        title(sprintf('Tx_data 时域 - 虚部\n总长度: %d 个采样点, 时长: %.2f μs (fs=%.2f MHz)\n显示: 前 %d 个采样点', ...
            sig_length, sig_duration_us_original, fs_original/1e6, display_samples));
    else
        plot(time_axis(1:display_samples), zeros(1, display_samples), 'r-', 'LineWidth', 0.5);
        title(sprintf('Tx_data 时域 - 虚部（实数信号，虚部为0）\n总长度: %d 个采样点, 时长: %.2f μs (fs=%.2f MHz)', ...
            sig_length, sig_duration_us_original, fs_original/1e6));
    end
    xlabel('时间 (μs)');
    ylabel('幅度');
    grid on;
    
    % 频域信号 - 使用FFT（双边频谱）
    % 注意：使用Tx_data的原始采样频率fs_original，而不是被覆盖后的fs
    subplot(2, 2, 3);
    N_fft = 2^nextpow2(sig_length);
    freq_axis_neg = (-N_fft/2:N_fft/2-1) * fs_original / N_fft / 1e6;  % 使用原始采样频率，转换为MHz，包含负频率
    fft_sig = fft(sig_baseband, N_fft);
    fft_sig_shifted = fftshift(fft_sig);  % 零频率居中
    fft_mag = 20*log10(abs(fft_sig_shifted) + eps);  % 加eps避免log(0)
    fft_mag = fft_mag - max(fft_mag);  % 归一化到最大值
    plot(freq_axis_neg, fft_mag, 'b-', 'LineWidth', 0.5);
    xlabel('频率 (MHz)');
    ylabel('幅度 (dB)');
    title(sprintf('Tx_data 频域 - 双边频谱\n(使用原始采样频率 %.2f MHz)', fs_original/1e6));
    grid on;
    hold on;
    % 标记零频率
    plot([0, 0], [min(fft_mag), max(fft_mag)], 'k--', 'LineWidth', 1, ...
        'DisplayName', '零频率 (基带)');
    % 标记理论带宽范围（基带，中心在0）
    bandwidth_base_mhz = B_ideal / 1e6;
    plot([-bandwidth_base_mhz/2, -bandwidth_base_mhz/2], [min(fft_mag), max(fft_mag)], ...
        'g--', 'LineWidth', 1.5, 'DisplayName', sprintf('理论带宽边界 -%.2f MHz', bandwidth_base_mhz/2));
    plot([bandwidth_base_mhz/2, bandwidth_base_mhz/2], [min(fft_mag), max(fft_mag)], ...
        'g--', 'LineWidth', 1.5, 'DisplayName', sprintf('理论带宽边界 +%.2f MHz', bandwidth_base_mhz/2));
    legend('Location', 'best');
    xlim([-fs_original/2/1e6, fs_original/2/1e6]);
    
    % 频域信号 - 单边频谱（只显示正频率）
    subplot(2, 2, 4);
    freq_axis_pos = (0:N_fft/2-1) * fs_original / N_fft / 1e6;  % 使用原始采样频率，转换为MHz，只显示正频率
    fft_mag_pos = 20*log10(abs(fft_sig(1:N_fft/2)) + eps);
    fft_mag_pos = fft_mag_pos - max(fft_mag_pos);  % 归一化到最大值
    plot(freq_axis_pos, fft_mag_pos, 'b-', 'LineWidth', 0.5);
    xlabel('频率 (MHz)');
    ylabel('幅度 (dB)');
    title(sprintf('Tx_data 频域 - 单边频谱\n(使用原始采样频率 %.2f MHz)', fs_original/1e6));
    grid on;
    hold on;
    % 标记理论带宽上界（基带，从0开始）
    plot([bandwidth_base_mhz, bandwidth_base_mhz], [min(fft_mag_pos), max(fft_mag_pos)], ...
        'g--', 'LineWidth', 1.5, 'DisplayName', sprintf('理论带宽 %.2f MHz', bandwidth_base_mhz));
    legend('Location', 'best');
    xlim([0, fs_original/2/1e6]);
    
    fprintf('  ✓ 图形生成完成\n');
    fprintf('----------------------------------------\n\n');
end

%**************************************************************************
% 绘制经过处理后的基带信号（基于Tx_data，经过上变频、瑞利衰落和加噪声）
%**************************************************************************
if plot_flag
    fprintf('【步骤4】绘制处理后信号（Rayleigh信道）\n');
    fprintf('----------------------------------------\n');
    fprintf('  注意: Tx_data的原始采样频率: %.2f MHz\n', fs_original/1e6);
    fprintf('        Rayleigh信道仿真采样频率: %.2f MHz\n', fs/1e6);
    snr_plot_signal = 20;  % 用于绘图的SNR值
    
    % 处理信号：与PSD_OFDM_rayleigh函数使用相同的处理流程
    % 为了与PSD_OFDM_rayleigh函数保持一致，需要先将信号重采样到fs（40MHz）
    
    % 如果采样频率不同，需要重采样到仿真采样频率
    if abs(fs_original - fs) > 1  % 采样频率不同（容差1Hz）
        fprintf('  将sig_baseband重采样到 %.2f MHz（用于Rayleigh信道仿真）...\n', fs/1e6);
        % 计算重采样比例，使用最大公约数简化分数
        % resample(x, p, q) 将信号重采样到 p/q 倍的采样率
        % 需要找到 fs/fs_original 的最简分数形式
        ratio = fs / fs_original;
        % 使用目标样本数方式重采样（更简单可靠）
        target_length = round(length(sig_baseband) * ratio);
        sig_baseband_resampled_plot = resample(sig_baseband, target_length);
        fprintf('  重采样完成：%d -> %d 样本 (比例: %.4f)\n', ...
            length(sig_baseband), length(sig_baseband_resampled_plot), ratio);
        sig_baseband_for_plot = sig_baseband_resampled_plot;
        fs_for_plot = fs;  % 使用重采样后的采样频率
    else
        % 采样频率相同，直接使用
        sig_baseband_for_plot = sig_baseband;
        fs_for_plot = fs_original;
    end
    
    % 1. 上变频（将基带信号上变频到载波频率fc）
    % 使用重采样后的采样频率（与PSD_OFDM_rayleigh保持一致）
    sig_processed = (sig_baseband_for_plot.*exp(1j*2*pi*fc/fs_for_plot*(0:length(sig_baseband_for_plot)-1)));
    
    % 2. 通过瑞利衰落信道（使用仿真采样频率fs）
    sig_processed = real(MUL_RAYLEIGH(sig_processed,itau,power,itn,length(itau),length(sig_processed),1/fs_for_plot,fmax,0));
    
    % 3. 加AWGN噪声
    sig_processed = awgn(sig_processed, snr_plot_signal, 'measured');

    
    % 计算信号长度信息（使用重采样后的采样频率）
    signal_length = length(sig_processed);  % 采样点数
    signal_duration = signal_length / fs_for_plot;  % 信号时长（秒），使用重采样后的采样频率
    signal_duration_us = signal_duration * 1e6;  % 转换为微秒
    signal_duration_ms = signal_duration * 1e3;  % 转换为毫秒
    
    fprintf('处理后信号参数（基于Tx_data，经过Rayleigh衰落信道）：\n');
    fprintf('  采样点数: %d 个\n', signal_length);
    fprintf('  采样频率: %.2f MHz (重采样后，用于Rayleigh信道仿真)\n', fs_for_plot/1e6);
    fprintf('  信号时长: %.4f 秒 = %.2f 毫秒 = %.2f 微秒\n', ...
        signal_duration, signal_duration_ms, signal_duration_us);
    fprintf('  采样间隔: %.2f 纳秒\n', 1/fs_for_plot*1e9);
    fprintf('  信道类型: Rayleigh多径衰落 + AWGN\n');
    fprintf('  信噪比: %.1f dB\n', snr_plot_signal);
    fprintf('  多径数: %d 径\n', length(itau));
    fprintf('  最大多普勒频率: %.1f Hz\n\n', fmax);
    
    % 绘制时域和频谱图
    figure('Name', '基于Tx_data的处理后信号（Rayleigh信道）', 'Position', [100, 100, 1400, 600]);

    % 时域信号（完整信号）
    subplot(1, 2, 1);
    time_axis = (0:signal_length-1) / fs_for_plot * 1e6;  % 使用重采样后的采样频率，转换为微秒
    plot(time_axis, sig_processed);
    xlabel('时间 (μs)');
    ylabel('幅度');
    title(sprintf('基于Tx_data的处理后信号时域波形 (SNR = %.1f dB)\n总长度: %d 个采样点, 时长: %.2f μs (fs=%.2f MHz)', ...
        snr_plot_signal, signal_length, signal_duration_us, fs_for_plot/1e6));
    grid on;
    xlim([time_axis(1), time_axis(end)]);
    
    % 频谱信号（FFT）
    % 注意：使用fs_for_plot计算频率轴，与上变频和Rayleigh信道处理保持一致
    subplot(1, 2, 2);
    N_fft = 2^nextpow2(length(sig_processed));
    freq_axis = (0:N_fft/2-1) * fs_for_plot / N_fft / 1e6;  % 使用重采样后的采样频率，转换为MHz
    fft_sig = fft(sig_processed, N_fft);
    fft_mag = 20*log10(abs(fft_sig(1:N_fft/2)) + eps);
    fft_mag = fft_mag - max(fft_mag);  % 归一化到最大值
    plot(freq_axis, fft_mag);
    xlabel('频率 (MHz)');
    ylabel('幅度 (dB)');
    title(sprintf('基于Tx_data的处理后信号频谱 (SNR = %.1f dB)\n载波频率: %.1f MHz, 理论带宽: %.2f MHz (fs=%.2f MHz)', ...
        snr_plot_signal, fc/1e6, B_ideal/1e6, fs_for_plot/1e6));
    grid on;
    hold on;
    % 标记载波频率
    fc_mhz = fc/1e6;
    plot([fc_mhz, fc_mhz], [min(fft_mag), max(fft_mag)], 'r--', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('载波频率 %.1f MHz', fc_mhz));
    % 标记理论带宽范围（中心在fc）
    bandwidth_base_mhz_processed = B_ideal / 1e6;
    bandwidth_lower = (fc - bandwidth_base_mhz_processed*1e6/2)/1e6;
    bandwidth_upper = (fc + bandwidth_base_mhz_processed*1e6/2)/1e6;
    plot([bandwidth_lower, bandwidth_lower], [min(fft_mag), max(fft_mag)], 'g--', ...
        'LineWidth', 1, 'DisplayName', sprintf('带宽下边界 %.2f MHz', bandwidth_lower));
    plot([bandwidth_upper, bandwidth_upper], [min(fft_mag), max(fft_mag)], 'g--', ...
        'LineWidth', 1, 'DisplayName', sprintf('带宽上边界 %.2f MHz', bandwidth_upper));
    legend('Location', 'best');
    xlim([0, fs_for_plot/2/1e6]);  % 使用重采样后的采样频率的一半作为最大频率
    fprintf('  ✓ 图形生成完成\n');
    fprintf('----------------------------------------\n\n');
end

%**************************************************************************
% Rayleigh信道下的带宽检测率计算（基于signalGenerate.m的Tx_data）
% 注意：sig_baseband 来自 signalGenerate.m 的 Tx_data（基带信号）
%**************************************************************************
fprintf('【步骤5】计算Rayleigh信道带宽检测率\n');
fprintf('----------------------------------------\n');
fprintf('  使用基带信号: signalGenerate.m 的 Tx_data (sig_baseband)\n');
fprintf('  Tx_data采样频率: %.2f MHz\n', fs_original/1e6);
fprintf('  Rayleigh信道仿真采样频率: %.2f MHz\n', fs/1e6);

% 如果采样频率不同，需要重采样到仿真采样频率
if abs(fs_original - fs) > 1  % 采样频率不同（容差1Hz）
    fprintf('  检测到采样频率不同，将sig_baseband重采样到 %.2f MHz...\n', fs/1e6);
    % 重采样到仿真采样频率
    % 使用目标样本数方式重采样
    ratio = fs / fs_original;
    target_length = round(length(sig_baseband) * ratio);
    sig_baseband_resampled = resample(sig_baseband, target_length);
    fprintf('  重采样完成：%d -> %d 样本 (比例: %.4f)\n', ...
        length(sig_baseband), length(sig_baseband_resampled), ratio);
    sig_baseband_for_rayleigh = sig_baseband_resampled;
else
    % 采样频率相同，直接使用
    sig_baseband_for_rayleigh = sig_baseband;
end

rayleigh_start = tic;
[B_rate_welch_rayleigh,B_rate_ar_rayleigh ] = Bandwidth_rate_rayleigh(sig_baseband_for_rayleigh,fc,fs,snr,itau,power,fmax,itn,B_ideal);
rayleigh_time = toc(rayleigh_start);
fprintf('  ✓ 带宽检测率计算完成，耗时: %.2f 秒\n', rayleigh_time);
fprintf('----------------------------------------\n\n');

%**************************************************************************
% 生成PSD图形并计算单次带宽估计值
% 注意：sig_baseband 来自 signalGenerate.m 的 Tx_data（基带信号）
%**************************************************************************
fprintf('【步骤6】生成PSD图形并计算带宽估计值\n');
fprintf('----------------------------------------\n');
fprintf('  使用基带信号: signalGenerate.m 的 Tx_data (sig_baseband)\n');
fprintf('  Tx_data采样频率: %.2f MHz\n', fs_original/1e6);
fprintf('  Rayleigh信道仿真采样频率: %.2f MHz\n', fs/1e6);

% 如果采样频率不同，需要重采样到仿真采样频率
% 注意：如果之前已经重采样过，这里可以重用
if abs(fs_original - fs) > 1  % 采样频率不同（容差1Hz）
    if ~exist('sig_baseband_for_rayleigh', 'var')
        fprintf('  将sig_baseband重采样到 %.2f MHz...\n', fs/1e6);
        % 使用目标样本数方式重采样
        ratio = fs / fs_original;
        target_length = round(length(sig_baseband) * ratio);
        sig_baseband_resampled = resample(sig_baseband, target_length);
        fprintf('  重采样完成：%d -> %d 样本 (比例: %.4f)\n', ...
            length(sig_baseband), length(sig_baseband_resampled), ratio);
        sig_baseband_for_rayleigh = sig_baseband_resampled;
    else
        fprintf('  使用之前重采样后的信号\n');
    end
else
    % 采样频率相同，直接使用
    sig_baseband_for_rayleigh = sig_baseband;
end

psd_start = tic;
snr_plot = 20;  % 用于生成 PSD 图形的 SNR 值（标量）
[B_welch,B_ar] = PSD_OFDM_rayleigh(sig_baseband_for_rayleigh,fc,fs,snr_plot,plot_flag,itau,power,fmax,itn,B_ideal);
psd_time = toc(psd_start);
fprintf('  ✓ PSD估计完成，耗时: %.2f 秒\n', psd_time);
fprintf('  ✓ Welch算法估计带宽: %.2f MHz\n', B_welch/1e6);
fprintf('  ✓ AR模型估计带宽: %.2f MHz\n', B_ar/1e6);
fprintf('----------------------------------------\n\n');

%**************************************************************************
% 绘制带宽检测率对比图
%**************************************************************************
fprintf('【步骤7】绘制结果图形\n');
fprintf('----------------------------------------\n');
plot_start = tic;
figure('Name', 'Rayleigh信道下的带宽检测率对比')
plot(snr,B_rate_ar_rayleigh,'r-o', 'LineWidth', 1.5, 'MarkerSize', 8);
hold on
plot(snr,B_rate_welch_rayleigh,'b--s', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('SNR (dB)', 'FontSize', 12);
ylabel('带宽检测率 (%)', 'FontSize', 12);
legend('AR模型法','Welch算法', 'Location', 'best', 'FontSize', 11);
title('Rayleigh信道下的带宽检测率对比', 'FontSize', 13);
grid on;
plot_time = toc(plot_start);
fprintf('  ✓ 图形绘制完成，耗时: %.2f 秒\n', plot_time);
fprintf('----------------------------------------\n\n');

% 输出总仿真时间和结果总结
total_time = toc(sim_start_time);
fprintf('========================================\n');
fprintf('【仿真完成】\n');
fprintf('========================================\n');
fprintf('结束时间: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf('总仿真时间: %.2f 秒 (%.2f 分钟)\n', total_time, total_time/60);
fprintf('\n【结果总结】\n');
fprintf('  理论带宽: %.2f MHz\n', B_ideal/1e6);
fprintf('  Welch算法估计: %.2f MHz (SNR=%.1f dB时)\n', B_welch/1e6, snr_plot);
fprintf('  AR模型估计: %.2f MHz (SNR=%.1f dB时)\n', B_ar/1e6, snr_plot);
fprintf('  Welch相对误差: %.2f%%\n', abs(B_welch - B_ideal) / B_ideal * 100);
fprintf('  AR模型相对误差: %.2f%%\n', abs(B_ar - B_ideal) / B_ideal * 100);
fprintf('========================================\n');

%**************************************************************************
%AWGN信道下循环自相关估计有效数据长度和总长度以及循环前缀长度
%**************************************************************************
% sig_a = sig.*exp(1i*2*pi*fc/fs*(0:length(sig)-1));  %上变频到载频 
% [Tu_128,Ts_128,Tg_128 ] = effectivelength(sig_a,fs,20,N,K);
% [ Tu_rate_128,Ts_rate_128,Tg_rate_128] = Length_rate(sig_a,fs,snr,N,1 );
% 
% %**************************************************************************
% %Rayleigh信道下估计有效数据长度和总长度以及循环前缀长度
% %**************************************************************************
% [Tu_128_rayleigh,Ts_128_rayleigh,Tg_128_rayleigh ] = effectivelength_rayleigh(sig_a,fs,10,N,1,itau,power,fmax,itn);
% [Tu_rate_128_rayleigh,Ts_rate_128_rayleigh,Tg_rate_128_rayleigh] = Length_rate_rayleigh(sig_a,fs,snr,N,0,itau,power,fmax,itn );
% figure
% plot(snr,Tu_rate_128_rayleigh,'r-o');
% hold on

% plot(snr,Tu_rate_128,'k-x');
% xlabel('snr/db');
% ylabel('percentage/%');
% legend('Rayleigh','Awgn');
% title('不同信道对有效数据长度估计的影响')
% figure
% plot(snr,Ts_rate_128_rayleigh,'r-o');
% hold on
% plot(snr,Ts_rate_128,'k-x');
% legend('Rayleigh','Awgn');
% xlabel('snr/db');
% ylabel('percentage/%');
% title('不同信道对符号总长度估计的影响')
% figure
% plot(snr,Tg_rate_128_rayleigh,'r-o');
% hold on
% plot(snr,Tg_rate_128,'k-x');
% legend('Rayleigh','Awgn');
% xlabel('snr/db');
% ylabel('percentage/%');
% title('不同信道对循环前缀长度估计的影响')
% 
% %**************************************************************************
% %AWGN信道下子载波数目对信号参数估计性能的影响
% %**************************************************************************
% % sig_256 = ofdm(N,256,ratio);
% % sig256 = sig_256.*exp(1i*2*pi*fc/fs*(0:length(sig_256)-1));  %上变频到载频 
% % [ Tu_rate_256,Ts_rate_256,Tg_rate_256] = Length_rate_difpara(sig256,fs,snr,N,0 );
% % figure
% % plot(snr,Tu_rate_256,'r-o');
% % hold on
% % plot(snr,Tu_rate_128,'k-x');
% % xlabel('snr/db');
% % ylabel('percentage/%');
% % legend('256个子载波','128个子载波');
% % title('子载波数目对有效数据长度估计的影响')
% % figure
% % plot(snr,Ts_rate_256,'r-o');
% % hold on
% % plot(snr,Ts_rate_128,'k-x');
% % legend('256个子载波','128个子载波');
% % xlabel('snr/db');
% % ylabel('percentage/%');
% % title('子载波数目对符号总长度估计的影响')
% % figure
% % plot(snr,Tg_rate_256,'r-o');
% % hold on
% % plot(snr,Tg_rate_128,'k-x');
% % legend('256个子载波','128个子载波');
% % xlabel('snr/db');
% % ylabel('percentage/%');
% % title('子载波数目对循环前缀长度估计的影响')
% 
% %**************************************************************************
% %AWGN信道下符号个数对信号参数估计性能的影响
% %**************************************************************************
% sig_10 = ofdm(10,para,ratio);
% sig10 = sig_10.*exp(1i*2*pi*fc/fs*(0:length(sig_10)-1));  %上变频到载频 
% [ Tu_rate_10,Ts_rate_10,Tg_rate_10] = Length_rate(sig10,fs,snr,10,0 );
% sig_30 = ofdm(30,para,ratio);
% sig30 = sig_30.*exp(1i*2*pi*fc/fs*(0:length(sig_30)-1));  %上变频到载频 
% [Tu_rate_30,Ts_rate_30,Tg_rate_30] = Length_rate(sig30,fs,snr,30,0 );
% figure
% plot(snr,Tu_rate_10,'r-o');
% hold on
% plot(snr,Tu_rate_128,'k-x');
% hold on
% plot(snr,Tu_rate_30,'b-*');
% legend('10个符号','20个符号','30个符号');
% xlabel('snr/db');
% ylabel('percentage/%');
% title('符号个数对有效数据长度估计的影响')
% figure
% plot(snr,Ts_rate_10,'r-o');
% hold on
% plot(snr,Ts_rate_128,'k-x');
% hold on
% plot(snr,Ts_rate_30,'b-*');
% legend('10个符号','20个符号','30个符号');
% xlabel('snr/db');
% ylabel('percentage/%');
% title('符号个数对符号总长度估计的影响')
% figure
% plot(snr,Tg_rate_10,'r-o');
% hold on
% plot(snr,Tg_rate_128,'k-x');
% hold on
% plot(snr,Tg_rate_30,'b-*');
% legend('10个符号','20个符号','30个符号');
% xlabel('snr/db');
% ylabel('percentage/%');
% title('符号个数对循环前缀长度估计的影响')
% 
% %**************************************************************************
% %AWGN信道下循环前缀比例对信号参数估计性能的影响
% %**************************************************************************
% sig_ratio = ofdm(N,para,1/8);
% sig_r = sig_ratio.*exp(1i*2*pi*fc/fs*(0:length(sig_ratio)-1));  %上变频到载频 
% [ Tu_rate_ratio,Ts_rate_ratio,Tg_rate_ratio] = Length_rate_difratio(sig_r,fs,snr,N,0 );
% 
% sig_ratio1 = ofdm(N,para,3/16);
% sig_r1= sig_ratio1.*exp(1i*2*pi*fc/fs*(0:length(sig_ratio1)-1));  %上变频到载频 
% [Tu_rate_ratio1,Ts_rate_ratio1,Tg_rate_ratio1] = Length_rate_difratio(sig_r1,fs,snr,N,0 );
% figure
% plot(snr,Tu_rate_ratio,'r-o');
% hold on
% plot(snr,Tu_rate_ratio1,'b-*');
% hold on
% plot(snr,Tu_rate_128,'k-x');
% legend('ratio=1/8','ratio=3/16','ratio=1/4');
% xlabel('snr/db');
% ylabel('percentage/%');
% title('循环前缀比例对有效数据长度估计的影响')
% figure
% plot(snr,Ts_rate_ratio,'r-o');
% hold on
% plot(snr,Ts_rate_ratio1,'b-*');
% hold on
% plot(snr,Ts_rate_128,'k-x');
% legend('ratio=1/8','ratio=3/16','ratio=1/4');
% xlabel('snr/db');
% ylabel('percentage/%');
% title('循环前缀比例对符号总长度估计的影响')
% figure
% plot(snr,Tg_rate_ratio,'r-o');
% hold on
% plot(snr,Tg_rate_ratio1,'b-*');
% hold on
% plot(snr,Tg_rate_128,'k-x');
% legend('ratio=1/8','ratio=3/16','ratio=1/4');
% xlabel('snr/db');
% ylabel('percentage/%');
% title('循环前缀比例对循环前缀长度估计的影响')
% 
% %**************************************************************************
% %AWGN信道下OFDM子载波数目估计
% %**************************************************************************
% [carry_num] = carrier_number(sig_a ,N,25,K);
% [Carry_num_rate_128 ] = Rate_Carrynum(sig_a,snr,para,N,K);
% 
% %**************************************************************************
% %AWGN信道下每帧的符号个数对子载波数目估计的影响
% %**************************************************************************
% [Carry_num_rate_10 ] = Rate_Carrynum(sig_10,snr,para,10,0);
% %[Carry_num_rate_30 ] = Rate_Carrynum(sig_30,snr,para,30,0);
% Carry_num_rate_30 = Tu_rate_30 ;
% figure
% plot(snr,Carry_num_rate_10,'r-o');
% hold on
% plot(snr,Carry_num_rate_128,'k-x');
% hold on
% plot(snr,Carry_num_rate_30,'b-*');
% xlabel('snr/db');
% ylabel('percentage/%');
% legend('10个符号','20个符号','30个符号');
% title('符号个数对子载波数目估计的影响')
% 
% %**************************************************************************
% %AWGN信道下循环前缀的比例对子载波数目估计的影响
% %**************************************************************************
% [Carry_num_rate_ratio ] = Rate_Carrynum(sig_ratio,snr,para,N,0);
% [Carry_num_rate_ratio1 ] = Rate_Carrynum(sig_ratio1,snr,para,N,0);
% figure
% plot(snr,Carry_num_rate_ratio,'r-o');
% hold on
% plot(snr,Carry_num_rate_ratio1,'b-*');
% hold on
% plot(snr,Carry_num_rate_128,'k-x');
% xlabel('snr/db');
% ylabel('percentage/%');
% legend('ratio=1/8','ratio=3/16','ratio=1/4');
% title('循环前缀比例对子载波数目估计的影响')
% 
% %**************************************************************************
% %AWGN信道下子载波数目对子载波数目估计的影响
% %**************************************************************************
% %[Carry_num_rate_256 ] = Rate_Carrynum(sig_256,snr,256,N,0);
% % Carry_num_rate_256 = Tu_rate_256;
% % figure
% % plot(snr,Carry_num_rate_256,'r-o');
% % hold on
% % plot(snr,Carry_num_rate_128,'k-x');
% % legend('256个子载波','128个子载波');
% % xlabel('snr/db');
% % ylabel('percentage/%');
% % title('子载波数目对子载波数目估计的影响')
% 
% %**************************************************************************
% %Rayleigh信道下OFDM子载波数目估计准确率
% %**************************************************************************
% [Carry_num_rate_rayleigh ]= Rate_Carrynum_rayleigh(sig_a, snr,N ,para,itau,power,fmax,fs,itn);
% figure
% plot(snr,Carry_num_rate_rayleigh,'r-o');
% hold on
% plot(snr,Carry_num_rate_128,'k-x');
% legend('Rayleigh','Awgn');
% xlabel('snr/db');
% ylabel('percentage/%');
% title('不同信道对子载波数目估计准确率');
% 
% %**************************************************************************
% %利用带宽估计子载波数量
% %**************************************************************************
% [BB, Carrynum_B] =solve_carrynum_rate(sig_baseband,fc,fs,snr,N );    %利用带宽的方式估计子载波数量（基于Tx_data）

end  % 主函数 estimateBandwidth 结束

%**************************************************************************
% 以下为内部函数定义，使test.m可以独立运行
%**************************************************************************

% ==================== PSD_generate ====================
function [ sig_chnl ] = PSD_generate( N)
%生成用于测试的ofdm信号
%N一个帧结构中OFDM信号的个数，每次发送的ofdm信号的个数
%使用窗函数实现符号间平滑过渡，避免频谱泄露

para=128;%设置带内数据子载波数量,有效数据长度
M=16;
Signal=randi([0,M-1],1,para*N);
QAM_out=qammod(Signal,M);  % 使用新的 qammod 函数替代已删除的 modem.qammod
x=reshape(QAM_out,para,N);      %矩阵转换
L=length(x(:,1));
for i=1:N
      x1(:,i)  = [x( 1 : end / 2,i) ;zeros(4*L,1) ;x( end / 2 + 1 : end,i)];  %4倍过采样
end

yy=ifft(x1);                         %对调制后的数据进行ifft变换
yy1=[yy(end-round(length(yy)/4)+1:end,:);yy];%添加循环前缀，比例为1/4

% 移除补零，使用窗函数实现平滑过渡
% 设置过渡窗长度（符号长度的5%）
symbol_len = length(yy1(:,1));
transition_len = max(round(symbol_len * 0.05), 20);  % 过渡区域长度，最小20个点
if transition_len > symbol_len / 4
    transition_len = round(symbol_len / 4);  % 不超过符号长度的1/4
end

% 初始化输出信号
sig_chnl = [];

for i = 1:N
    symbol = yy1(:, i);  % 当前符号
    
    % 创建窗函数实现平滑过渡
    window = ones(length(symbol), 1);
    
    if i == 1
        % 第一个符号：只在末尾加窗（平滑过渡到下一个符号）
        fade_out = hanning(2*transition_len);
        fade_out = fade_out(transition_len+1:end);  % 取后半部分（从0.5到1）
        window(end-transition_len+1:end) = fade_out;
    elseif i == N
        % 最后一个符号：只在开头加窗（从上一个符号平滑过渡）
        fade_in = hanning(2*transition_len);
        fade_in = fade_in(1:transition_len);  % 取前半部分（从0到0.5）
        window(1:transition_len) = fade_in;
    else
        % 中间符号：两端都加窗
        % 开头：从0平滑上升到1（与上一个符号的末尾重叠相加）
        fade_in = hanning(2*transition_len);
        fade_in = fade_in(1:transition_len);
        window(1:transition_len) = fade_in;
        % 末尾：从1平滑下降到0.5（为下一个符号的上升做准备）
        fade_out = hanning(2*transition_len);
        fade_out = fade_out(transition_len+1:end);
        window(end-transition_len+1:end) = fade_out;
    end
    
    % 应用窗函数
    symbol_windowed = symbol .* window;
    
    % 连接符号（平滑过渡，无补零）
    if i == 1
        sig_chnl = symbol_windowed;
    else
        % 重叠相加：将过渡区域重叠相加，实现平滑连接
        overlap_start = length(sig_chnl) - transition_len + 1;
        if overlap_start > 0 && overlap_start <= length(sig_chnl)
            % 重叠区域：上一个符号的末尾 + 当前符号的开头
            sig_chnl(overlap_start:end) = sig_chnl(overlap_start:end) + symbol_windowed(1:transition_len);
            % 添加当前符号的剩余部分
            sig_chnl = [sig_chnl; symbol_windowed(transition_len+1:end)];
        else
            % 如果没有重叠，直接连接
            sig_chnl = [sig_chnl; symbol_windowed];
        end
    end
end

% 转换为行向量
sig_chnl = sig_chnl(:)';
P0=std(sig_chnl);
sig_chnl=sig_chnl/P0;
end

% ==================== Bandwidth_rate_rayleigh ====================
function [ B_rate_welch,B_rate_ar ] = Bandwidth_rate_rayleigh(sig1_chnl,fc,fs, snr,itau,power,fmax,itn,B_ideal)
 %不同信道下的带宽检测率
 %瑞利衰落信道
 % B_ideal: 理想带宽值（Hz），从主函数传入，避免硬编码
 numb = 1; %蒙特卡洛仿真的次数
 LL = length(snr);
 B_welch = zeros(1,LL);
 B_ar = zeros(1,LL);
 B_rate_welch = zeros(1,LL);
 B_rate_ar = zeros(1,LL);
 fprintf('  总进度: 0/%d SNR 值\n', LL);
 fprintf('  参数: 每个SNR值重复 %d 次\n', numb);
 fprintf('  Rayleigh信道参数: fmax=%.1f Hz, 多径数=%d\n', fmax, length(itau));
 fprintf('----------------------------------------\n');
 for i = 1:LL
     snr_start = tic;
     fprintf('  [%d/%d] 处理 SNR = %.1f dB...\n', i, LL, snr(i));
     for j = 1:numb
         [b_welch(i,j),b_ar(i,j)]=PSD_OFDM_rayleigh(sig1_chnl,fc,fs,snr(i),0,itau,power,fmax,itn,B_ideal);
     end
     B_welch(i) = sum(b_welch(i,:))/numb;
     B_ar(i) = sum(b_ar(i,:))/numb;
     B_rate_welch(i) = 1-abs((B_welch(i)-B_ideal))/B_ideal;
     B_rate_ar(i) = 1-abs((B_ar(i)-B_ideal))/B_ideal;
     
     % 计算估计误差
     error_welch_abs = abs(B_welch(i) - B_ideal);
     error_ar_abs = abs(B_ar(i) - B_ideal);
     error_welch_rel = error_welch_abs / B_ideal * 100;
     error_ar_rel = error_ar_abs / B_ideal * 100;
     
     snr_time = toc(snr_start);
     fprintf('    ✓ 完成，耗时: %.2f 秒\n', snr_time);
     fprintf('      Welch算法: %.2f MHz (误差: %.2f MHz, %.2f%%)\n', ...
         B_welch(i)/1e6, error_welch_abs/1e6, error_welch_rel);
     fprintf('      AR模型:    %.2f MHz (误差: %.2f MHz, %.2f%%)\n', ...
         B_ar(i)/1e6, error_ar_abs/1e6, error_ar_rel);
 end
 fprintf('----------------------------------------\n');
 
 % 输出估计误差总结
 fprintf('\n【误差总结】\n');
 fprintf('  理想带宽: %.2f MHz\n', B_ideal/1e6);
 fprintf('  %-8s | %-12s | %-15s | %-15s | %-15s | %-15s\n', ...
     'SNR(dB)', 'Welch估计(MHz)', 'Welch绝对误差(MHz)', 'Welch相对误差(%)', ...
     'AR估计(MHz)', 'AR相对误差(%)');
 fprintf('  %-8s-+-%-12s-+-%-15s-+-%-15s-+-%-15s-+-%-15s\n', ...
     repmat('-',1,8), repmat('-',1,12), repmat('-',1,15), repmat('-',1,15), ...
     repmat('-',1,15), repmat('-',1,15));
 for i = 1:LL
     error_welch_abs = abs(B_welch(i) - B_ideal);
     error_ar_abs = abs(B_ar(i) - B_ideal);
     error_welch_rel = error_welch_abs / B_ideal * 100;
     error_ar_rel = error_ar_abs / B_ideal * 100;
     fprintf('  %-8.1f | %-12.4f | %-15.4f | %-15.2f | %-15.4f | %-15.2f\n', ...
         snr(i), B_welch(i)/1e6, error_welch_abs/1e6, error_welch_rel, ...
         B_ar(i)/1e6, error_ar_rel);
 end
end

% ==================== PSD_OFDM_rayleigh ====================
function [B_welch,B_ar] =PSD_OFDM_rayleigh(sig1_chnl,fc,fs,snr,k,itau,power,fmax,itn,B_ideal)
%k=1时绘制PSD，k=0时不绘制
%B_ideal: 理想带宽值（Hz），从主函数传入（carrier_count × subcarrier_spacing）
sig3_chnl = (sig1_chnl.*exp(1j*2*pi*fc/fs*(0:length(sig1_chnl)-1)));
sig3_chnl = real(MUL_RAYLEIGH(sig3_chnl,itau,power,itn,length(itau),length(sig3_chnl),1/fs,fmax,0));
sig3_chnl = awgn(sig3_chnl,snr,'measured');
%q= Burg(sig1_chnl,'AIC');
%P=mean(abs(sig1_chnl));           %计算信号的标准差，使用函数std来计算，mean和abs函数计算的是不同的标准差
%sig1_chnl=sig1_chnl/P;            %P1和P虽然相同，但一个是除以标准差，一个是除以模型的均值,一个是除以平方和的开方
% [Pxx1,f]=pburg(sig3_chnl,100,4096*2,fs);  %使用AR模型方法进行功率谱估计
[Pxx1, f, p] = Burg(sig3_chnl,fs, 'AIC');   %使用AR模型方法进行功率谱估计
[Pxx2,f1]=pwelch(sig3_chnl,hanning(100),55,4096*2,fs);  %使用welch算法进行功率谱估计

Frc=0:fs/(length(sig3_chnl)):fs/2-1;
OfdmSymComput = 20 * log10(abs(fft(sig3_chnl)));
OfdmSymPSDy = fftshift(OfdmSymComput) - max(OfdmSymComput);

Pxx22 = Pxx2;
Pxx22=Pxx22/min(Pxx22);%找出功率密度中的最小值，归一化处理 
Pxx22=10*log10(Pxx22);%将功率密度转换为dB单位
Pxx22=Pxx22-max(Pxx22);
if k==1
  figure('Name', 'Rayleigh信道-AR模型功率谱密度估计')
  plot(f,Pxx1);
  grid on;
  xlabel('频率 f (Hz)'); 
  ylabel('PSD (dB)');
  title('Rayleigh信道-AR模型方法的功率谱密度估计'); 
  
  figure('Name', 'Rayleigh信道-Welch算法功率谱密度估计')
  plot(f1,Pxx22);
  grid on;
  xlabel('频率 f (Hz)');
  ylabel('PSD (dB)');
  title('Rayleigh信道-Welch算法估计的功率谱密度估计'); 
  
  figure('Name', 'Rayleigh信道-OFDM信号频谱')
  %plot(Frc,OfdmSymPSDy(1,1:end/2));
  plot(Frc,OfdmSymPSDy(1,1:end/2));
  xlabel('频率 f (Hz)');
  ylabel('PSD (dB)');
  title('Rayleigh信道-OFDM信号频谱'); 
end
%****************************
%计算信号的带宽
%****************************

% Welch算法带宽计算：根据SNR自适应选择阈值
L1=ceil(length(Pxx22)/2);
P1=Pxx22(1:L1,1);
P2=Pxx22(L1:end,1);

% 根据SNR选择阈值（与AR模型保持一致）
if snr > 4
    welch_threshold = -6;  % 高SNR使用-6dB
elseif snr > 0 && snr <= 4
    welch_threshold = -5;  % 中等SNR使用-5dB
else
    welch_threshold = -3;  % 低SNR使用-3dB
end

[as1, as11]=Proximate(welch_threshold,P1);  %取最接近阈值处信号的f值
band1=f1(as1);
[as2, as22]=Proximate(welch_threshold,P2);
band2=f1(as2+L1-1);
B_welch =abs(band1-band2);

L2=ceil(length(Pxx1)/2);
P3=Pxx1(1:L2,1);
P4=Pxx1(L2:end,1);
if snr>4
    [as3, as33]=Proximate(-6,P3);  %取最接近-6dB处信号的f值
    band3=f(as3);
    [as4, as44]=Proximate(-6,P4);
    band4=f(as4+L2-1);
    B_ar =abs(band4-band3);
elseif (snr>0&&snr<=4)
    [as3, as33]=Proximate(-5,P3);  %取最接近-5dB处信号的f值
    band3=f(as3);
    [as4, as44]=Proximate(-5,P4);
    band4=f(as4+L2-1);
    B_ar =abs(band4-band3);
else
    [as3, as33]=Proximate(-3,P3);  %取最接近-3dB处信号的f值
    band3=f(as3);
    [as4, as44]=Proximate(-3,P4);
    band4=f(as4+L2-1);
    B_ar =abs(band4-band3);
end
end

% ==================== MUL_RAYLEIGH ====================
function [xout]=MUL_RAYLEIGH(x,itau,dlvl,itn,n1,nsamp,tstp,fd,flat)
%****************** variables *************************
% x  input Ich baseband data     
% yout   output Qch baseband data
% itau   : Delay time for each multipath fading
% dlvl   : Attenuation level for each multipath fading
% itn    : Fading counter for each multipath fading
% n1     : Number of summation for direct and delayed waves 
% nsamp   : Total number od symbols
% tstp   : Mininum time resolution
% fd   : Maxmum doppler frequency
% flat     flat fading or not 
% (1->flat (only amplitude is fluctuated),0->nomal(phase and amplitude are fluctutated)   
%******************************************************
n0 = 25;      % n0     : Number of waves in order to generate each multipath fading
xout = zeros(1,nsamp);
total_attn = sum(10 .^(  dlvl ./ 20.0));

for k = 1 : n1 
    
	atts = 10.^ (dlvl(k)/20.0);

	if dlvl(k) >= 40.0 



    
	       atts = 0.0;
	end

	[xtmp] = delay ( x, nsamp , itau(k));
	[xtmp3] = siglfade (xtmp,nsamp,tstp,fd,n0,itn(k),flat);% single fade
	
  xout = xout + atts .* xtmp3 ./ sqrt(total_attn);

end
end

% ==================== delay (MUL_RAYLEIGH的子函数) ====================
function [xout] = delay(x,nsamp,idel )
% Gives delay to input signal
%****************** variables *************************
% x      input Ich data     
% xout   output Qch data
% nsamp   Number of samples to be simulated 
% idel   Number of samples to be delayed
%******************************************************

xout=zeros(1,nsamp);
if idel ~= 0 
  xout(1:idel) = zeros(1,idel);
end

xout(idel+1:nsamp) = x(1:nsamp-idel);
end

% ==================== siglfade (MUL_RAYLEIGH的子函数) ====================
function [xout]=siglfade(x,nsamp,tstp,fd,no,counter,flat)
% Generate Rayleigh fading
% %****************** variables *************************
% x  : input Ich data     
% xout   : output Qch data
% nsamp  : Number of samples to be simulated       
% tstp   : Minimum time resolution                    
% fd     : maximum doppler frequency               
% no     : number of waves in order to generate fading   
% counter  : fading counter                          
% flat     : flat fading or not 
% (1->flat (only amplitude is fluctuated),0->nomal(phase and amplitude are fluctutated)    
%******************************************************

if fd ~= 0.0  
    ac0 = sqrt(1.0 ./ (2.0.*(no + 1)));   % power normalized constant(ich)
    as0 = sqrt(1.0 ./ (2.0.*no));         % power normalized constant(qch)
    ic0 = counter;                        % fading counter
 
    pai = 3.14159265;   
    wm = 2.0.*pai.*fd;
    n = 4.*no + 2;
    ts = tstp;
    wmts = wm.*ts;
    paino = pai./no;                        

    xc=zeros(1,nsamp);
    xs=zeros(1,nsamp);
    ic=[1:nsamp]+ic0;

  for nn = 1: no
	  cwn = cos( cos(2.0.*pai.*nn./n).*ic.*wmts );
	  xc = xc + cos(paino.*nn).*cwn;
	  xs = xs + sin(paino.*nn).*cwn;
  end

  cwmt = sqrt(2.0).*cos(ic.*wmts);
  xc = (2.0.*xc + cwmt).*ac0;
  xs = 2.0.*xs.*as0;

  ramp=sqrt(xc.^2+xs.^2);   

  if flat ==1
    xout = sqrt(xc.^2+xs.^2).*x;    % output signal
  else
    xout = x .*(xc+1i*xs);
  end

else  
  xout = x;
end
end

% ==================== Burg ====================
function [psdviaBurg, f, p] = Burg(x, Fs, varargin)
%MYBURG      使用burg算法实现的AR模型功率谱估计
% psdviaBurg 使用burg算法计算的功率谱值
% f          频率采样点
% p          模型阶数
% x          输入信号
% Fs         采样频率
% varargin   可以为数值型，即为AR模型阶数
%            可以为字符串，即为准则准则AR模型阶数由准则确定
%
% 根据输入参数类型判断
if strcmp(class(varargin{1}), 'double')
    p = varargin{1};
elseif ischar(varargin{1})
    criterion = varargin{1};
else
    error('第2个参数必须为数值型或字符串');
end
x = x(:);
N = length(x);
% 模型参数估计
if exist('p', 'var') % p变量是否存在，如果存在则不需要估计，直接使用p值
    [a, E] = computeARpara(x, p);
else % p不存在，需要估计，根据准则criterion
    % 限制最大AR模型阶数，避免对长信号计算量过大
    p_max = 2000;  % 最大AR模型阶数（可根据需要调整）
    p = min(ceil(N/3), p_max); % 阶数一般不超出信号长度的1/3，但不超过p_max
    
    % 计算1到p阶的误差
    [a, E] = computeARpara(x, p);
    
    % 计算目标函数的最小值
    kc = 1:p + 1;
    switch criterion
        case 'FPE'
            goalF = E.*(N + (kc + 1))./(N - (kc + 1));
        case 'AIC'
            goalF = N.*log(E) + 2.*kc;
    end
    [minF, p] = min(goalF); % p是目标函数最小值的位值，也即准则准则确定的阶数
    
    % 使用p值重新计算AR模型参数
    [a, E] = computeARpara(x, p);
end

% 优化：根据信号长度和采样频率动态计算合理的频率点数
% 原代码使用20e5=200,000个点，对于带宽估计来说过度采样
% 使用2^nextpow2(length(x))*4，但不超过50,000点，保证精度和效率平衡
N_freq = min(2^nextpow2(length(x)) * 4, 50000);
[h, f] = freqz(1, a, N_freq, Fs);
psdviaBurg = E(end)*abs(h).^2./Fs;
psdviaBurg=psdviaBurg/abs(max(psdviaBurg));
psdviaBurg=(10*log10(abs(psdviaBurg)));
end

% ==================== computeARpara (Burg的子函数) ====================
function [a, E] = computeARpara(x, p)
% 根据输入信号x和阶数p计算AR模型参数估计
N = length(x);
% 初始值
ef = x; % 前向预测误差
eb = x; % 后向预测误差
a  = 1; % 初始模型参数
E  = x'*x/N; % 初始误差
k  = zeros(1, p); % 为反射系数预分配空间，加快循环速度
E  = [E k]; % 为误差预分配空间，加快速度
for m = 1:p
    % 按照burg算法步骤，首先计算m阶的反射系数
    efm = ef(2:end); % 前一阶次的前向预测误差
    ebm = eb(1:end - 1); % 前一阶次的后向预测误差
    num = -2.*ebm'*efm;  % 反射系数的分子项
    den = efm'*efm + ebm'*ebm; % 反射系数的分母项
    k(m) = num./den; % 当前阶次的反射系数
    
    % 更新前后向预测误差
    ef = efm + k(m)*ebm;
    eb = ebm + conj(k(m))*efm;
    
    % 更新模型系数a
    a = [a; 0] + k(m)*[0; conj(flipud(a))];
    
    % 当前阶次的误差功率
    E(m + 1) = (1 - conj(k(m))*k(m))*E(m);
end
end

% ==================== Proximate ====================
function [as1, as11]=Proximate(b,aa)
%%% b是要找的参数值
%%% aa是要找的数组

%%判断aa是几维数组，统一转换为一维
a=aa(:);                                            %%将数组转换为化为一维数组
ab=(a(:)-b)';                                     %%计算数组a和b的差值
abc=abs(ab);                  
abc=sort(abc);                                  %%绝对值取最小值，排序

%%%  二维数组as存储b值最接近的值，在abc中的位置 
%    找到与b值最接近的第一个数组元素的位置（最接近b的）
[as1 ,as11] =find(abs((a(:)-b))==abc(1,1));
end
