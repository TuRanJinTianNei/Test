function [B_welch, B_ar, B_ideal, results] = estimateBandwidth(Tx_data, fs, carrier_count, subcarrier_spacing, varargin)
%===============================================================================
% estimateBandwidth copy.m - OFDM信号带宽估计函数（使用发送信号Tx_data）
% 
% 功能说明:
%   使用Welch算法和AR模型法估计OFDM信号的带宽
%   - 使用signalGenerate.m生成的发送信号（Tx_data）
%   - 对发送信号进行上变频和加噪声处理，然后进行带宽估计
%   - 可独立运行，自动调用signalGenerate.m生成信号
%
% 输入参数:
%   Tx_data          - OFDM发送信号（基带复数信号，行向量）
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
% 依赖关系:
%   - 如果工作空间没有Tx_data，会自动调用signalGenerate.m生成信号
%   - 需要相关函数（Burg.m, Proximate.m等）
%   - 需要Signal Processing Toolbox（用于pwelch、awgn等函数）
%
% 信号处理流程:
%   1. 获取Tx_data（发送信号，基带复数信号）
%   2. 上变频：real(Tx_data * exp(j*2*pi*fc/fs*t))
%   3. 加噪声：awgn(sig_processed, snr, 'measured')
%   4. 使用已处理的信号（sig_processed）直接进行功率谱估计和带宽计算
%      - 使用Welch算法和AR模型法进行功率谱估计
%      - 计算-3dB/-4dB/-5dB带宽（根据SNR选择阈值）
%
% 使用示例:
%   1. 无参数调用（自动调用signalGenerate.m生成Tx_data）:
%      estimateBandwidth
%
%   2. 完整参数调用:
%      [B_welch, B_ar, B_ideal, results] = estimateBandwidth(Tx_data, fs, ...
%          carrier_count, subcarrier_spacing, 'fc', 0, 'snr', 20, 'plot', true);
%
% 创建日期: 2025.12.10
% 基于: signalGenerate.m的带宽估计部分
% 修改说明: 使用signalGenerate.m的Tx_data（发送信号），进行上变频和加噪声后估计带宽
%===============================================================================

% 处理无参数调用情况：从工作空间获取参数或自动调用signalGenerate
if nargin == 0
    fprintf('检测到无参数调用，尝试从工作空间获取参数...\n');
    
    % 尝试从基础工作空间获取Tx_data
    tx_data_found = false;
    try
        Tx_data = evalin('base', 'Tx_data');
        if ~isempty(Tx_data)
            tx_data_found = true;
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
            end
        catch
            % Tx_data不存在
        end
    end
    
    % 如果仍然没有找到Tx_data，自动调用signalGenerate.m生成信号
    if ~tx_data_found
        fprintf('工作空间中未找到Tx_data变量，自动调用signalGenerate.m生成信号...\n');
        
        % 检查signalGenerate.m是否存在
        if exist('signalGenerate.m', 'file') == 2
            fprintf('正在运行signalGenerate.m...\n');
            fprintf('========================================\n');
            
            % 保存当前工作空间状态（避免signalGenerate清除变量）
            % 注意：signalGenerate.m中有clear all，会清除所有变量
            % 我们需要在调用后重新获取变量
            
            % 运行signalGenerate.m（在基础工作空间中）
            try
                fprintf('正在调用signalGenerate.m生成发送信号（Tx_data）...\n');
                % 使用evalin在基础工作空间执行signalGenerate
                evalin('base', 'signalGenerate');
                
                % 等待一下确保signalGenerate执行完成
                pause(0.1);
                
                % 从基础工作空间获取生成的发送信号（Tx_data）
                % 注意：signalGenerate.m总是会生成Tx_data，无论use_channel设置如何
                Tx_data = evalin('base', 'Tx_data');
                if isempty(Tx_data)
                    error('Tx_data为空');
                end
                
                % 同时获取其他参数（signalGenerate会生成这些变量）
                try
                    fs = evalin('base', 'fs');
                    if isempty(fs)
                        error('fs为空');
                    end
                catch
                    fprintf('警告：未找到fs变量，使用默认值15.36e6 Hz\n');
                    fs = 15.36e6;  % 默认值
                end
                
                try
                    carrier_count = evalin('base', 'carrier_count');
                    if isempty(carrier_count)
                        error('carrier_count为空');
                    end
                catch
                    fprintf('警告：未找到carrier_count变量，使用默认值300\n');
                    carrier_count = 300;  % 默认值
                end
                
                try
                    subcarrier_spacing = evalin('base', 'subcarrier_spacing');
                    if isempty(subcarrier_spacing)
                        error('subcarrier_spacing为空');
                    end
                catch
                    fprintf('警告：未找到subcarrier_spacing变量，使用默认值15e3 Hz\n');
                    subcarrier_spacing = 15e3;  % 默认值
                end
                
                fprintf('信号生成完成！\n');
                fprintf('  - 使用信号: Tx_data（发送信号，来自signalGenerate.m）\n');
                fprintf('  - Tx_data长度: %d 样本\n', length(Tx_data));
                fprintf('  - fs: %.2e Hz\n', fs);
                fprintf('  - carrier_count: %d\n', carrier_count);
                fprintf('  - subcarrier_spacing: %.1f Hz\n', subcarrier_spacing);
                fprintf('========================================\n\n');
                
                % 标记参数已获取，跳过后续的参数获取步骤
                skip_param_fetch = true;
            catch ME
                error('错误：自动调用signalGenerate.m失败。\n错误信息：%s\n\n请手动运行signalGenerate.m生成信号，或提供Tx_data参数。', ME.message);
            end
        else
            error('错误：工作空间中未找到Tx_data变量，且signalGenerate.m文件不存在。\n请先运行signalGenerate.m生成信号，或提供Tx_data参数。');
        end
    else
        fprintf('成功从工作空间获取Tx_data\n');
    end
    
    % 如果自动调用了signalGenerate，参数已经获取，跳过后续步骤
    if ~exist('skip_param_fetch', 'var') || ~skip_param_fetch
        % 尝试从工作空间获取fs（优先从基础工作空间获取，因为signalGenerate在基础工作空间运行）
        try
            fs = evalin('base', 'fs');
            if isempty(fs)
                error('fs为空');
            end
        catch
            try
                fs = evalin('caller', 'fs');
                if isempty(fs)
                    error('fs为空');
                end
            catch
                fprintf('警告：工作空间中未找到fs变量，使用默认值15.36e6 Hz\n');
                fs = 15.36e6;
            end
        end
        
        % 尝试从工作空间获取carrier_count
        try
            carrier_count = evalin('base', 'carrier_count');
            if isempty(carrier_count)
                error('carrier_count为空');
            end
        catch
            try
                carrier_count = evalin('caller', 'carrier_count');
                if isempty(carrier_count)
                    error('carrier_count为空');
                end
            catch
                fprintf('警告：工作空间中未找到carrier_count变量，使用默认值300\n');
                carrier_count = 300;
            end
        end
        
        % 尝试从工作空间获取subcarrier_spacing
        try
            subcarrier_spacing = evalin('base', 'subcarrier_spacing');
            if isempty(subcarrier_spacing)
                error('subcarrier_spacing为空');
            end
        catch
            try
                subcarrier_spacing = evalin('caller', 'subcarrier_spacing');
                if isempty(subcarrier_spacing)
                    error('subcarrier_spacing为空');
                end
            catch
                fprintf('警告：工作空间中未找到subcarrier_spacing变量，使用默认值15e3 Hz\n');
                subcarrier_spacing = 15e3;
            end
        end
        
        fprintf('成功从工作空间获取参数：\n');
        fprintf('  - Tx_data长度: %d\n', length(Tx_data));
        fprintf('  - fs: %.2e Hz\n', fs);
        fprintf('  - carrier_count: %d\n', carrier_count);
        fprintf('  - subcarrier_spacing: %.1f Hz\n\n', subcarrier_spacing);
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

% 计算理论带宽
B_ideal = carrier_count * subcarrier_spacing;

% 输出仿真参数
fprintf('\n========================================\n');
fprintf('OFDM信号带宽估计仿真开始\n');
fprintf('========================================\n');
fprintf('开始时间: %s\n', start_time_str);
fprintf('\n【输入参数】\n');
fprintf('  信号长度: %d 样本\n', length(Tx_data));
fprintf('  采样频率: %.2e Hz (%.2f MHz)\n', fs, fs/1e6);
fprintf('  有效数据子载波数: %d\n', carrier_count);
fprintf('  子载波间隔: %.1f Hz (%.3f kHz)\n', subcarrier_spacing, subcarrier_spacing/1e3);
fprintf('  载波频率: %.6e Hz (%.3f MHz)\n', fc, fc/1e6);
fprintf('  信噪比: %.1f dB\n', snr_estimate);
fprintf('  是否绘图: %s\n', mat2str(plot_flag));
fprintf('\n【理论参数】\n');
fprintf('  理论带宽: %.6e Hz (%.3f MHz)\n', B_ideal, B_ideal/1e6);
fprintf('========================================\n\n');

% 处理信号：上变频和加噪声（参考PSD_OFDM.m的处理方式）
fprintf('步骤1: 处理信号（上变频和加噪声）...\n');
fprintf('  - 输入信号: Tx_data（基带复数信号，来自signalGenerate.m）\n');
fprintf('  - 信号长度: %d 样本\n', length(Tx_data));
process_start = tic;

fprintf('  1.1: 上变频处理...\n');
upconversion_start = tic;
sig_processed = real(Tx_data.*exp(1j*2*pi*fc/fs*(0:length(Tx_data)-1)));
upconversion_time = toc(upconversion_start);
fprintf('      - 载波频率: %.6e Hz (%.3f MHz)\n', fc, fc/1e6);
fprintf('      - 采样频率: %.2e Hz (%.2f MHz)\n', fs, fs/1e6);
if fc == 0
    fprintf('      - 说明: 基带信号（fc=0），上变频后仍为基带\n');
else
    fprintf('      - 说明: 已上变频到载波频率\n');
end
fprintf('      - 处理时间: %.3f 秒\n', upconversion_time);

fprintf('  1.2: 加噪声处理...\n');
noise_start = tic;
sig_processed = awgn(sig_processed, snr_estimate, 'measured');
noise_time = toc(noise_start);
fprintf('      - 信噪比: %.1f dB\n', snr_estimate);
fprintf('      - 噪声类型: 加性高斯白噪声（AWGN）\n');
fprintf('      - 功率测量: 基于信号功率自动测量\n');
fprintf('      - 处理时间: %.3f 秒\n', noise_time);

process_time = toc(process_start);
fprintf('\n  步骤1完成！\n');
fprintf('  - 输出信号: sig_processed（已上变频和加噪声的实信号）\n');
fprintf('  - 信号长度: %d 样本\n', length(sig_processed));
fprintf('  - 总处理时间: %.2f 秒\n', process_time);
fprintf('    * 上变频: %.3f 秒\n', upconversion_time);
fprintf('    * 加噪声: %.3f 秒\n', noise_time);
fprintf('\n');

% 使用已处理的信号（sig_processed）进行带宽估计
fprintf('步骤2: 使用Welch算法和AR模型法估计带宽...\n');
fprintf('  - 使用已上变频和加噪声的信号（sig_processed）进行估计\n');
fprintf('  - 信号长度: %d 样本\n', length(sig_processed));
estimation_start = tic;
try
    % 使用已处理的信号（sig_processed）进行功率谱估计
    % sig_processed已经完成了上变频和加噪声处理
    
    fprintf('  2.1: AR模型功率谱估计...\n');
    ar_start = tic;
    [Pxx1, f, p] = Burg(sig_processed, fs, 'AIC');
    ar_time = toc(ar_start);
    fprintf('      - AR模型阶数: %d (AIC准则)\n', p);
    fprintf('      - 功率谱点数: %d\n', length(Pxx1));
    fprintf('      - 估计时间: %.3f 秒\n', ar_time);
    
    fprintf('  2.2: Welch算法功率谱估计...\n');
    welch_start = tic;
    [Pxx2, f1] = pwelch(sig_processed, hanning(100), 55, 4096*2, fs);
    welch_time = toc(welch_start);
    fprintf('      - 窗函数: Hanning(100)\n');
    fprintf('      - 重叠: 55样本\n');
    fprintf('      - FFT点数: %d\n', 4096*2);
    fprintf('      - 功率谱点数: %d\n', length(Pxx2));
    fprintf('      - 估计时间: %.3f 秒\n', welch_time);
    
    fprintf('  2.3: 归一化功率谱...\n');
    % 归一化Welch功率谱
    Pxx22 = Pxx2;
    Pxx22 = Pxx22 / min(Pxx22);  % 归一化处理
    Pxx22 = 10 * log10(Pxx22);   % 转换为dB单位
    Pxx22 = Pxx22 - max(Pxx22);
    fprintf('      - 完成归一化（dB单位）\n');
    
    fprintf('  2.4: 分析功率谱特性并自适应选择阈值...\n');
    % 分析功率谱特性，自适应选择阈值
    % 对于Welch算法功率谱
    psd_max_welch = max(Pxx22);
    psd_min_welch = min(Pxx22);
    psd_range_welch = psd_max_welch - psd_min_welch;
    % 计算功率谱的动态范围（峰值到谷值的差值）
    % 如果动态范围较大（>15dB），说明信号有明显的主瓣和旁瓣
    % 如果动态范围较小（<10dB），说明信号功率分布较均匀
    
    % 对于AR模型功率谱
    psd_max_ar = max(Pxx1);
    psd_min_ar = min(Pxx1);
    psd_range_ar = psd_max_ar - psd_min_ar;
    
    fprintf('      - Welch功率谱范围: %.2f dB 到 %.2f dB (动态范围: %.2f dB)\n', ...
        psd_min_welch, psd_max_welch, psd_range_welch);
    fprintf('      - AR功率谱范围: %.2f dB 到 %.2f dB (动态范围: %.2f dB)\n', ...
        psd_min_ar, psd_max_ar, psd_range_ar);
    
    % 自适应选择Welch算法阈值
    % 根据功率谱的动态范围和SNR选择阈值
    if psd_range_welch > 15
        % 动态范围大，使用较宽松的阈值（-6dB或-10dB）
        if snr_estimate > 10
            welch_threshold = -6;
        elseif snr_estimate > 5
            welch_threshold = -5;
        else
            welch_threshold = -4;
        end
    elseif psd_range_welch > 10
        % 中等动态范围，使用中等阈值
        if snr_estimate > 10
            welch_threshold = -5;
        else
            welch_threshold = -4;
        end
    else
        % 动态范围小，使用标准阈值
        welch_threshold = -3;
    end
    fprintf('      - Welch算法自适应阈值: %d dB (基于动态范围%.2f dB和SNR %.1f dB)\n', ...
        welch_threshold, psd_range_welch, snr_estimate);
    
    % 自适应选择AR模型阈值
    % 根据功率谱的动态范围和SNR选择阈值
    if psd_range_ar > 15
        % 动态范围大，使用较宽松的阈值
        if snr_estimate > 10
            ar_threshold = -6;
        elseif snr_estimate > 5
            ar_threshold = -5;
        else
            ar_threshold = -4;
        end
    elseif psd_range_ar > 10
        % 中等动态范围
        if snr_estimate > 10
            ar_threshold = -5;
        elseif snr_estimate > 4
            ar_threshold = -4;
        else
            ar_threshold = -3;
        end
    else
        % 动态范围小，根据SNR选择
        if snr_estimate > 4
            ar_threshold = -5;
        elseif snr_estimate > 0
            ar_threshold = -4;
        else
            ar_threshold = -3;
        end
    end
    fprintf('      - AR模型自适应阈值: %d dB (基于动态范围%.2f dB和SNR %.1f dB)\n', ...
        ar_threshold, psd_range_ar, snr_estimate);
    
    fprintf('  2.5: 计算OFDM信号频谱（用于绘图）...\n');
    % 计算OFDM信号频谱（用于绘图）
    Frc = 0:fs/(length(sig_processed)):fs/2-1;
    OfdmSymComput = 20 * log10(abs(fft(sig_processed)));
    OfdmSymPSDy = fftshift(OfdmSymComput) - max(OfdmSymComput);
    fprintf('      - 完成\n');
    
    % 绘制功率谱图（如果需要）
    if plot_flag
        fprintf('  2.6: 绘制功率谱图...\n');
        figure('Name', 'AR模型功率谱密度估计')
        plot(f, Pxx1);
        grid on;
        xlabel('频率 f (Hz)');
        ylabel('PSD (dB)');
        title(sprintf('AR模型方法的功率谱密度估计 (阈值: %d dB)', ar_threshold));
        hold on;
        plot([f(1), f(end)], [ar_threshold, ar_threshold], 'r--', 'LineWidth', 1.5);
        legend('功率谱', sprintf('阈值 (%d dB)', ar_threshold), 'Location', 'best');
        hold off;
        
        figure('Name', 'Welch算法功率谱密度估计')
        plot(f1, Pxx22);
        grid on;
        xlabel('频率 f (Hz)');
        ylabel('PSD (dB)');
        title(sprintf('Welch算法估计的功率谱密度估计 (阈值: %d dB)', welch_threshold));
        hold on;
        plot([f1(1), f1(end)], [welch_threshold, welch_threshold], 'r--', 'LineWidth', 1.5);
        legend('功率谱', sprintf('阈值 (%d dB)', welch_threshold), 'Location', 'best');
        hold off;
        
        figure('Name', 'OFDM信号频谱')
        plot(Frc, OfdmSymPSDy(1, 1:end/2));
        xlabel('频率 f (Hz)');
        ylabel('PSD (dB)');
        title('OFDM信号频谱');
        fprintf('      - 已生成3个功率谱图（包含阈值线）\n');
    end
    
    fprintf('  2.7: 计算Welch算法带宽（自适应阈值）...\n');
    % 计算Welch算法带宽（使用自适应阈值）
    L1 = ceil(length(Pxx22) / 2);
    P1 = Pxx22(1:L1, 1);
    P2 = Pxx22(L1:end, 1);
    [as1, as11] = Proximate(welch_threshold, P1);  % 取最接近阈值的f值
    band1 = f1(as1);
    [as2, as22] = Proximate(welch_threshold, P2);
    band2 = f1(as2 + L1 - 1);
    B_welch = abs(band1 - band2);
    fprintf('      - 阈值: %d dB (自适应选择)\n', welch_threshold);
    fprintf('      - 下边带频率: %.6e Hz (%.3f MHz)\n', band1, band1/1e6);
    fprintf('      - 上边带频率: %.6e Hz (%.3f MHz)\n', band2, band2/1e6);
    fprintf('      - 估计带宽: %.6e Hz (%.3f MHz)\n', B_welch, B_welch/1e6);
    
    fprintf('  2.8: 计算AR模型带宽（自适应阈值）...\n');
    % 计算AR模型带宽（使用自适应阈值）
    L2 = ceil(length(Pxx1) / 2);
    P3 = Pxx1(1:L2, 1);
    P4 = Pxx1(L2:end, 1);
    [as3, as33] = Proximate(ar_threshold, P3);  % 取最接近阈值的f值
    band3 = f(as3);
    [as4, as44] = Proximate(ar_threshold, P4);
    band4 = f(as4 + L2 - 1);
    B_ar = abs(band4 - band3);
    threshold_db = ar_threshold;  % 保存用于后续输出
    fprintf('      - 阈值: %d dB (自适应选择)\n', ar_threshold);
    fprintf('      - 下边带频率: %.6e Hz (%.3f MHz)\n', band3, band3/1e6);
    fprintf('      - 上边带频率: %.6e Hz (%.3f MHz)\n', band4, band4/1e6);
    fprintf('      - 估计带宽: %.6e Hz (%.3f MHz)\n', B_ar, B_ar/1e6);
    
    estimation_time = toc(estimation_start);
    
    % 计算估计误差
    error_welch_abs = abs(B_welch - B_ideal);
    error_welch_rel = (error_welch_abs / B_ideal) * 100;
    
    error_ar_abs = abs(B_ar - B_ideal);
    error_ar_rel = (error_ar_abs / B_ideal) * 100;
    
    fprintf('\n  步骤2完成！\n');
    fprintf('  - 总估计时间: %.2f 秒\n', estimation_time);
    fprintf('    * AR模型估计: %.3f 秒\n', ar_time);
    fprintf('    * Welch算法估计: %.3f 秒\n', welch_time);
    fprintf('    * 带宽计算: %.3f 秒\n', estimation_time - ar_time - welch_time);
    fprintf('\n');
    
    % 计算总仿真时间
    total_time = toc(sim_start_time);
    end_time_str = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    
    % 输出估计结果
    fprintf('\n========================================\n');
    fprintf('带宽估计结果\n');
    fprintf('========================================\n');
    fprintf('结束时间: %s\n', end_time_str);
    fprintf('总仿真时间: %.2f 秒 (%.2f 分钟)\n', total_time, total_time/60);
    fprintf('\n【输入参数总结】\n');
    fprintf('  原始信号长度: %d 样本 (Tx_data)\n', length(Tx_data));
    fprintf('  处理后信号长度: %d 样本 (sig_processed)\n', length(sig_processed));
    fprintf('  采样频率: %.2e Hz (%.2f MHz)\n', fs, fs/1e6);
    fprintf('  有效数据子载波数: %d\n', carrier_count);
    fprintf('  子载波间隔: %.1f Hz (%.3f kHz)\n', subcarrier_spacing, subcarrier_spacing/1e3);
    fprintf('  载波频率: %.6e Hz (%.3f MHz)\n', fc, fc/1e6);
    fprintf('  信噪比: %.1f dB\n', snr_estimate);
    
    fprintf('\n【信号处理说明】\n');
    fprintf('  1. 使用signalGenerate.m生成的发送信号（Tx_data，基带复数信号）\n');
    fprintf('  2. 上变频处理: real(Tx_data * exp(j*2*pi*fc/fs*t))\n');
    fprintf('  3. 加噪声处理: awgn(sig_processed, SNR, ''measured'')\n');
    fprintf('  4. 使用处理后的信号（sig_processed）进行带宽估计\n');
    
    fprintf('\n【理论带宽】\n');
    fprintf('  理论带宽: %.6e Hz (%.3f MHz)\n', B_ideal, B_ideal/1e6);
    fprintf('  计算公式: carrier_count * subcarrier_spacing\n');
    
    fprintf('\n【功率谱特性分析】\n');
    fprintf('  Welch功率谱动态范围: %.2f dB (峰值: %.2f dB, 谷值: %.2f dB)\n', ...
        psd_range_welch, psd_max_welch, psd_min_welch);
    fprintf('  AR功率谱动态范围: %.2f dB (峰值: %.2f dB, 谷值: %.2f dB)\n', ...
        psd_range_ar, psd_max_ar, psd_min_ar);
    fprintf('  说明: 根据功率谱动态范围和SNR自适应选择阈值\n');
    
    fprintf('\n【Welch算法估计结果】\n');
    fprintf('  估计带宽: %.6e Hz (%.3f MHz)\n', B_welch, B_welch/1e6);
    fprintf('  带宽阈值: %d dB (自适应选择，基于动态范围%.2f dB和SNR %.1f dB)\n', ...
        welch_threshold, psd_range_welch, snr_estimate);
    fprintf('  下边带频率: %.6e Hz (%.3f MHz)\n', band1, band1/1e6);
    fprintf('  上边带频率: %.6e Hz (%.3f MHz)\n', band2, band2/1e6);
    fprintf('  绝对误差: %.6e Hz (%.3f MHz)\n', error_welch_abs, error_welch_abs/1e6);
    fprintf('  相对误差: %.2f%%\n', error_welch_rel);
    % 评估估计精度
    if error_welch_rel < 5
        welch_precision = '优秀 (<5%%)';
    elseif error_welch_rel < 10
        welch_precision = '良好 (<10%%)';
    else
        welch_precision = '一般 (>=10%%)';
    end
    fprintf('  估计精度: %s\n', welch_precision);
    
    fprintf('\n【AR模型法估计结果】\n');
    fprintf('  估计带宽: %.6e Hz (%.3f MHz)\n', B_ar, B_ar/1e6);
    fprintf('  带宽阈值: %d dB (自适应选择，基于动态范围%.2f dB和SNR %.1f dB)\n', ...
        threshold_db, psd_range_ar, snr_estimate);
    fprintf('  下边带频率: %.6e Hz (%.3f MHz)\n', band3, band3/1e6);
    fprintf('  上边带频率: %.6e Hz (%.3f MHz)\n', band4, band4/1e6);
    fprintf('  AR模型阶数: %d (AIC准则)\n', p);
    fprintf('  绝对误差: %.6e Hz (%.3f MHz)\n', error_ar_abs, error_ar_abs/1e6);
    fprintf('  相对误差: %.2f%%\n', error_ar_rel);
    % 评估估计精度
    if error_ar_rel < 5
        ar_precision = '优秀 (<5%%)';
    elseif error_ar_rel < 10
        ar_precision = '良好 (<10%%)';
    else
        ar_precision = '一般 (>=10%%)';
    end
    fprintf('  估计精度: %s\n', ar_precision);
    
    fprintf('\n【方法对比】\n');
    if error_welch_abs < error_ar_abs
        fprintf('  Welch算法估计更准确（误差更小）\n');
        fprintf('  误差差异: %.6e Hz (%.2f%%)\n', ...
            error_ar_abs - error_welch_abs, error_ar_rel - error_welch_rel);
        better_method = 'Welch';
    elseif error_ar_abs < error_welch_abs
        fprintf('  AR模型法估计更准确（误差更小）\n');
        fprintf('  误差差异: %.6e Hz (%.2f%%)\n', ...
            error_welch_abs - error_ar_abs, error_welch_rel - error_ar_rel);
        better_method = 'AR';
    else
        fprintf('  两种方法估计精度相同\n');
        better_method = 'Equal';
    end
    
    fprintf('\n【时间统计】\n');
    fprintf('  信号处理时间: %.2f 秒\n', process_time);
    fprintf('  带宽估计时间: %.2f 秒\n', estimation_time);
    fprintf('  总仿真时间: %.2f 秒 (%.2f 分钟)\n', total_time, total_time/60);
    fprintf('========================================\n\n');
    
    % 构建结果结构体
    results.B_ideal = B_ideal;
    results.B_welch = B_welch;
    results.B_ar = B_ar;
    results.error_welch_abs = error_welch_abs;
    results.error_welch_rel = error_welch_rel;
    results.error_ar_abs = error_ar_abs;
    results.error_ar_rel = error_ar_rel;
    results.better_method = better_method;
    results.fc = fc;
    results.snr = snr_estimate;
    results.fs = fs;
    results.carrier_count = carrier_count;
    results.subcarrier_spacing = subcarrier_spacing;
    results.signal_length = length(Tx_data);
    results.process_time = process_time;
    results.estimation_time = estimation_time;
    results.total_time = total_time;
    results.start_time = start_time_str;
    results.end_time = end_time_str;
    % 自适应阈值信息
    results.welch_threshold = welch_threshold;
    results.ar_threshold = ar_threshold;
    results.psd_range_welch = psd_range_welch;
    results.psd_range_ar = psd_range_ar;
    results.psd_max_welch = psd_max_welch;
    results.psd_min_welch = psd_min_welch;
    results.psd_max_ar = psd_max_ar;
    results.psd_min_ar = psd_min_ar;
    results.band1 = band1;
    results.band2 = band2;
    results.band3 = band3;
    results.band4 = band4;
    
catch ME
    estimation_time = toc(estimation_start);
    fprintf('  - 估计时间: %.2f 秒\n', estimation_time);
    fprintf('  完成！\n\n');
    
    fprintf('错误：带宽估计失败！\n');
    fprintf('错误信息：%s\n', ME.message);
    fprintf('可能原因：\n');
    fprintf('  1. 缺少必要的函数（Burg.m, Proximate.m等）\n');
    fprintf('  2. 缺少必要的工具箱（Signal Processing Toolbox等）\n');
    fprintf('  3. 信号处理过程中出现错误\n');
    fprintf('\n尝试使用简化方法进行带宽估计...\n');
    
    % 如果PSD_OFDM不可用，使用简化的方法
    try
        simplified_start = tic;
        % 使用Welch算法
        [Pxx2, f1] = pwelch(sig_processed, hanning(100), 55, 4096*2, fs);
        Pxx22 = Pxx2;
        Pxx22 = Pxx22 / min(Pxx22);
        Pxx22 = 10 * log10(Pxx22);
        Pxx22 = Pxx22 - max(Pxx22);
        
        % 找到-3dB点
        L1 = ceil(length(Pxx22) / 2);
        P1 = Pxx22(1:L1, 1);
        P2 = Pxx22(L1:end, 1);
        
        % 找到最接近-3dB的点
        [~, idx1] = min(abs(P1 - (-3)));
        [~, idx2] = min(abs(P2 - (-3)));
        band1 = f1(idx1);
        band2 = f1(idx2 + L1 - 1);
        B_welch = abs(band2 - band1);
        B_ar = [];  % AR模型需要Burg函数，无法使用简化方法
        
        simplified_time = toc(simplified_start);
        fprintf('简化Welch算法估计带宽: %.6e Hz (%.3f MHz)\n', ...
            B_welch, B_welch/1e6);
        fprintf('简化方法执行时间: %.2f 秒\n', simplified_time);
        
        % 计算估计误差
        error_welch_abs = abs(B_welch - B_ideal);
        error_welch_rel = (error_welch_abs / B_ideal) * 100;
        
        % 计算总仿真时间
        total_time = toc(sim_start_time);
        end_time_str = datestr(now, 'yyyy-mm-dd HH:MM:SS');
        
        % 输出简化结果
        fprintf('\n========================================\n');
        fprintf('带宽估计结果（简化方法）\n');
        fprintf('========================================\n');
        fprintf('结束时间: %s\n', end_time_str);
        fprintf('总仿真时间: %.2f 秒 (%.2f 分钟)\n', total_time, total_time/60);
        fprintf('\n【输入参数总结】\n');
        fprintf('  信号长度: %d 样本\n', length(Tx_data));
        fprintf('  采样频率: %.2e Hz (%.2f MHz)\n', fs, fs/1e6);
        fprintf('  有效数据子载波数: %d\n', carrier_count);
        fprintf('  子载波间隔: %.1f Hz (%.3f kHz)\n', subcarrier_spacing, subcarrier_spacing/1e3);
        fprintf('  载波频率: %.6e Hz (%.3f MHz)\n', fc, fc/1e6);
        fprintf('  信噪比: %.1f dB\n', snr_estimate);
        fprintf('\n【理论带宽】\n');
        fprintf('  理论带宽: %.6e Hz (%.3f MHz)\n', B_ideal, B_ideal/1e6);
        fprintf('\n【Welch算法估计结果（简化方法）】\n');
        fprintf('  估计带宽: %.6e Hz (%.3f MHz)\n', B_welch, B_welch/1e6);
        fprintf('  绝对误差: %.6e Hz (%.3f MHz)\n', error_welch_abs, error_welch_abs/1e6);
        fprintf('  相对误差: %.2f%%\n', error_welch_rel);
        fprintf('\n【时间统计】\n');
        fprintf('  信号处理时间: %.2f 秒\n', process_time);
        fprintf('  带宽估计时间: %.2f 秒\n', estimation_time);
        fprintf('  简化方法时间: %.2f 秒\n', simplified_time);
        fprintf('  总仿真时间: %.2f 秒 (%.2f 分钟)\n', total_time, total_time/60);
        fprintf('========================================\n\n');
        
        % 构建结果结构体（简化版本）
        results.B_ideal = B_ideal;
        results.B_welch = B_welch;
        results.B_ar = [];
        results.error_welch_abs = error_welch_abs;
        results.error_welch_rel = error_welch_rel;
        results.error_ar_abs = [];
        results.error_ar_rel = [];
        results.better_method = 'Welch (simplified)';
        results.fc = fc;
        results.snr = snr_estimate;
        results.total_time = total_time;
    catch
        total_time = toc(sim_start_time);
        end_time_str = datestr(now, 'yyyy-mm-dd HH:MM:SS');
        fprintf('简化方法也失败，请检查MATLAB工具箱\n');
        fprintf('\n========================================\n');
        fprintf('仿真失败总结\n');
        fprintf('========================================\n');
        fprintf('开始时间: %s\n', start_time_str);
        fprintf('结束时间: %s\n', end_time_str);
        fprintf('总仿真时间: %.2f 秒 (%.2f 分钟)\n', total_time, total_time/60);
        fprintf('========================================\n\n');
        B_welch = [];
        B_ar = [];
        results = [];
    end
end

% 如果成功完成，输出最终总结
if exist('B_welch', 'var') && ~isempty(B_welch) && exist('end_time_str', 'var')
    fprintf('========================================\n');
    fprintf('带宽估计完成！\n');
    fprintf('========================================\n');
    fprintf('开始时间: %s\n', start_time_str);
    fprintf('结束时间: %s\n', end_time_str);
    fprintf('总仿真时间: %.2f 秒 (%.2f 分钟)\n', total_time, total_time/60);
    fprintf('========================================\n\n');
end

end

