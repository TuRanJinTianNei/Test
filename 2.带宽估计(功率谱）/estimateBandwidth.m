function [B_welch, B_ar, B_ideal, results] = estimateBandwidth(sig_processed, fs, ...
    carrier_count, subcarrier_spacing, varargin)
%===============================================================================
% estimateBandwidth.m - 完整的带宽估计流程
% 
% 功能说明:
%   完整的带宽估计流程，包括：
%   1. 调用estimatePSD.m生成两个功率谱（Welch和AR）
%   2. 使用estimateBandwidthMethods.m的方法对两个功率谱进行带宽估计
%   3. 对比两种PSD估计方法的带宽估计结果
%
% 依赖关系:
%   - estimatePSD.m: 功率谱密度估计
%   - estimateBandwidthMethods.m: 带宽估计方法
%
% 输入参数:
%   sig_processed      - 输入信号（实数或复数信号，行向量或列向量）
%                        可以是：
%                        - signalGenerate.m生成的Rx_data（推荐）
%                        - signalGenerate.m生成的Tx_data
%                        - 其他已处理的OFDM基带信号
%   fs                 - 采样频率（Hz）
%                        对于signalGenerate.m：fs = 3.84e6 Hz（LTE 3MHz标准）
%   carrier_count      - 子载波数（用于计算理论带宽）
%   subcarrier_spacing - 子载波间隔（Hz）
%                        对于signalGenerate.m：subcarrier_spacing = 15e3 Hz
%   varargin           - 可选参数（名称-值对）:
%                        'snr'            - 信噪比（dB），用于显示，默认20
%                        'fc'             - 载波频率（Hz），默认0（基带信号）
%                        'plot'           - 是否绘制结果图，默认true
%                        'method'         - 带宽估计方法，可选：
%                                          'threshold' - 阈值法（默认）
%                                          'energy'    - 能量百分比法
%                                          'peak_drop' - 峰值下降法
%                                          'rms'       - RMS带宽法
%                        'threshold'      - 阈值法阈值（dB），默认-3
%                        'energy_percent' - 能量百分比（0-1），默认0.9（90%）
%                        'peak_drop_db'   - 峰值下降dB值，默认3
%                        'welch_window'   - Welch算法窗长度，默认100
%                        'welch_overlap'  - Welch算法重叠样本数，默认55
%                        'welch_nfft'     - Welch算法FFT点数，默认8192
%                        'ar_criterion'   - AR模型阶数选择准则，默认'AIC'
%
% 输出参数:
%   B_welch            - Welch算法PSD的带宽估计结果（Hz）
%                        如果method='all'，返回结构体，包含所有方法的估计值
%                        否则返回标量
%   B_ar               - AR模型PSD的带宽估计结果（Hz）
%                        如果method='all'，返回结构体，包含所有方法的估计值
%                        否则返回标量
%   B_ideal             - 理论带宽（Hz）= carrier_count * subcarrier_spacing
%   results             - 详细结果结构体，包含：
%                        .welch
%                                          'integral'  - 积分法
%                                          'all'       - 使用所有方法 - Welch算法PSD的带宽估计结果
%                        .ar    - AR模型PSD的带宽估计结果
%                        .psd   - PSD估计结果（Pxx_welch, f_welch, Pxx_ar, f_ar）
%
% 使用示例:
%   % 示例1：使用接收信号进行带宽估计
%   signalGenerate;  % 生成Tx_data和Rx_data
%   [B_w, B_a, B_i, results] = estimateBandwidth(Rx_data, fs, ...
%       carrier_count, subcarrier_spacing, 'snr', targetSNRdB, 'plot', true);
%
%   % 示例2：使用所有方法进行带宽估计
%   [B_w, B_a, B_i, results] = estimateBandwidth(Rx_data, fs, ...
%       carrier_count, subcarrier_spacing, 'method', 'all', 'plot', true);
%
%   % 示例3：使用能量百分比法
%   [B_w, B_a, B_i, results] = estimateBandwidth(Rx_data, fs, ...
%       carrier_count, subcarrier_spacing, 'method', 'energy', ...
%       'energy_percent', 0.9);
%
% 创建日期: 2025.12.23
%===============================================================================

%**************************************************************************
% 独立运行模式：如果没有输入参数，自动检测或生成信号
%**************************************************************************
if nargin == 0
    fprintf('========================================\n');
    fprintf('带宽估计程序（独立运行模式）\n');
    fprintf('========================================\n\n');
    
    % 检查工作空间中是否已有信号数据和系统参数
    if exist('Rx_data', 'var') && exist('fs', 'var') && ...
       exist('carrier_count', 'var') && exist('subcarrier_spacing', 'var')
        fprintf('检测到工作空间中已有信号数据和系统参数，直接使用...\n');
        fprintf('  - 使用变量: Rx_data\n');
        fprintf('  - 采样频率: %.2f MHz\n', fs/1e6);
        fprintf('  - 子载波数: %d\n', carrier_count);
        fprintf('  - 子载波间隔: %.1f kHz\n', subcarrier_spacing/1e3);
        if exist('targetSNRdB', 'var')
            fprintf('  - SNR: %.1f dB\n\n', targetSNRdB);
        else
            fprintf('  - SNR: 20 dB（默认值，未找到targetSNRdB）\n\n');
        end
        sig_processed = Rx_data;
    elseif exist('Tx_data', 'var') && exist('fs', 'var') && ...
           exist('carrier_count', 'var') && exist('subcarrier_spacing', 'var')
        fprintf('检测到工作空间中已有发送信号和系统参数，使用Tx_data...\n');
        fprintf('  - 使用变量: Tx_data\n');
        fprintf('  - 采样频率: %.2f MHz\n', fs/1e6);
        fprintf('  - 子载波数: %d\n', carrier_count);
        fprintf('  - 子载波间隔: %.1f kHz\n', subcarrier_spacing/1e3);
        if exist('targetSNRdB', 'var')
            fprintf('  - SNR: %.1f dB\n\n', targetSNRdB);
        else
            fprintf('  - SNR: 20 dB（默认值，未找到targetSNRdB）\n\n');
        end
        sig_processed = Tx_data;
    else
        fprintf('未找到信号数据或系统参数，正在调用signalGenerate.m生成信号...\n\n');
        % 调用signalGenerate.m生成信号
        signalGenerate;
        
        % 检查是否生成了接收信号和系统参数
        if exist('Rx_data', 'var') && exist('fs', 'var') && ...
           exist('carrier_count', 'var') && exist('subcarrier_spacing', 'var')
            fprintf('\n使用生成的接收信号（Rx_data）进行带宽估计...\n');
            sig_processed = Rx_data;
            fprintf('  采样频率: %.2f MHz\n', fs/1e6);
            fprintf('  子载波数: %d\n', carrier_count);
            fprintf('  子载波间隔: %.1f kHz\n', subcarrier_spacing/1e3);
            if exist('targetSNRdB', 'var')
                fprintf('  SNR: %.1f dB\n\n', targetSNRdB);
            else
                fprintf('  SNR: 20 dB（默认值）\n\n');
            end
        elseif exist('Tx_data', 'var') && exist('fs', 'var') && ...
               exist('carrier_count', 'var') && exist('subcarrier_spacing', 'var')
            fprintf('\n使用生成的发送信号（Tx_data）进行带宽估计...\n');
            sig_processed = Tx_data;
            fprintf('  采样频率: %.2f MHz\n', fs/1e6);
            fprintf('  子载波数: %d\n', carrier_count);
            fprintf('  子载波间隔: %.1f kHz\n', subcarrier_spacing/1e3);
            if exist('targetSNRdB', 'var')
                fprintf('  SNR: %.1f dB\n\n', targetSNRdB);
            else
                fprintf('  SNR: 20 dB（默认值）\n\n');
            end
        else
            error('错误：signalGenerate.m未能生成信号数据或缺少系统参数（fs, carrier_count, subcarrier_spacing）');
        end
    end
    
    % 独立运行模式：设置默认参数
    % 如果工作空间中有targetSNRdB，使用它；否则使用默认值20
    if exist('targetSNRdB', 'var')
        varargin = {'plot', true, 'snr', targetSNRdB};
    else
        varargin = {'plot', true, 'snr', 20};
    end
    fprintf('独立运行模式：默认启用绘图\n\n');
end

% 输入验证（函数调用模式）
if nargin > 0 && nargin < 4
    error('错误：函数调用模式需要至少4个输入参数：sig_processed, fs, carrier_count, subcarrier_spacing\n或者无参数运行以进入独立运行模式');
end

% 解析可选参数
p = inputParser;
% 如果工作空间中有targetSNRdB，使用它作为默认值；否则使用20
if exist('targetSNRdB', 'var')
    default_snr = targetSNRdB;
else
    default_snr = 20;
end
addParameter(p, 'snr', default_snr, @isnumeric);
addParameter(p, 'fc', 0, @isnumeric);
addParameter(p, 'plot', true, @islogical);
addParameter(p, 'method', 'threshold', @(x) ismember(x, ...
    {'threshold', 'energy', 'peak_drop', 'rms', 'integral', 'all'}));
addParameter(p, 'threshold', -3, @isnumeric);
addParameter(p, 'energy_percent', 0.9, @(x) isnumeric(x) && x > 0 && x <= 1);
addParameter(p, 'peak_drop_db', 3, @isnumeric);
addParameter(p, 'welch_window', 100, @isnumeric);
addParameter(p, 'welch_overlap', 55, @isnumeric);
addParameter(p, 'welch_nfft', 8192, @isnumeric);
addParameter(p, 'ar_criterion', 'AIC', @(x) ischar(x) && ...
    (strcmpi(x, 'AIC') || strcmpi(x, 'FPE')));
parse(p, varargin{:});

snr = p.Results.snr;
fc = p.Results.fc;
plot_flag = p.Results.plot;
method = p.Results.method;
threshold = p.Results.threshold;
energy_percent = p.Results.energy_percent;
peak_drop_db = p.Results.peak_drop_db;
welch_window = p.Results.welch_window;
welch_overlap = p.Results.welch_overlap;
welch_nfft = p.Results.welch_nfft;
ar_criterion = p.Results.ar_criterion;

% 计算理论带宽
B_ideal = carrier_count * subcarrier_spacing;

fprintf('========================================\n');
fprintf('带宽估计流程\n');
fprintf('========================================\n');
fprintf('输入参数：\n');
fprintf('  - 信号长度: %d 样本\n', length(sig_processed));
fprintf('  - 采样频率: %.2f MHz\n', fs/1e6);
fprintf('  - SNR: %.1f dB\n', snr);
fprintf('  - 载波频率: %.3f MHz\n', fc/1e6);
fprintf('  - 子载波数: %d\n', carrier_count);
fprintf('  - 子载波间隔: %.1f kHz\n', subcarrier_spacing/1e3);
fprintf('  - 理论带宽: %.3f MHz (%.0f kHz)\n', B_ideal/1e6, B_ideal/1e3);
fprintf('  - 带宽估计方法: %s\n', method);
fprintf('========================================\n\n');

%===============================================================================
% 步骤1: 估计功率谱密度（调用estimatePSD）
%===============================================================================
fprintf('步骤1: 估计功率谱密度...\n');
fprintf('----------------------------------------\n');

[Pxx_welch, f_welch, Pxx_ar, f_ar] = estimatePSD(sig_processed, fs, snr, ...
    'plot', false, ...
    'welch_window', welch_window, ...
    'welch_overlap', welch_overlap, ...
    'welch_nfft', welch_nfft, ...
    'ar_criterion', ar_criterion);

fprintf('  - Welch算法：频率点数 = %d，频率范围 = %.3f MHz 到 %.3f MHz\n', ...
    length(f_welch), min(f_welch)/1e6, max(f_welch)/1e6);
fprintf('  - AR模型法：频率点数 = %d，频率范围 = %.3f MHz 到 %.3f MHz\n', ...
    length(f_ar), min(f_ar)/1e6, max(f_ar)/1e6);
fprintf('  完成！\n\n');

%===============================================================================
% 步骤2: 对Welch算法PSD进行带宽估计（调用estimateBandwidthMethods）
%===============================================================================
fprintf('步骤2: 对Welch算法PSD进行带宽估计...\n');
fprintf('----------------------------------------\n');

[B_welch, results_welch] = estimateBandwidthMethods(Pxx_welch, f_welch, ...
    'method', method, ...
    'threshold', threshold, ...
    'energy_percent', energy_percent, ...
    'peak_drop_db', peak_drop_db, ...
    'plot', false);

printBandwidthResults('Welch算法', B_welch, results_welch);
fprintf('  完成！\n\n');

%===============================================================================
% 步骤3: 对AR模型PSD进行带宽估计（调用estimateBandwidthMethods）
%===============================================================================
fprintf('步骤3: 对AR模型PSD进行带宽估计...\n');
fprintf('----------------------------------------\n');

[B_ar, results_ar] = estimateBandwidthMethods(Pxx_ar, f_ar, ...
    'method', method, ...
    'threshold', threshold, ...
    'energy_percent', energy_percent, ...
    'peak_drop_db', peak_drop_db, ...
    'plot', false);

printBandwidthResults('AR模型法', B_ar, results_ar);
fprintf('  完成！\n\n');

%===============================================================================
% 步骤4: 结果汇总
%===============================================================================
fprintf('========================================\n');
fprintf('带宽估计结果汇总\n');
fprintf('========================================\n');
fprintf('理论带宽: %.3f MHz (%.0f kHz)\n', B_ideal/1e6, B_ideal/1e3);
fprintf('\n');

if isstruct(B_welch)
    % 使用所有方法
    fprintf('Welch算法估计结果（所有方法）：\n');
    methods = {'threshold', 'energy', 'peak_drop', 'rms', 'integral'};
    method_names = {'阈值法', '能量法', '峰值下降法', 'RMS法', '积分法'};
    B_welch_values = [B_welch.threshold, B_welch.energy, ...
        B_welch.peak_drop, B_welch.rms, B_welch.integral];
    for i = 1:length(methods)
        error_abs = abs(B_welch_values(i) - B_ideal);
        error_rel = error_abs / B_ideal * 100;
        fprintf('  %s: %.3f MHz, 误差: %.3f MHz (%.2f%%)\n', ...
            method_names{i}, B_welch_values(i)/1e6, error_abs/1e6, error_rel);
    end
    fprintf('  平均带宽: %.3f MHz, 标准差: %.3f MHz\n', ...
        B_welch.mean/1e6, B_welch.std/1e6);
else
    % 单一方法
    error_abs_w = abs(B_welch - B_ideal);
    error_rel_w = error_abs_w / B_ideal * 100;
    fprintf('Welch算法估计结果：\n');
    fprintf('  估计带宽: %.3f MHz\n', B_welch/1e6);
    fprintf('  误差: %.3f MHz (%.2f%%)\n', error_abs_w/1e6, error_rel_w);
end

fprintf('\n');

if isstruct(B_ar)
    % 使用所有方法
    fprintf('AR模型法估计结果（所有方法）：\n');
    methods = {'threshold', 'energy', 'peak_drop', 'rms', 'integral'};
    method_names = {'阈值法', '能量法', '峰值下降法', 'RMS法', '积分法'};
    B_ar_values = [B_ar.threshold, B_ar.energy, ...
        B_ar.peak_drop, B_ar.rms, B_ar.integral];
    for i = 1:length(methods)
        error_abs = abs(B_ar_values(i) - B_ideal);
        error_rel = error_abs / B_ideal * 100;
        fprintf('  %s: %.3f MHz, 误差: %.3f MHz (%.2f%%)\n', ...
            method_names{i}, B_ar_values(i)/1e6, error_abs/1e6, error_rel);
    end
    fprintf('  平均带宽: %.3f MHz, 标准差: %.3f MHz\n', ...
        B_ar.mean/1e6, B_ar.std/1e6);
else
    % 单一方法
    error_abs_a = abs(B_ar - B_ideal);
    error_rel_a = error_abs_a / B_ideal * 100;
    fprintf('AR模型法估计结果：\n');
    fprintf('  估计带宽: %.3f MHz\n', B_ar/1e6);
    fprintf('  误差: %.3f MHz (%.2f%%)\n', error_abs_a/1e6, error_rel_a);
end

fprintf('========================================\n\n');

%===============================================================================
% 步骤5: 组织输出结果
%===============================================================================
results = struct();
results.welch = results_welch;
results.ar = results_ar;
results.psd = struct();
results.psd.Pxx_welch = Pxx_welch;
results.psd.f_welch = f_welch;
results.psd.Pxx_ar = Pxx_ar;
results.psd.f_ar = f_ar;
results.B_ideal = B_ideal;
results.method = method;

%===============================================================================
% 步骤6: 绘制对比图（可选）
%===============================================================================
if plot_flag
    fprintf('正在绘制对比图...\n');
    
    figure('Name', '带宽估计结果对比（Welch vs AR）', 'Position', [100, 100, 1600, 1000]);
    
    % 子图1：Welch算法PSD和带宽估计
    subplot(2, 3, 1);
    plot(f_welch/1e6, Pxx_welch, 'b-', 'LineWidth', 1.5);
    hold on;
    % 检查是否是'all'模式的结果（包含多个方法）
    if isstruct(results_welch) && isfield(results_welch, 'threshold')
        % 所有方法的结果
        methods = {'threshold', 'energy', 'peak_drop', 'rms', 'integral'};
        colors = lines(length(methods));
        for i = 1:length(methods)
            if isfield(results_welch, methods{i}) && isstruct(results_welch.(methods{i}))
                r = results_welch.(methods{i});
                plot([r.lower_freq/1e6, r.lower_freq/1e6], ylim, '--', ...
                    'Color', colors(i,:), 'LineWidth', 1.5);
                plot([r.upper_freq/1e6, r.upper_freq/1e6], '--', ...
                    'Color', colors(i,:), 'LineWidth', 1.5);
            end
        end
    elseif isstruct(results_welch) && isfield(results_welch, 'lower_freq')
        % 单一方法的结果
        plot([results_welch.lower_freq/1e6, results_welch.lower_freq/1e6], ...
            ylim, 'r--', 'LineWidth', 2, 'DisplayName', '下边界');
        plot([results_welch.upper_freq/1e6, results_welch.upper_freq/1e6], ...
            ylim, 'r--', 'LineWidth', 2, 'DisplayName', '上边界');
    end
    grid on;
    xlabel('频率 (MHz)');
    ylabel('PSD (dB)');
    title(sprintf('Welch算法PSD和带宽估计\n(SNR=%.1f dB)', snr));
    legend('Location', 'best');
    hold off;
    
    % 子图2：AR模型PSD和带宽估计
    subplot(2, 3, 2);
    plot(f_ar/1e6, Pxx_ar, 'r-', 'LineWidth', 1.5);
    hold on;
    % 检查是否是'all'模式的结果（包含多个方法）
    if isstruct(results_ar) && isfield(results_ar, 'threshold')
        % 所有方法的结果
        methods = {'threshold', 'energy', 'peak_drop', 'rms', 'integral'};
        colors = lines(length(methods));
        for i = 1:length(methods)
            if isfield(results_ar, methods{i}) && isstruct(results_ar.(methods{i}))
                r = results_ar.(methods{i});
                plot([r.lower_freq/1e6, r.lower_freq/1e6], ylim, '--', ...
                    'Color', colors(i,:), 'LineWidth', 1.5);
                plot([r.upper_freq/1e6, r.upper_freq/1e6], '--', ...
                    'Color', colors(i,:), 'LineWidth', 1.5);
            end
        end
    elseif isstruct(results_ar) && isfield(results_ar, 'lower_freq')
        % 单一方法的结果
        plot([results_ar.lower_freq/1e6, results_ar.lower_freq/1e6], ...
            ylim, 'r--', 'LineWidth', 2, 'DisplayName', '下边界');
        plot([results_ar.upper_freq/1e6, results_ar.upper_freq/1e6], ...
            ylim, 'r--', 'LineWidth', 2, 'DisplayName', '上边界');
    end
    grid on;
    xlabel('频率 (MHz)');
    ylabel('PSD (dB)');
    title(sprintf('AR模型PSD和带宽估计\n(SNR=%.1f dB)', snr));
    legend('Location', 'best');
    hold off;
    
    % 子图3：Welch和AR PSD对比
    subplot(2, 3, 3);
    plot(f_welch/1e6, Pxx_welch, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Welch');
    hold on;
    plot(f_ar/1e6, Pxx_ar, 'r-', 'LineWidth', 1.5, 'DisplayName', 'AR模型');
    grid on;
    xlabel('频率 (MHz)');
    ylabel('PSD (dB)');
    title('PSD估计方法对比');
    legend('Location', 'best');
    hold off;
    
    % 子图4：带宽估计值对比
    subplot(2, 3, 4);
    if isstruct(B_welch) && isstruct(B_ar)
        % 所有方法的结果
        methods = {'threshold', 'energy', 'peak_drop', 'rms', 'integral'};
        method_names = {'阈值法', '能量法', '峰值下降法', 'RMS法', '积分法'};
        B_welch_values = [B_welch.threshold, B_welch.energy, ...
            B_welch.peak_drop, B_welch.rms, B_welch.integral];
        B_ar_values = [B_ar.threshold, B_ar.energy, ...
            B_ar.peak_drop, B_ar.rms, B_ar.integral];
        
        bar([B_welch_values/1e6; B_ar_values/1e6]');
        set(gca, 'XTickLabel', method_names, 'XTickLabelRotation', 45);
        ylabel('带宽 (MHz)');
        title('各方法带宽估计值对比');
        legend('Welch算法', 'AR模型法', 'Location', 'best');
        hold on;
        plot([0.5, length(methods)+0.5], [B_ideal/1e6, B_ideal/1e6], ...
            'k--', 'LineWidth', 2, 'DisplayName', '理论值');
        legend('Location', 'best');
        grid on;
        hold off;
    else
        % 单一方法的结果
        bar([B_welch/1e6; B_ar/1e6]);
        set(gca, 'XTickLabel', {'Welch算法', 'AR模型法'});
        ylabel('带宽 (MHz)');
        title('带宽估计值对比');
        hold on;
        plot([0.5, 2.5], [B_ideal/1e6, B_ideal/1e6], ...
            'r--', 'LineWidth', 2, 'DisplayName', '理论值');
        legend('Location', 'best');
        grid on;
        hold off;
    end
    
    % 子图5：误差对比
    subplot(2, 3, 5);
    if isstruct(B_welch) && isstruct(B_ar)
        % 所有方法的结果
        method_names = {'阈值法', '能量法', '峰值下降法', 'RMS法', '积分法'};
        B_welch_values = [B_welch.threshold, B_welch.energy, ...
            B_welch.peak_drop, B_welch.rms, B_welch.integral];
        B_ar_values = [B_ar.threshold, B_ar.energy, ...
            B_ar.peak_drop, B_ar.rms, B_ar.integral];
        
        errors_welch = abs(B_welch_values - B_ideal) / B_ideal * 100;
        errors_ar = abs(B_ar_values - B_ideal) / B_ideal * 100;
        
        bar([errors_welch; errors_ar]');
        set(gca, 'XTickLabel', method_names, 'XTickLabelRotation', 45);
        ylabel('相对误差 (%)');
        title('各方法估计误差对比');
        legend('Welch算法', 'AR模型法', 'Location', 'best');
        grid on;
    else
        % 单一方法的结果
        error_welch = abs(B_welch - B_ideal) / B_ideal * 100;
        error_ar = abs(B_ar - B_ideal) / B_ideal * 100;
        
        bar([error_welch; error_ar]);
        set(gca, 'XTickLabel', {'Welch算法', 'AR模型法'});
        ylabel('相对误差 (%)');
        title('估计误差对比');
        grid on;
    end
    
    % 子图6：结果汇总文本
    subplot(2, 3, 6);
    axis off;
    text(0.1, 0.9, '带宽估计结果汇总', 'FontSize', 14, ...
        'FontWeight', 'bold', 'Units', 'normalized');
    text(0.1, 0.8, sprintf('理论带宽: %.3f MHz', B_ideal/1e6), ...
        'FontSize', 12, 'Units', 'normalized');
    
    if isstruct(B_welch)
        text(0.1, 0.7, sprintf('Welch平均: %.3f MHz', B_welch.mean/1e6), ...
            'FontSize', 11, 'Units', 'normalized', 'Color', 'b');
        text(0.1, 0.65, sprintf('AR平均: %.3f MHz', B_ar.mean/1e6), ...
            'FontSize', 11, 'Units', 'normalized', 'Color', 'r');
    else
        text(0.1, 0.7, sprintf('Welch估计: %.3f MHz', B_welch/1e6), ...
            'FontSize', 11, 'Units', 'normalized', 'Color', 'b');
        text(0.1, 0.65, sprintf('AR估计: %.3f MHz', B_ar/1e6), ...
            'FontSize', 11, 'Units', 'normalized', 'Color', 'r');
    end
    
    fprintf('图形绘制完成！\n\n');
end

fprintf('========================================\n');
fprintf('带宽估计完成！\n');
fprintf('========================================\n');

end

%===============================================================================
% 辅助函数：打印带宽估计结果
%===============================================================================
function printBandwidthResults(method_name, B_result, results_struct)
% 打印带宽估计结果的辅助函数
% 输入参数:
%   method_name  - 方法名称（如'Welch算法'或'AR模型法'）
%   B_result     - 带宽估计结果（标量或结构体）
%   results_struct - 详细结果结构体

if isstruct(B_result) && isfield(B_result, 'threshold')
    % 所有方法的结果
    fprintf('  %s - 各方法估计结果：\n', method_name);
    fprintf('    阈值法: %.3f MHz\n', B_result.threshold/1e6);
    fprintf('    能量法: %.3f MHz\n', B_result.energy/1e6);
    fprintf('    峰值下降法: %.3f MHz\n', B_result.peak_drop/1e6);
    fprintf('    RMS法: %.3f MHz\n', B_result.rms/1e6);
    fprintf('    积分法: %.3f MHz\n', B_result.integral/1e6);
    fprintf('    平均带宽: %.3f MHz (标准差: %.3f MHz)\n', ...
        B_result.mean/1e6, B_result.std/1e6);
else
    % 单一方法的结果
    fprintf('  %s估计带宽: %.3f MHz\n', method_name, B_result/1e6);
    fprintf('    下边界: %.3f MHz, 上边界: %.3f MHz\n', ...
        results_struct.lower_freq/1e6, results_struct.upper_freq/1e6);
end

end

