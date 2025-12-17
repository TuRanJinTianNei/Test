function [B_estimated, results] = estimateBandwidthFromPSD(Pxx, f, varargin)
%===============================================================================
% estimateBandwidthFromPSD.m - 从功率谱密度估计带宽
% 
% 功能说明:
%   提供多种方法从功率谱密度（PSD）估计信号带宽
%   支持的方法：
%   1. 阈值法（固定阈值或自适应阈值）
%   2. 能量百分比法（如90%能量带宽）
%   3. 峰值下降法（从峰值下降固定dB）
%   4. RMS带宽法（二阶矩法）
%   5. 积分法（累积功率法）
%
% 输入参数:
%   Pxx         - 功率谱密度（dB，归一化，行向量或列向量）
%   f           - 对应的频率向量（Hz，行向量或列向量）
%   varargin    - 可选参数（名称-值对）:
%                 'method'      - 估计方法，可选：
%                                 'threshold' - 阈值法（默认）
%                                 'energy'    - 能量百分比法
%                                 'peak_drop' - 峰值下降法
%                                 'rms'       - RMS带宽法
%                                 'integral'  - 积分法
%                                 'all'       - 使用所有方法
%                 'threshold'   - 阈值法阈值（dB），默认-3
%                 'energy_percent' - 能量百分比（0-1），默认0.9（90%）
%                 'peak_drop_db'  - 峰值下降dB值，默认3
%                 'plot'         - 是否绘制结果图，默认false
%                 'peak_freq'    - 峰值频率（Hz），如果未提供则自动检测
%
% 输出参数:
%   B_estimated - 估计的带宽（Hz）
%                  - 如果method='all'，返回结构体，包含所有方法的估计值
%                  - 否则返回标量
%   results     - 详细结果结构体，包含：
%                 .method - 使用的方法
%                 .bandwidth - 估计带宽（Hz）
%                 .lower_freq - 下边界频率（Hz）
%                 .upper_freq - 上边界频率（Hz）
%                 .peak_freq - 峰值频率（Hz）
%                 .peak_power - 峰值功率（dB）
%
% 使用示例:
%   % 方法1：阈值法
%   [B, results] = estimateBandwidthFromPSD(Pxx, f, 'method', 'threshold', ...
%       'threshold', -3);
%
%   % 方法2：能量百分比法
%   [B, results] = estimateBandwidthFromPSD(Pxx, f, 'method', 'energy', ...
%       'energy_percent', 0.9);
%
%   % 方法3：使用所有方法
%   [B_all, results_all] = estimateBandwidthFromPSD(Pxx, f, 'method', 'all');
%
% 创建日期: 2025.12.10
%===============================================================================

% 解析可选参数
p = inputParser;
addParameter(p, 'method', 'threshold', @(x) ismember(x, ...
    {'threshold', 'energy', 'peak_drop', 'rms', 'integral', 'all'}));
addParameter(p, 'threshold', -3, @isnumeric);
addParameter(p, 'energy_percent', 0.9, @(x) isnumeric(x) && x > 0 && x <= 1);
addParameter(p, 'peak_drop_db', 3, @isnumeric);
addParameter(p, 'plot', false, @islogical);
addParameter(p, 'peak_freq', [], @isnumeric);
parse(p, varargin{:});

method = p.Results.method;
threshold = p.Results.threshold;
energy_percent = p.Results.energy_percent;
peak_drop_db = p.Results.peak_drop_db;
plot_flag = p.Results.plot;
peak_freq = p.Results.peak_freq;

% 确保输入为列向量
Pxx = Pxx(:);
f = f(:);

% 检查输入有效性
if length(Pxx) ~= length(f)
    error('Pxx和f的长度必须相同');
end

% 转换为线性功率（从dB）
Pxx_linear = 10.^(Pxx / 10);

% 检测峰值频率和峰值功率
if isempty(peak_freq)
    [peak_power, peak_idx] = max(Pxx);
    peak_freq = f(peak_idx);
else
    [~, peak_idx] = min(abs(f - peak_freq));
    peak_power = Pxx(peak_idx);
end

% 根据方法选择执行
if strcmpi(method, 'all')
    % 使用所有方法
    methods = {'threshold', 'energy', 'peak_drop', 'rms', 'integral'};
    B_estimated = struct();
    results = struct();
    
    for i = 1:length(methods)
        [B_temp, results_temp] = estimateBandwidthFromPSD(Pxx, f, ...
            'method', methods{i}, ...
            'threshold', threshold, ...
            'energy_percent', energy_percent, ...
            'peak_drop_db', peak_drop_db, ...
            'peak_freq', peak_freq, ...
            'plot', false);
        
        B_estimated.(methods{i}) = B_temp;
        results.(methods{i}) = results_temp;
    end
    
    % 计算平均带宽（排除异常值）
    B_values = [B_estimated.threshold, B_estimated.energy, ...
        B_estimated.peak_drop, B_estimated.rms, B_estimated.integral];
    B_mean = mean(B_values);
    B_std = std(B_values);
    
    % 移除超出2倍标准差的值后重新计算平均
    valid_idx = abs(B_values - B_mean) <= 2 * B_std;
    if sum(valid_idx) > 0
        B_estimated.mean = mean(B_values(valid_idx));
        B_estimated.std = std(B_values(valid_idx));
    else
        B_estimated.mean = B_mean;
        B_estimated.std = B_std;
    end
    
    if plot_flag
        plotBandwidthResults(Pxx, f, results, B_estimated);
    end
    
else
    % 使用单一方法
    switch lower(method)
        case 'threshold'
            [B_estimated, results] = method_threshold(Pxx, f, threshold, ...
                peak_freq, peak_power, peak_idx);
            
        case 'energy'
            [B_estimated, results] = method_energy(Pxx_linear, f, ...
                energy_percent, peak_freq, peak_power, peak_idx);
            
        case 'peak_drop'
            [B_estimated, results] = method_peak_drop(Pxx, f, peak_drop_db, ...
                peak_freq, peak_power, peak_idx);
            
        case 'rms'
            [B_estimated, results] = method_rms(Pxx_linear, f, ...
                peak_freq, peak_power, peak_idx);
            
        case 'integral'
            [B_estimated, results] = method_integral(Pxx_linear, f, ...
                peak_freq, peak_power, peak_idx);
            
        otherwise
            error('未知的方法: %s', method);
    end
    
    results.method = method;
    
    if plot_flag
        plotBandwidthResults(Pxx, f, results, B_estimated);
    end
end

end

%===============================================================================
% 方法1：阈值法
%===============================================================================
function [B, results] = method_threshold(Pxx, f, threshold, ...
    peak_freq, peak_power, peak_idx)

% 在峰值左侧（低频侧）找阈值点
left_side = Pxx(1:peak_idx);
left_indices = find(left_side <= threshold);
if ~isempty(left_indices)
    left_idx = left_indices(end);
    lower_freq = f(left_idx);
else
    lower_freq = f(1);
end

% 在峰值右侧（高频侧）找阈值点
right_side = Pxx(peak_idx:end);
right_indices = find(right_side <= threshold);
if ~isempty(right_indices)
    right_idx = peak_idx + right_indices(1) - 1;
    upper_freq = f(right_idx);
else
    upper_freq = f(end);
end

B = abs(upper_freq - lower_freq);

results = struct();
results.bandwidth = B;
results.lower_freq = lower_freq;
results.upper_freq = upper_freq;
results.peak_freq = peak_freq;
results.peak_power = peak_power;
results.threshold = threshold;

end

%===============================================================================
% 方法2：能量百分比法（90%能量带宽等）
%===============================================================================
function [B, results] = method_energy(Pxx_linear, f, energy_percent, ...
    peak_freq, peak_power, peak_idx)

% 计算总能量
total_energy = trapz(f, Pxx_linear);

% 目标能量
target_energy = energy_percent * total_energy;

% 从峰值向两侧扩展，找到包含目标能量的频率范围
% 方法：从峰值开始，向两侧累积能量

% 初始化
current_energy = Pxx_linear(peak_idx) * (f(2) - f(1));  % 峰值处的能量
lower_idx = peak_idx;
upper_idx = peak_idx;

% 向两侧扩展，直到累积能量达到目标
while current_energy < target_energy && ...
        (lower_idx > 1 || upper_idx < length(f))
    
    % 计算向两侧扩展的能量增量
    if lower_idx > 1
        energy_left = Pxx_linear(lower_idx - 1) * (f(lower_idx) - f(lower_idx - 1));
    else
        energy_left = 0;
    end
    
    if upper_idx < length(f)
        energy_right = Pxx_linear(upper_idx + 1) * (f(upper_idx + 1) - f(upper_idx));
    else
        energy_right = 0;
    end
    
    % 选择能量增量更大的一侧扩展
    if energy_left > energy_right && lower_idx > 1
        lower_idx = lower_idx - 1;
        current_energy = current_energy + energy_left;
    elseif upper_idx < length(f)
        upper_idx = upper_idx + 1;
        current_energy = current_energy + energy_right;
    else
        break;
    end
end

lower_freq = f(lower_idx);
upper_freq = f(upper_idx);
B = abs(upper_freq - lower_freq);

results = struct();
results.bandwidth = B;
results.lower_freq = lower_freq;
results.upper_freq = upper_freq;
results.peak_freq = peak_freq;
results.peak_power = peak_power;
results.energy_percent = energy_percent;
results.actual_energy_percent = current_energy / total_energy;

end

%===============================================================================
% 方法3：峰值下降法
%===============================================================================
function [B, results] = method_peak_drop(Pxx, f, peak_drop_db, ...
    peak_freq, peak_power, peak_idx)

threshold = peak_power - peak_drop_db;

% 在峰值左侧（低频侧）找阈值点
left_side = Pxx(1:peak_idx);
left_indices = find(left_side <= threshold);
if ~isempty(left_indices)
    left_idx = left_indices(end);
    lower_freq = f(left_idx);
else
    lower_freq = f(1);
end

% 在峰值右侧（高频侧）找阈值点
right_side = Pxx(peak_idx:end);
right_indices = find(right_side <= threshold);
if ~isempty(right_indices)
    right_idx = peak_idx + right_indices(1) - 1;
    upper_freq = f(right_idx);
else
    upper_freq = f(end);
end

B = abs(upper_freq - lower_freq);

results = struct();
results.bandwidth = B;
results.lower_freq = lower_freq;
results.upper_freq = upper_freq;
results.peak_freq = peak_freq;
results.peak_power = peak_power;
results.peak_drop_db = peak_drop_db;
results.threshold = threshold;

end

%===============================================================================
% 方法4：RMS带宽法（二阶矩法）
%===============================================================================
function [B, results] = method_rms(Pxx_linear, f, ...
    peak_freq, peak_power, peak_idx)

% 计算归一化功率谱密度
Pxx_norm = Pxx_linear / trapz(f, Pxx_linear);

% 计算一阶矩（中心频率）
f_center = trapz(f, f .* Pxx_norm);

% 计算二阶矩（RMS带宽）
f_squared = trapz(f, (f - f_center).^2 .* Pxx_norm);
B_rms = 2 * sqrt(f_squared);  % RMS带宽

% 计算带宽边界（中心频率 ± RMS带宽/2）
lower_freq = f_center - B_rms / 2;
upper_freq = f_center + B_rms / 2;

% 限制在有效频率范围内
lower_freq = max(lower_freq, f(1));
upper_freq = min(upper_freq, f(end));

B = B_rms;

results = struct();
results.bandwidth = B;
results.lower_freq = lower_freq;
results.upper_freq = upper_freq;
results.peak_freq = peak_freq;
results.peak_power = peak_power;
results.center_freq = f_center;
results.rms_bandwidth = B_rms;

end

%===============================================================================
% 方法5：积分法（累积功率法）
%===============================================================================
function [B, results] = method_integral(Pxx_linear, f, ...
    peak_freq, peak_power, peak_idx)

% 计算累积功率
cumulative_power = cumtrapz(f, Pxx_linear);
total_power = cumulative_power(end);

% 找到峰值位置对应的累积功率
peak_cumulative = cumulative_power(peak_idx);

% 定义带宽边界（例如峰值两侧各占一定百分比）
% 这里使用峰值两侧各占45%的功率（共90%）
power_lower = peak_cumulative - 0.45 * total_power;
power_upper = peak_cumulative + 0.45 * total_power;

% 限制在有效范围内
power_lower = max(power_lower, cumulative_power(1));
power_upper = min(power_upper, cumulative_power(end));

% 插值找到对应的频率
lower_freq = interp1(cumulative_power(1:peak_idx), f(1:peak_idx), ...
    power_lower, 'linear', f(1));
upper_freq = interp1(cumulative_power(peak_idx:end), f(peak_idx:end), ...
    power_upper, 'linear', f(end));

B = abs(upper_freq - lower_freq);

results = struct();
results.bandwidth = B;
results.lower_freq = lower_freq;
results.upper_freq = upper_freq;
results.peak_freq = peak_freq;
results.peak_power = peak_power;

end

%===============================================================================
% 绘图函数
%===============================================================================
function plotBandwidthResults(Pxx, f, results, B_estimated)

if isstruct(B_estimated) && isfield(B_estimated, 'threshold')
    % 多个方法的结果
    figure('Name', '带宽估计结果对比（所有方法）', 'Position', [100, 100, 1400, 800]);
    
    methods = {'threshold', 'energy', 'peak_drop', 'rms', 'integral'};
    colors = lines(length(methods));
    
    subplot(2, 2, 1);
    plot(f/1e6, Pxx, 'k-', 'LineWidth', 1);
    hold on;
    for i = 1:length(methods)
        r = results.(methods{i});
        plot([r.lower_freq/1e6, r.lower_freq/1e6], ylim, '--', ...
            'Color', colors(i,:), 'LineWidth', 1.5, ...
            'DisplayName', sprintf('%s下边界', methods{i}));
        plot([r.upper_freq/1e6, r.upper_freq/1e6], '--', ...
            'Color', colors(i,:), 'LineWidth', 1.5, ...
            'DisplayName', sprintf('%s上边界', methods{i}));
    end
    grid on;
    xlabel('频率 (MHz)');
    ylabel('PSD (dB)');
    title('所有方法的带宽估计结果');
    legend('Location', 'best');
    hold off;
    
    subplot(2, 2, 2);
    B_values = [B_estimated.threshold, B_estimated.energy, ...
        B_estimated.peak_drop, B_estimated.rms, B_estimated.integral];
    bar(B_values/1e6);
    set(gca, 'XTickLabel', methods);
    ylabel('带宽 (MHz)');
    title('各方法估计的带宽值');
    grid on;
    
    subplot(2, 2, 3);
    plot(f/1e6, Pxx, 'k-', 'LineWidth', 1);
    hold on;
    r_mean = results.threshold;  % 使用threshold方法作为参考
    plot([r_mean.lower_freq/1e6, r_mean.lower_freq/1e6], ylim, 'r--', ...
        'LineWidth', 2, 'DisplayName', '下边界');
    plot([r_mean.upper_freq/1e6, r_mean.upper_freq/1e6], ylim, 'r--', ...
        'LineWidth', 2, 'DisplayName', '上边界');
    plot([r_mean.peak_freq/1e6, r_mean.peak_freq/1e6], ylim, 'g--', ...
        'LineWidth', 1.5, 'DisplayName', '峰值');
    grid on;
    xlabel('频率 (MHz)');
    ylabel('PSD (dB)');
    title(sprintf('阈值法带宽估计 (%.3f MHz)', B_estimated.threshold/1e6));
    legend('Location', 'best');
    hold off;
    
    subplot(2, 2, 4);
    text(0.1, 0.8, sprintf('平均带宽: %.3f MHz', B_estimated.mean/1e6), ...
        'FontSize', 12, 'Units', 'normalized');
    text(0.1, 0.6, sprintf('标准差: %.3f MHz', B_estimated.std/1e6), ...
        'FontSize', 12, 'Units', 'normalized');
    for i = 1:length(methods)
        text(0.1, 0.8 - i*0.15, ...
            sprintf('%s: %.3f MHz', methods{i}, B_estimated.(methods{i})/1e6), ...
            'FontSize', 10, 'Units', 'normalized');
    end
    axis off;
    title('带宽估计结果汇总');
    
else
    % 单一方法的结果
    figure('Name', sprintf('带宽估计结果（%s方法）', results.method), ...
        'Position', [100, 100, 1200, 600]);
    
    subplot(1, 2, 1);
    plot(f/1e6, Pxx, 'b-', 'LineWidth', 1.5);
    hold on;
    plot([results.lower_freq/1e6, results.lower_freq/1e6], ylim, ...
        'r--', 'LineWidth', 2, 'DisplayName', '下边界');
    plot([results.upper_freq/1e6, results.upper_freq/1e6], ylim, ...
        'r--', 'LineWidth', 2, 'DisplayName', '上边界');
    plot([results.peak_freq/1e6, results.peak_freq/1e6], ylim, ...
        'g--', 'LineWidth', 1.5, 'DisplayName', '峰值');
    
    % 如果方法有阈值，绘制阈值线
    if isfield(results, 'threshold')
        plot([results.lower_freq/1e6, results.upper_freq/1e6], ...
            [results.threshold, results.threshold], 'k:', ...
            'LineWidth', 1, 'DisplayName', sprintf('阈值 (%.1f dB)', results.threshold));
    end
    
    grid on;
    xlabel('频率 (MHz)');
    ylabel('PSD (dB)');
    title(sprintf('%s方法带宽估计\n估计带宽: %.3f MHz', ...
        results.method, B_estimated/1e6));
    legend('Location', 'best');
    hold off;
    
    subplot(1, 2, 2);
    text(0.1, 0.9, sprintf('方法: %s', results.method), ...
        'FontSize', 12, 'Units', 'normalized', 'FontWeight', 'bold');
    text(0.1, 0.75, sprintf('估计带宽: %.3f MHz', B_estimated/1e6), ...
        'FontSize', 11, 'Units', 'normalized');
    text(0.1, 0.65, sprintf('下边界: %.3f MHz', results.lower_freq/1e6), ...
        'FontSize', 10, 'Units', 'normalized');
    text(0.1, 0.55, sprintf('上边界: %.3f MHz', results.upper_freq/1e6), ...
        'FontSize', 10, 'Units', 'normalized');
    text(0.1, 0.45, sprintf('峰值频率: %.3f MHz', results.peak_freq/1e6), ...
        'FontSize', 10, 'Units', 'normalized');
    text(0.1, 0.35, sprintf('峰值功率: %.2f dB', results.peak_power), ...
        'FontSize', 10, 'Units', 'normalized');
    axis off;
    title('估计结果详情');
end

end

