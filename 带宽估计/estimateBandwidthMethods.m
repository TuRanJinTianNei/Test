function [band_lower, band_upper, B_estimated, details] = estimateBandwidthMethods(Pxx, f, method, varargin)
%===============================================================================
% estimateBandwidthMethods.m - 多种带宽估计方法的综合工具箱
% 
% 功能说明:
%   提供多种带宽估计方法，适用于OFDM等信号的带宽估计
%
% 输入参数:
%   Pxx      - 功率谱密度估计值（dB，归一化，行向量或列向量）
%   f        - 对应的频率向量（Hz，行向量或列向量）
%   method   - 估计方法，可选：
%              'threshold'      - 阈值法（默认-3dB）
%              'diff_interp'    - 差分内插法
%              'energy_percent' - 能量百分比法（默认90%）
%              'peak_drop'      - 峰值下降法（默认-3dB）
%              'rms'            - RMS带宽法（二阶矩）
%              'adaptive'       - 自适应阈值法
%              'gradient'       - 基于梯度的边缘检测法
%   varargin - 可选参数（名称-值对）:
%              'threshold_db'      - 阈值法的阈值（dB），默认-3
%              'energy_percent'    - 能量百分比法的百分比，默认90
%              'peak_drop_db'      - 峰值下降法的下降量（dB），默认-3
%              'adaptive_factor'   - 自适应阈值因子，默认0.1
%
% 输出参数:
%   band_lower   - 带宽下边界频率（Hz）
%   band_upper   - 带宽上边界频率（Hz）
%   B_estimated  - 估计的带宽（Hz）
%   details      - 结构体，包含详细的计算信息
%
% 使用示例:
%   % 阈值法
%   [bl, bu, B, d] = estimateBandwidthMethods(Pxx, f, 'threshold', 'threshold_db', -3);
%   
%   % 能量百分比法
%   [bl, bu, B, d] = estimateBandwidthMethods(Pxx, f, 'energy_percent', 'energy_percent', 90);
%
% 创建日期: 2025.12.10
%===============================================================================

% 解析输入参数
p = inputParser;
addParameter(p, 'threshold_db', -3, @isnumeric);
addParameter(p, 'energy_percent', 90, @isnumeric);
addParameter(p, 'peak_drop_db', -3, @isnumeric);
addParameter(p, 'adaptive_factor', 0.1, @isnumeric);
parse(p, varargin{:});

% 确保输入为列向量
Pxx = Pxx(:);
f = f(:);

% 检查输入长度
if length(Pxx) ~= length(f)
    error('Pxx和f的长度必须相同');
end

% 将dB转换为线性功率（用于能量计算）
Pxx_linear = 10.^(Pxx / 10);

% 根据方法选择不同的估计策略
switch lower(method)
    case 'threshold'
        % 方法1: 阈值法
        [band_lower, band_upper, B_estimated, details] = ...
            estimateBandwidthThreshold(Pxx, f, p.Results.threshold_db);
        
    case 'diff_interp'
        % 方法2: 差分内插法（已移除）
        error('差分内插法已移除，请使用其他方法');
        
    case 'energy_percent'
        % 方法3: 能量百分比法
        [band_lower, band_upper, B_estimated, details] = ...
            estimateBandwidthEnergyPercent(Pxx_linear, f, p.Results.energy_percent);
        
    case 'peak_drop'
        % 方法4: 峰值下降法
        [band_lower, band_upper, B_estimated, details] = ...
            estimateBandwidthPeakDrop(Pxx, f, p.Results.peak_drop_db);
        
    case 'rms'
        % 方法5: RMS带宽法（二阶矩）
        [band_lower, band_upper, B_estimated, details] = ...
            estimateBandwidthRMS(Pxx_linear, f);
        
    case 'adaptive'
        % 方法6: 自适应阈值法
        [band_lower, band_upper, B_estimated, details] = ...
            estimateBandwidthAdaptive(Pxx, f, p.Results.adaptive_factor);
        
    case 'gradient'
        % 方法7: 基于梯度的边缘检测法
        [band_lower, band_upper, B_estimated, details] = ...
            estimateBandwidthGradient(Pxx, f);
        
    otherwise
        error('未知的方法: %s', method);
end

details.method = method;

end

%===============================================================================
% 辅助函数：各种带宽估计方法的具体实现
%===============================================================================

% 方法1: 阈值法
function [band_lower, band_upper, B_estimated, details] = estimateBandwidthThreshold(Pxx, f, threshold_db)
    [~, peak_idx] = max(Pxx);
    peak_freq = f(peak_idx);
    
    % 在峰值左侧找阈值点
    left_side = Pxx(1:peak_idx);
    left_indices = find(left_side <= threshold_db);
    if ~isempty(left_indices)
        left_idx = left_indices(end);
        band_lower = f(left_idx);
    else
        band_lower = f(1);
    end
    
    % 在峰值右侧找阈值点
    right_side = Pxx(peak_idx:end);
    right_indices = find(right_side <= threshold_db);
    if ~isempty(right_indices)
        right_idx = peak_idx + right_indices(1) - 1;
        band_upper = f(right_idx);
    else
        band_upper = f(end);
    end
    
    B_estimated = abs(band_upper - band_lower);
    details.threshold_db = threshold_db;
    details.peak_idx = peak_idx;
    details.peak_freq = peak_freq;
end

% 方法3: 能量百分比法
function [band_lower, band_upper, B_estimated, details] = estimateBandwidthEnergyPercent(Pxx_linear, f, energy_percent)
    % 计算总能量
    total_energy = trapz(f, Pxx_linear);
    target_energy = total_energy * energy_percent / 100;
    
    % 找到峰值位置
    [~, peak_idx] = max(Pxx_linear);
    
    % 从峰值向两侧扩展，直到累积能量达到目标百分比
    % 左侧
    cum_energy_left = 0;
    left_idx = peak_idx;
    for i = peak_idx:-1:1
        if i < length(f)
            segment_energy = trapz(f(i:i+1), Pxx_linear(i:i+1));
            cum_energy_left = cum_energy_left + segment_energy;
            if cum_energy_left >= target_energy / 2
                left_idx = i;
                break;
            end
        end
    end
    
    % 右侧
    cum_energy_right = 0;
    right_idx = peak_idx;
    for i = peak_idx:length(f)-1
        segment_energy = trapz(f(i:i+1), Pxx_linear(i:i+1));
        cum_energy_right = cum_energy_right + segment_energy;
        if cum_energy_right >= target_energy / 2
            right_idx = i + 1;
            break;
        end
    end
    
    band_lower = f(left_idx);
    band_upper = f(right_idx);
    B_estimated = abs(band_upper - band_lower);
    
    details.energy_percent = energy_percent;
    details.total_energy = total_energy;
    details.target_energy = target_energy;
    details.left_idx = left_idx;
    details.right_idx = right_idx;
end

% 方法4: 峰值下降法
function [band_lower, band_upper, B_estimated, details] = estimateBandwidthPeakDrop(Pxx, f, peak_drop_db)
    [peak_value, peak_idx] = max(Pxx);
    peak_freq = f(peak_idx);
    threshold = peak_value + peak_drop_db;  % 峰值下降peak_drop_db dB
    
    % 在峰值左侧找阈值点
    left_side = Pxx(1:peak_idx);
    left_indices = find(left_side <= threshold);
    if ~isempty(left_indices)
        left_idx = left_indices(end);
        band_lower = f(left_idx);
    else
        band_lower = f(1);
    end
    
    % 在峰值右侧找阈值点
    right_side = Pxx(peak_idx:end);
    right_indices = find(right_side <= threshold);
    if ~isempty(right_indices)
        right_idx = peak_idx + right_indices(1) - 1;
        band_upper = f(right_idx);
    else
        band_upper = f(end);
    end
    
    B_estimated = abs(band_upper - band_lower);
    details.peak_value = peak_value;
    details.peak_idx = peak_idx;
    details.peak_freq = peak_freq;
    details.threshold = threshold;
    details.peak_drop_db = peak_drop_db;
end

% 方法5: RMS带宽法（二阶矩）
function [band_lower, band_upper, B_estimated, details] = estimateBandwidthRMS(Pxx_linear, f)
    % 归一化PSD
    Pxx_norm = Pxx_linear / trapz(f, Pxx_linear);
    
    % 计算一阶矩（中心频率）
    f_mean = trapz(f, f .* Pxx_norm);
    
    % 计算二阶矩（RMS带宽）
    f_var = trapz(f, (f - f_mean).^2 .* Pxx_norm);
    rms_bandwidth = 2 * sqrt(f_var);  % 2倍标准差
    
    % 带宽边界
    band_lower = f_mean - rms_bandwidth / 2;
    band_upper = f_mean + rms_bandwidth / 2;
    
    % 限制在频率范围内
    band_lower = max(band_lower, f(1));
    band_upper = min(band_upper, f(end));
    
    B_estimated = abs(band_upper - band_lower);
    
    details.f_mean = f_mean;
    details.f_var = f_var;
    details.rms_bandwidth = rms_bandwidth;
end

% 方法6: 自适应阈值法
function [band_lower, band_upper, B_estimated, details] = estimateBandwidthAdaptive(Pxx, f, adaptive_factor)
    % 计算PSD的动态范围
    peak_value = max(Pxx);
    noise_floor = prctile(Pxx, 10);  % 使用10%分位数作为噪声底
    dynamic_range = peak_value - noise_floor;
    
    % 自适应阈值 = 峰值 - 动态范围 * 因子
    threshold = peak_value - dynamic_range * adaptive_factor;
    
    [~, peak_idx] = max(Pxx);
    
    % 在峰值左侧找阈值点
    left_side = Pxx(1:peak_idx);
    left_indices = find(left_side <= threshold);
    if ~isempty(left_indices)
        left_idx = left_indices(end);
        band_lower = f(left_idx);
    else
        band_lower = f(1);
    end
    
    % 在峰值右侧找阈值点
    right_side = Pxx(peak_idx:end);
    right_indices = find(right_side <= threshold);
    if ~isempty(right_indices)
        right_idx = peak_idx + right_indices(1) - 1;
        band_upper = f(right_idx);
    else
        band_upper = f(end);
    end
    
    B_estimated = abs(band_upper - band_lower);
    
    details.peak_value = peak_value;
    details.noise_floor = noise_floor;
    details.dynamic_range = dynamic_range;
    details.adaptive_factor = adaptive_factor;
    details.threshold = threshold;
end

% 方法7: 基于梯度的边缘检测法
function [band_lower, band_upper, B_estimated, details] = estimateBandwidthGradient(Pxx, f)
    % 计算梯度（一阶差分）
    gradient = diff(Pxx);
    gradient_abs = abs(gradient);
    
    % 找到峰值位置
    [~, peak_idx] = max(Pxx);
    
    % 在峰值左侧找最大梯度（上升沿）
    left_gradient = gradient_abs(1:peak_idx-1);
    if ~isempty(left_gradient)
        [~, left_max_grad_idx] = max(left_gradient);
        band_lower = f(left_max_grad_idx);
    else
        band_lower = f(1);
    end
    
    % 在峰值右侧找最大梯度（下降沿）
    if peak_idx < length(gradient_abs)
        right_gradient = gradient_abs(peak_idx:end);
        [~, right_max_grad_idx] = max(right_gradient);
        band_upper = f(peak_idx + right_max_grad_idx - 1);
    else
        band_upper = f(end);
    end
    
    B_estimated = abs(band_upper - band_lower);
    
    details.gradient = gradient;
    details.peak_idx = peak_idx;
    details.left_max_grad_idx = left_max_grad_idx;
    details.right_max_grad_idx = peak_idx + right_max_grad_idx - 1;
end
