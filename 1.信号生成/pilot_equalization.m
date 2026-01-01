function [Rx_carriers, Rx_carriers_before_eq, H_est_full] = pilot_equalization(...
    Y1, pilot_carriers, data_carriers, pilot_conjugate_carriers, ...
    data_conjugate_carriers, pilot_symbol, interpolation_method, ...
    IFFT_bin_length, symbols_per_carrier, use_pilot_equalization)
%PILOT_EQUALIZATION 导频辅助信道估计与频域均衡
%   实现基于导频的信道估计和零迫均衡（Zero-Forcing Equalization）
%
%   输入参数：
%       Y1                      - FFT后的接收信号（所有子载波），矩阵 [symbols_per_carrier × IFFT_bin_length]
%       pilot_carriers          - 导频子载波位置（绝对bin索引）
%       data_carriers           - 数据子载波位置（绝对bin索引）
%       pilot_conjugate_carriers - 导频共轭位置（绝对bin索引）
%       data_conjugate_carriers  - 数据共轭位置（绝对bin索引）
%       pilot_symbol            - 导频符号（已知的复数值）
%       interpolation_method     - 插值方法：'linear'（线性插值）或'spline'（样条插值）
%       IFFT_bin_length         - IFFT/FFT点数
%       symbols_per_carrier     - OFDM符号数量
%       use_pilot_equalization  - 是否使用导频均衡（true=使用，false=不使用）
%
%   输出参数：
%       Rx_carriers            - 均衡后的数据子载波信号，矩阵 [symbols_per_carrier × data_carrier_count]
%       Rx_carriers_before_eq  - 均衡前的数据子载波信号（用于对比），矩阵 [symbols_per_carrier × data_carrier_count]
%       H_est_full             - 完整的信道估计（所有子载波），矩阵 [symbols_per_carrier × IFFT_bin_length]
%                                如果未使用均衡，返回空矩阵
%
%   功能说明：
%       1. 信道估计：基于导频位置的信道估计 H_est = Y_pilot / X_pilot
%       2. 插值估计：对非导频位置的信道响应进行插值（线性或样条插值）
%       3. 频域均衡：零迫均衡（ZF）：Y_eq = Y / H_est
%       4. 除零保护：当信道估计过小时使用最小阈值（1e-10）避免数值问题
%
%   使用示例：
%       [Rx_eq, Rx_before, H_est] = pilot_equalization(Y1, pilot_carriers, ...
%           data_carriers, pilot_conjugate_carriers, data_conjugate_carriers, ...
%           pilot_symbol, 'linear', 1024, 100, true);

% 提取数据子载波位置的接收信号（均衡前）
Y_data_before_eq = Y1(:, data_carriers);

% 计算导频数量
pilot_count = length(pilot_carriers);

if use_pilot_equalization && pilot_count > 0
    %--------------------------------------------------------------------------
    % 步骤1: 提取导频位置的接收信号
    %--------------------------------------------------------------------------
    Y_pilot = Y1(:, pilot_carriers);  % 导频位置的接收信号
    
    %--------------------------------------------------------------------------
    % 步骤2: 信道估计：H_est = Y_pilot / X_pilot（在导频位置）
    %--------------------------------------------------------------------------
    % X_pilot是已知的导频符号
    H_est_pilot = Y_pilot ./ repmat(pilot_symbol, symbols_per_carrier, pilot_count);
    
    %--------------------------------------------------------------------------
    % 步骤3: 对非导频位置的信道响应进行插值
    %--------------------------------------------------------------------------
    % 为每个符号分别进行插值
    H_est_full = zeros(symbols_per_carrier, IFFT_bin_length);
    
    for sym_idx = 1:symbols_per_carrier
        % 提取当前符号在导频位置的信道估计
        H_pilot_sym = H_est_pilot(sym_idx, :);
        
        % 创建插值点：导频位置
        pilot_bins = pilot_carriers;
        % 创建查询点：所有数据子载波位置
        data_bins = data_carriers;
        
        % 使用插值方法估计数据子载波位置的信道响应
        if strcmp(interpolation_method, 'linear')
            H_data_sym = interp1(pilot_bins, H_pilot_sym, data_bins, 'linear', 'extrap');
        else  % spline
            H_data_sym = interp1(pilot_bins, H_pilot_sym, data_bins, 'spline', 'extrap');
        end
        
        % 存储估计的信道响应
        H_est_full(sym_idx, pilot_carriers) = H_pilot_sym;
        H_est_full(sym_idx, data_carriers) = H_data_sym;
        
        % 共轭对称：负频率位置的信道响应
        H_est_full(sym_idx, pilot_conjugate_carriers) = conj(H_pilot_sym);
        H_est_full(sym_idx, data_conjugate_carriers) = conj(H_data_sym);
    end
    
    %--------------------------------------------------------------------------
    % 步骤4: 提取数据子载波位置的接收信号和信道估计
    %--------------------------------------------------------------------------
    Y_data = Y1(:, data_carriers);
    H_est_data = H_est_full(:, data_carriers);
    
    %--------------------------------------------------------------------------
    % 步骤5: 频域均衡：Y_eq = Y / H_est（零迫均衡，Zero-Forcing）
    %--------------------------------------------------------------------------
    % 避免除零：当信道估计很小时，使用最小阈值
    H_est_data_safe = H_est_data;
    H_est_data_safe(abs(H_est_data_safe) < 1e-10) = 1e-10;
    Y_eq = Y_data ./ H_est_data_safe;
    
    % 均衡后的数据子载波
    Rx_carriers = Y_eq;
    
    % 保存均衡前的数据用于对比显示
    Rx_carriers_before_eq = Y_data_before_eq;
else
    % 不使用均衡，直接提取数据子载波
    Rx_carriers = Y1(:, data_carriers);
    Rx_carriers_before_eq = Y_data_before_eq;  % 均衡前后相同（无均衡）
    H_est_full = [];  % 未使用均衡，返回空矩阵
end

end

