function [Pxx_welch, f_welch, Pxx_ar, f_ar] = estimateBandwidthPSD(sig_processed, fs, snr, varargin)
%===============================================================================
% estimateBandwidthPSD.m - 使用Welch算法和AR模型法估计信号功率谱密度
% 
% 功能说明:
%   对已处理的信号（上变频、Rayleigh信道、加噪声）进行功率谱密度估计
%   使用两种方法：
%   1. Welch算法（周期图平均法）
%   2. AR模型法（自回归模型，使用Burg算法）
%
% 输入参数:
%   sig_processed  - 已处理的信号（实数信号，行向量）
%                    通常已经过：上变频、Rayleigh衰落信道、AWGN噪声
%   fs             - 采样频率（Hz）
%   snr            - 信噪比（dB），用于绘图标题显示
%   varargin       - 可选参数（名称-值对）:
%                    'plot'      - 是否绘制PSD图，默认false
%                    'welch_window' - Welch算法窗长度，默认100
%                    'welch_overlap' - Welch算法重叠样本数，默认55
%                    'welch_nfft' - Welch算法FFT点数，默认8192
%                    'ar_criterion' - AR模型阶数选择准则，默认'AIC'
%                    'ar_max_order_ratio' - AR模型最大阶数比例，默认0.5
%                                           (阶数上限 = ceil(N * ar_max_order_ratio))
%                                           增加此值可获得更平滑的PSD，但可能过拟合
%
% 输出参数:
%   Pxx_welch      - Welch算法功率谱密度估计值（dB，归一化）
%   f_welch        - Welch算法对应的频率向量（Hz）
%   Pxx_ar         - AR模型功率谱密度估计值（dB，归一化）
%   f_ar           - AR模型对应的频率向量（Hz）
%
% 使用示例:
%   % 基本用法
%   [Pxx_w, f_w, Pxx_a, f_a] = estimateBandwidthPSD(sig_processed, fs, snr);
%
%   % 绘制PSD图
%   [Pxx_w, f_w, Pxx_a, f_a] = estimateBandwidthPSD(...
%       sig_processed, fs, snr, 'plot', true);
%
%   % 增加AR模型阶数上限以获得更平滑的PSD（默认0.5，即N/2）
%   [Pxx_w, f_w, Pxx_a, f_a] = estimateBandwidthPSD(...
%       sig_processed, fs, snr, 'ar_max_order_ratio', 0.6);
%
% 创建日期: 2025.12.10
% 基于: method.m 中的 PSD_OFDM_rayleigh 函数
% 修改日期: 2025.12.10 - 移除带宽计算功能，仅保留功率谱密度估计
%           2025.12.16 - 增加AR模型阶数上限参数，支持更平滑的PSD估计
%===============================================================================

% 解析可选参数
p = inputParser;
addParameter(p, 'plot', false, @islogical);
addParameter(p, 'welch_window', 100, @isnumeric);
addParameter(p, 'welch_overlap', 55, @isnumeric);
addParameter(p, 'welch_nfft', 8192, @isnumeric);
addParameter(p, 'ar_criterion', 'AIC', @ischar);
addParameter(p, 'ar_max_order_ratio', 0.5, @isnumeric);  % 默认0.5，即N/2
parse(p, varargin{:});

plot_flag = p.Results.plot;
welch_window = p.Results.welch_window;
welch_overlap = p.Results.welch_overlap;
welch_nfft = p.Results.welch_nfft;
ar_criterion = p.Results.ar_criterion;
ar_max_order_ratio = p.Results.ar_max_order_ratio;

%**************************************************************************
% Welch算法功率谱估计
%**************************************************************************
% 注意：对于实信号，频谱是共轭对称的，即 X(-f) = X*(f)
%      pwelch函数对实信号自动返回单边谱（0到fs/2），并且已经自动将
%      正频率部分的幅度乘以2（除了DC和Nyquist频率），所以输出已经是单边谱
%      单边谱幅度 = 2 × 双边谱幅度（对于0 < f < fs/2）
%      DC（f=0）和Nyquist频率（f=fs/2）的幅度保持不变
[Pxx2, f1] = pwelch(sig_processed, hanning(welch_window), welch_overlap, welch_nfft, fs);

% 验证：确保频率范围正确（0到fs/2）
if max(f1) > fs/2 + eps
    warning('Welch算法返回的频率范围超出正频率范围，可能存在错误');
end

% 确保pwelch返回的是单边谱（0到fs/2）
% 如果频率范围不正确，需要手动提取正频率部分
if min(f1) < -eps || max(f1) > fs/2 + eps
    % 提取正频率部分（0到fs/2）
    positive_freq_idx = (f1 >= 0) & (f1 <= fs/2 + eps);
    f1 = f1(positive_freq_idx);
    Pxx2 = Pxx2(positive_freq_idx);
end

% 显式验证并确保单边谱转换正确
% 
% 【单边谱转换规则说明】
% 对于实信号，双边功率谱密度（PSD）是共轭对称的，即 P(-f) = P(f)
% 单边谱转换规则：
%   - DC分量（f=0）：幅度不变（P_single(0) = P_double(0)）
%   - 正频率（0 < f < fs/2）：幅度乘以2（P_single(f) = 2 × P_double(f)）
%   - Nyquist频率（f=fs/2）：幅度不变（P_single(fs/2) = P_double(fs/2)）
% 
% 【pwelch函数行为】
% MATLAB的pwelch函数对实信号自动返回单边谱（0到fs/2），并且已经自动将
% 正频率部分的功率谱密度乘以2（除了DC和Nyquist频率），所以输出已经是单边谱。
% 因此，我们不需要再次进行转换，只需要验证频率范围确实是0到fs/2。
% 
% 【验证步骤】
% 1. 验证频率范围：确保f1在[0, fs/2]范围内
% 2. 如果频率范围不正确，提取正频率部分
% 3. 确保输出是单边谱格式

% 归一化处理（在转换为dB之前）
Pxx22 = Pxx2;
Pxx22 = Pxx22 / min(Pxx22);  % 归一化到最小值
Pxx22 = 10*log10(Pxx22);     % 转换为dB单位
Pxx22 = Pxx22 - max(Pxx22);  % 归一化到最大值（峰值在0dB）

% 输出Welch功率谱（单边谱：0到fs/2）
Pxx_welch = Pxx22;
f_welch = f1;

%**************************************************************************
% AR模型功率谱估计（使用Burg算法）
%**************************************************************************
% 注意：对于实信号，频谱是共轭对称的
%      Burg函数内部会将双边谱转换为单边谱（0到fs/2），与Welch算法保持一致
%      单边谱幅度 = 2 × 双边谱幅度（对于0 < f < fs/2）
%      DC（f=0）和Nyquist频率（f=fs/2）的幅度保持不变
[Pxx1, f, p] = Burg(sig_processed, fs, ar_criterion, ar_max_order_ratio);

% 输出AR功率谱（正频率部分：0到fs/2）
Pxx_ar = Pxx1;
f_ar = f;

%**************************************************************************
% 绘制PSD图形（可选）
%**************************************************************************
if plot_flag
    % AR模型功率谱图（单边谱）
    figure('Name', 'AR模型功率谱密度估计（单边谱）')
    plot(f, Pxx1);
    grid on;
    xlabel('频率 f (Hz)');
    ylabel('PSD (dB)');
    title(sprintf('AR模型方法的功率谱密度估计（单边谱）(SNR = %.1f dB, 阶数 = %d)', snr, p));
    
    % Welch算法功率谱图（单边谱）
    figure('Name', 'Welch算法功率谱密度估计（单边谱）')
    plot(f1, Pxx22);
    grid on;
    xlabel('频率 f (Hz)');
    ylabel('PSD (dB)');
    title(sprintf('Welch算法估计的功率谱密度（单边谱）(SNR = %.1f dB)', snr));
    
    % OFDM信号频谱图（单边谱）
    N_sig = length(sig_processed);
    N_fft = N_sig;  % 使用信号长度作为FFT点数
    
    % FFT计算（双边谱）
    X_fft = fft(sig_processed, N_fft);
    X_mag_squared = abs(X_fft).^2;  % 功率谱（双边）
    
    % 提取正频率部分并转换为单边谱
    if mod(N_fft, 2) == 0
        % 偶数长度：正频率从 1 到 N_fft/2+1
        % 索引：1 (DC), 2:N_fft/2 (正频率), N_fft/2+1 (Nyquist)
        X_single_sided = X_mag_squared(1:N_fft/2+1);
        % 单边谱转换：正频率部分（除DC和Nyquist）乘以2
        X_single_sided(2:N_fft/2) = X_single_sided(2:N_fft/2) * 2;
        % 频率向量（0到fs/2）
        Frc = (0:N_fft/2) * fs / N_fft;
    else
        % 奇数长度：正频率从 1 到 (N_fft+1)/2
        % 索引：1 (DC), 2:(N_fft+1)/2 (正频率，无Nyquist频率)
        X_single_sided = X_mag_squared(1:(N_fft+1)/2);
        % 单边谱转换：正频率部分（除DC）乘以2
        X_single_sided(2:end) = X_single_sided(2:end) * 2;
        % 频率向量（0到略小于fs/2）
        Frc = (0:(N_fft-1)/2) * fs / N_fft;
    end
    
    % 转换为dB并归一化
    OfdmSymPSDy = 10 * log10(X_single_sided);
    OfdmSymPSDy = OfdmSymPSDy - max(OfdmSymPSDy);  % 归一化到峰值0dB
    
    figure('Name', 'OFDM信号频谱（单边谱）')
    plot(Frc/1e6, OfdmSymPSDy);
    xlabel('频率 f (MHz)');
    ylabel('PSD (dB)');
    title(sprintf('OFDM信号频谱（单边谱）(SNR = %.1f dB)', snr));
    grid on;
end

end

%**************************************************************************
% 以下为内部辅助函数
%**************************************************************************

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
%            如果为字符串，可额外提供最大阶数比例参数（默认0.33，即N/3）
%
% 根据输入参数类型判断
if isnumeric(varargin{1}) && isscalar(varargin{1})
    p = varargin{1};
elseif ischar(varargin{1})
    criterion = varargin{1};
    % 检查是否有最大阶数比例参数
    if length(varargin) >= 2 && isnumeric(varargin{2})
        max_order_ratio = varargin{2};
    else
        max_order_ratio = 0.33;  % 默认N/3（保持向后兼容）
    end
else
    error('第2个参数必须为数值型或字符串');
end
x = x(:);
N = length(x);
% 模型参数估计
if exist('p', 'var') % p变量是否存在，如果存在则不需要估计，直接使用p值
    [a, E] = computeARpara(x, p);
else % p不存在，需要估计，根据准则criterion
    p_max = ceil(N * max_order_ratio);  % 阶数上限，可配置
    
    % 计算到p_max阶的AR模型参数
    % 注意：computeARpara返回的E数组包含0到p_max阶的所有误差值
    [~, E] = computeARpara(x, p_max);
    
    % 计算目标函数的最小值（从0阶到p_max阶）
    kc = 0:p_max;  % 阶数从0到p_max
    switch criterion
        case 'FPE'
            goalF = E.*(N + (kc + 1))./(N - (kc + 1));
        case 'AIC'
            goalF = N.*log(E) + 2.*kc;
    end
    [~, p_idx] = min(goalF);  % 找到最小值的索引
    p = p_idx - 1;  % 转换为实际阶数（因为索引从1开始，阶数从0开始）
    
    % 使用最优阶数p重新计算AR模型参数（以获得最终的a和E）
    [a, E] = computeARpara(x, p);
end
[h, f] = freqz(1, a, 20e5, Fs);
% freqz返回的是双边谱（0到Fs），需要转换为单边谱
% 单边谱幅度规则：
%   - DC分量（f=0）：幅度不变
%   - 正频率（0 < f < fs/2）：幅度乘以2
%   - Nyquist频率（f=fs/2）：幅度不变（如果存在）
psdviaBurg = E(end)*abs(h).^2./Fs;

% 只保留正频率部分（0到Fs/2）
f_nyquist = Fs / 2;  % 奈奎斯特频率
positive_freq_idx = f <= f_nyquist + eps;  % 包含Nyquist频率
f = f(positive_freq_idx);
psdviaBurg = psdviaBurg(positive_freq_idx);

% 转换为单边谱：正频率部分（除DC和Nyquist）幅度乘以2
% 单边谱规则：
%   - DC分量（f=0）：幅度不变
%   - 正频率（0 < f < fs/2）：幅度乘以2
%   - Nyquist频率（f=fs/2）：幅度不变

% 找到DC和Nyquist频率的索引
dc_idx = find(abs(f) < eps, 1);  % DC频率（f=0）
nyquist_idx = find(abs(f - f_nyquist) < eps, 1);  % Nyquist频率（f=fs/2）

% 对正频率部分（0 < f < fs/2）乘以2
if ~isempty(dc_idx) && ~isempty(nyquist_idx)
    % 有DC和Nyquist频率：对中间部分乘以2
    if nyquist_idx > dc_idx + 1
        psdviaBurg(dc_idx+1:nyquist_idx-1) = psdviaBurg(dc_idx+1:nyquist_idx-1) * 2;
    end
elseif ~isempty(dc_idx)
    % 只有DC频率，没有Nyquist频率（奇数长度或频率分辨率问题）
    psdviaBurg(dc_idx+1:end) = psdviaBurg(dc_idx+1:end) * 2;
elseif ~isempty(nyquist_idx)
    % 只有Nyquist频率，没有DC频率（不太可能，但处理一下）
    psdviaBurg(1:nyquist_idx-1) = psdviaBurg(1:nyquist_idx-1) * 2;
else
    % 都没有，全部乘以2（不太可能，但处理一下）
    psdviaBurg = psdviaBurg * 2;
end

% 归一化和转换为dB
psdviaBurg = psdviaBurg/abs(max(psdviaBurg));
psdviaBurg = (10*log10(abs(psdviaBurg)));
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
