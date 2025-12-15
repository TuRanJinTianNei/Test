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
% 创建日期: 2025.12.10
% 基于: method.m 中的 PSD_OFDM_rayleigh 函数
% 修改日期: 2025.12.10 - 移除带宽计算功能，仅保留功率谱密度估计
%===============================================================================

% 解析可选参数
p = inputParser;
addParameter(p, 'plot', false, @islogical);
addParameter(p, 'welch_window', 100, @isnumeric);
addParameter(p, 'welch_overlap', 55, @isnumeric);
addParameter(p, 'welch_nfft', 8192, @isnumeric);
addParameter(p, 'ar_criterion', 'AIC', @ischar);
parse(p, varargin{:});

plot_flag = p.Results.plot;
welch_window = p.Results.welch_window;
welch_overlap = p.Results.welch_overlap;
welch_nfft = p.Results.welch_nfft;
ar_criterion = p.Results.ar_criterion;

%**************************************************************************
% Welch算法功率谱估计
%**************************************************************************
% 注意：pwelch对实信号只返回正频率部分（0到fs/2）
[Pxx2, f1] = pwelch(sig_processed, hanning(welch_window), welch_overlap, welch_nfft, fs);

% 归一化处理
Pxx22 = Pxx2;
Pxx22 = Pxx22 / min(Pxx22);  % 归一化到最小值
Pxx22 = 10*log10(Pxx22);     % 转换为dB单位
Pxx22 = Pxx22 - max(Pxx22);  % 归一化到最大值（峰值在0dB）

% 输出Welch功率谱（正频率部分：0到fs/2）
Pxx_welch = Pxx22;
f_welch = f1;

%**************************************************************************
% AR模型功率谱估计（使用Burg算法）
%**************************************************************************
% 注意：Burg函数内部已处理，只返回正频率部分（0到fs/2）
[Pxx1, f, p] = Burg(sig_processed, fs, ar_criterion);

% 输出AR功率谱（正频率部分：0到fs/2）
Pxx_ar = Pxx1;
f_ar = f;

%**************************************************************************
% 绘制PSD图形（可选）
%**************************************************************************
if plot_flag
    % AR模型功率谱图
    figure('Name', 'AR模型功率谱密度估计')
    plot(f, Pxx1);
    grid on;
    xlabel('频率 f (Hz)');
    ylabel('PSD (dB)');
    title(sprintf('AR模型方法的功率谱密度估计 (SNR = %.1f dB, 阶数 = %d)', snr, p));
    
    % Welch算法功率谱图
    figure('Name', 'Welch算法功率谱密度估计')
    plot(f1, Pxx22);
    grid on;
    xlabel('频率 f (Hz)');
    ylabel('PSD (dB)');
    title(sprintf('Welch算法估计的功率谱密度 (SNR = %.1f dB)', snr));
    
    % OFDM信号频谱图
    N_sig = length(sig_processed);
    N_fft = N_sig;  % 使用信号长度作为FFT点数
    OfdmSymComput = 20 * log10(abs(fft(sig_processed, N_fft)));
    OfdmSymPSDy = fftshift(OfdmSymComput) - max(OfdmSymComput);
    
    % 生成频率向量（只取正频率部分）
    N_half = floor(N_fft/2);
    Frc = (0:N_half-1) * fs / N_fft;  % 0 到 fs/2 的频率范围
    
    % 提取正频率部分的频谱（fftshift后，正频率在中间到末尾）
    if mod(N_fft, 2) == 0
        % 偶数长度：正频率从 N_fft/2+1 到 N_fft
        OfdmSymPSDy_positive = OfdmSymPSDy(N_fft/2+1:end);
    else
        % 奇数长度：正频率从 (N_fft+1)/2+1 到 N_fft
        OfdmSymPSDy_positive = OfdmSymPSDy((N_fft+1)/2+1:end);
    end
    
    % 确保长度匹配（取较小的长度）
    min_len = min(length(Frc), length(OfdmSymPSDy_positive));
    Frc = Frc(1:min_len);
    OfdmSymPSDy_positive = OfdmSymPSDy_positive(1:min_len);
    
    figure('Name', 'OFDM信号频谱')
    plot(Frc, OfdmSymPSDy_positive);
    xlabel('频率 f (Hz)');
    ylabel('PSD (dB)');
    title(sprintf('OFDM信号频谱 (SNR = %.1f dB)', snr));
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
%
% 根据输入参数类型判断
if isnumeric(varargin{1}) && isscalar(varargin{1})
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
    p = ceil(N/3); % 阶数一般不超出信号长度的1/3
    
    % 计算1到p阶的误差（用于选择最优阶数）
    % 注意：这里的a未使用，但E用于计算目标函数goalF
    [~, E] = computeARpara(x, p);
    
    % 计算目标函数的最小值
    kc = 1:p + 1;
    switch criterion
        case 'FPE'
            goalF = E.*(N + (kc + 1))./(N - (kc + 1));
        case 'AIC'
            goalF = N.*log(E) + 2.*kc;
    end
    [~, p] = min(goalF); % p是目标函数最小值的位值，也即准则准则确定的阶数
    
    % 使用p值重新计算AR模型参数
    [a, E] = computeARpara(x, p);
end
[h, f] = freqz(1, a, 20e5, Fs);
psdviaBurg = E(end)*abs(h).^2./Fs;
psdviaBurg = psdviaBurg/abs(max(psdviaBurg));
psdviaBurg = (10*log10(abs(psdviaBurg)));

% 只保留正频率部分（0到Fs/2），与Welch算法保持一致
% freqz返回的频率范围是0到Fs，对于实信号，负频率部分是正频率的镜像
f_nyquist = Fs / 2;  % 奈奎斯特频率
positive_freq_idx = f <= f_nyquist;
f = f(positive_freq_idx);
psdviaBurg = psdviaBurg(positive_freq_idx);
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
