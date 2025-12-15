function [B_welch, B_ar, Pxx_welch, f_welch, Pxx_ar, f_ar] = estimateBandwidthPSD(sig_processed, fs, snr, varargin)
%===============================================================================
% estimateBandwidthPSD.m - 使用Welch算法和AR模型法估计信号带宽
% 
% 功能说明:
%   对已处理的信号（上变频、Rayleigh信道、加噪声）进行功率谱估计和带宽计算
%   使用两种方法：
%   1. Welch算法（周期图平均法）
%   2. AR模型法（自回归模型，使用Burg算法）
%
% 输入参数:
%   sig_processed  - 已处理的信号（实数信号，行向量）
%                    通常已经过：上变频、Rayleigh衰落信道、AWGN噪声
%   fs             - 采样频率（Hz）
%   snr            - 信噪比（dB），用于自适应阈值选择
%   varargin       - 可选参数（名称-值对）:
%                    'plot'      - 是否绘制PSD图，默认false
%                    'welch_window' - Welch算法窗长度，默认100
%                    'welch_overlap' - Welch算法重叠样本数，默认55
%                    'welch_nfft' - Welch算法FFT点数，默认8192
%                    'ar_criterion' - AR模型阶数选择准则，默认'AIC'
%
% 输出参数:
%   B_welch        - Welch算法估计的带宽值（Hz）
%   B_ar           - AR模型法估计的带宽值（Hz）
%   Pxx_welch      - Welch算法功率谱密度估计值（dB，归一化）
%   f_welch        - Welch算法对应的频率向量（Hz）
%   Pxx_ar         - AR模型功率谱密度估计值（dB，归一化）
%   f_ar           - AR模型对应的频率向量（Hz）
%
% 带宽计算方法:
%   - Welch算法：使用固定-3dB阈值
%   - AR模型法：根据SNR自适应选择阈值（-6dB/-5dB/-3dB）
%
% 使用示例:
%   % 基本用法
%   [B_welch, B_ar] = estimateBandwidthPSD(sig_processed, fs, snr);
%
%   % 绘制PSD图
%   [B_welch, B_ar, Pxx_w, f_w, Pxx_a, f_a] = estimateBandwidthPSD(...
%       sig_processed, fs, snr, 'plot', true);
%
% 创建日期: 2025.12.10
% 基于: method.m 中的 PSD_OFDM_rayleigh 函数
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
[Pxx2, f1] = pwelch(sig_processed, hanning(welch_window), welch_overlap, welch_nfft, fs);

% 归一化处理
Pxx22 = Pxx2;
Pxx22 = Pxx22 / min(Pxx22);  % 归一化到最小值
Pxx22 = 10*log10(Pxx22);     % 转换为dB单位
Pxx22 = Pxx22 - max(Pxx22);  % 归一化到最大值（峰值在0dB）

% 输出Welch功率谱
Pxx_welch = Pxx22;
f_welch = f1;

%**************************************************************************
% AR模型功率谱估计（使用Burg算法）
%**************************************************************************
[Pxx1, f, p] = Burg(sig_processed, fs, ar_criterion);

% 输出AR功率谱
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
    Frc = 0:fs/(length(sig_processed)):fs/2-1;
    OfdmSymComput = 20 * log10(abs(fft(sig_processed)));
    OfdmSymPSDy = fftshift(OfdmSymComput) - max(OfdmSymComput);
    
    figure('Name', 'OFDM信号频谱')
    plot(Frc, OfdmSymPSDy(1, 1:end/2));
    xlabel('频率 f (Hz)');
    ylabel('PSD (dB)');
    title(sprintf('OFDM信号频谱 (SNR = %.1f dB)', snr));
end

%**************************************************************************
% 计算信号的带宽
%**************************************************************************

% Welch算法带宽计算：使用固定-3dB阈值
L1 = ceil(length(Pxx22)/2);
P1 = Pxx22(1:L1, 1);      % 前半部分（低频）
P2 = Pxx22(L1:end, 1);    % 后半部分（高频）

[as1, ~] = Proximate(-3, P1);  % 取最接近-3dB处信号的f值
band1 = f1(as1);
[as2, ~] = Proximate(-3, P2);
band2 = f1(as2 + L1 - 1);
B_welch = abs(band1 - band2);

% AR模型带宽计算：根据SNR自适应选择阈值
L2 = ceil(length(Pxx1)/2);
P3 = Pxx1(1:L2, 1);       % 前半部分（低频）
P4 = Pxx1(L2:end, 1);     % 后半部分（高频）

% 根据SNR选择阈值
if snr > 4
    threshold = -6;  % 高SNR：使用-6dB阈值
elseif snr > 0 && snr <= 4
    threshold = -5;  % 中等SNR：使用-5dB阈值
else
    threshold = -3;  % 低SNR：使用-3dB阈值
end

[as3, ~] = Proximate(threshold, P3);  % 取最接近阈值处信号的f值
band3 = f(as3);
[as4, ~] = Proximate(threshold, P4);
band4 = f(as4 + L2 - 1);
B_ar = abs(band4 - band3);

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
    p = ceil(N/3); % 阶数一般不超出信号长度的1/3
    
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
[h, f] = freqz(1, a, 20e5, Fs);
psdviaBurg = E(end)*abs(h).^2./Fs;
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

% ==================== Proximate ====================
function [as1, as11] = Proximate(b, aa)
%%% b是要找的参数值
%%% aa是要找的数组

%%判断aa是几维数组，统一转换为一维
a = aa(:);  % 将数组转换为一维数组
ab = (a(:) - b)';  % 计算数组a和b的差值
abc = abs(ab);
abc = sort(abc);  % 绝对值取最小值，排序

%%%  二维数组as存储b值最接近的值，在abc中的位置
%    找到与b值最接近的第一个数组元素的位置（最接近b的）
[as1, as11] = find(abs((a(:) - b)) == abc(1, 1));
end

