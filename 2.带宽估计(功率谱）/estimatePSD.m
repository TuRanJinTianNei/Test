function [Pxx_welch, f_welch, Pxx_ar, f_ar] = estimatePSD(sig_processed, fs, snr, varargin)
%===============================================================================
% estimatePSD.m - 使用Welch算法和AR模型法估计信号功率谱密度
% 
% 功能说明:
%   对已处理的信号进行功率谱密度估计，特别适用于signalGenerate.m生成的接收信号
%   使用两种方法：
%   1. Welch算法（周期图平均法）
%   2. AR模型法（自回归模型，使用Burg算法）
%
% 输入参数:
%   sig_processed  - 输入信号（实数信号，行向量或列向量）
%                    可以是：
%                    - signalGenerate.m生成的Rx_data（推荐）
%                    - signalGenerate.m生成的Tx_data
%                    - 其他已处理的OFDM基带信号
%                    注意：信号应为实数基带信号，已通过共轭对称映射
%   fs             - 采样频率（Hz）
%                    对于signalGenerate.m：fs = 3.84e6 Hz（LTE 3MHz标准）
%   snr            - 信噪比（dB），用于绘图标题显示
%                    对于signalGenerate.m：使用targetSNRdB变量
%   varargin       - 可选参数（名称-值对）:
%                    'plot'      - 是否绘制PSD图，默认false
%                    'welch_window' - Welch算法窗长度，默认100
%                    'welch_overlap' - Welch算法重叠样本数，默认55
%                    'welch_nfft' - Welch算法FFT点数，默认8192
%                    'ar_criterion' - AR模型阶数选择准则，默认'AIC'（可选'AIC'或'FPE'）
%
% 输出参数:
%   Pxx_welch      - Welch算法功率谱密度估计值（dB，归一化，峰值在0dB）
%   f_welch        - Welch算法对应的频率向量（Hz，单边谱：0到fs/2）
%   Pxx_ar         - AR模型功率谱密度估计值（dB，归一化，峰值在0dB）
%   f_ar           - AR模型对应的频率向量（Hz，单边谱：0到fs/2）
%
% 使用示例:
%   % 示例1：独立运行（自动生成信号）
%   estimatePSD  % 将自动调用signalGenerate.m生成信号，然后进行PSD估计并绘图
%
%   % 示例2：使用signalGenerate.m生成的接收信号
%   signalGenerate;  % 生成Tx_data和Rx_data
%   [Pxx_w, f_w, Pxx_a, f_a] = estimatePSD(Rx_data, fs, targetSNRdB);
%
%   % 示例3：绘制PSD图
%   [Pxx_w, f_w, Pxx_a, f_a] = estimatePSD(...
%       Rx_data, fs, targetSNRdB, 'plot', true);
%
%   % 示例4：自定义Welch算法参数
%   [Pxx_w, f_w, Pxx_a, f_a] = estimatePSD(...
%       Rx_data, fs, targetSNRdB, 'welch_window', 200, 'welch_nfft', 16384);
%
%   % 示例5：使用AR模型的FPE准则
%   [Pxx_w, f_w, Pxx_a, f_a] = estimatePSD(...
%       Rx_data, fs, targetSNRdB, 'ar_criterion', 'FPE');
%
% 创建日期: 2025.12.10
% 基于: method.m 中的 PSD_OFDM_rayleigh 函数
% 修改日期: 
%   2025.12.10 - 移除带宽计算功能，仅保留功率谱密度估计
%   2025.12.23 - 优化以适配signalGenerate.m的输出，添加输入验证和更好的注释
%   2025.12.23 - 添加独立运行模式，无输入参数时自动生成信号
%   2025.12.23 - 文件重命名，文件名与函数名一致
%===============================================================================

%**************************************************************************
% 独立运行模式：如果没有输入参数，自动生成信号
%**************************************************************************
if nargin == 0
    fprintf('========================================\n');
    fprintf('功率谱密度估计程序（独立运行模式）\n');
    fprintf('========================================\n\n');
    
    % 检查工作空间中是否已有信号数据
    if exist('Rx_data', 'var') && exist('fs', 'var') && exist('targetSNRdB', 'var')
        fprintf('检测到工作空间中已有信号数据，直接使用...\n');
        fprintf('  - 使用变量: Rx_data\n');
        fprintf('  - 采样频率: %.2f MHz\n', fs/1e6);
        fprintf('  - SNR: %.1f dB\n\n', targetSNRdB);
        sig_processed = Rx_data;
        snr = targetSNRdB;
    elseif exist('Tx_data', 'var') && exist('fs', 'var')
        fprintf('检测到工作空间中已有发送信号，使用Tx_data...\n');
        fprintf('  - 使用变量: Tx_data\n');
        fprintf('  - 采样频率: %.2f MHz\n', fs/1e6);
        if exist('targetSNRdB', 'var')
            snr = targetSNRdB;
        else
            snr = 20;  % 默认SNR
            fprintf('  - SNR: %.1f dB（默认值，未找到targetSNRdB）\n', snr);
        end
        fprintf('\n');
        sig_processed = Tx_data;
    else
        fprintf('未找到信号数据，正在调用signalGenerate.m生成信号...\n\n');
        % 调用signalGenerate.m生成信号
        signalGenerate;
        
        % 检查是否生成了接收信号
        if exist('Rx_data', 'var')
            fprintf('\n使用生成的接收信号（Rx_data）进行PSD估计...\n');
            sig_processed = Rx_data;
            snr = targetSNRdB;
        elseif exist('Tx_data', 'var')
            fprintf('\n使用生成的发送信号（Tx_data）进行PSD估计...\n');
            sig_processed = Tx_data;
            if exist('targetSNRdB', 'var')
                snr = targetSNRdB;
            else
                snr = 20;  % 默认SNR
            end
        else
            error('错误：signalGenerate.m未能生成信号数据');
        end
        fprintf('  采样频率: %.2f MHz\n', fs/1e6);
        fprintf('  SNR: %.1f dB\n\n', snr);
    end
    
    % 独立运行模式默认启用绘图
    varargin = {'plot', true};
    fprintf('独立运行模式：默认启用绘图\n\n');
end

% 输入验证（函数调用模式）
if nargin > 0 && nargin < 3
    error('错误：函数调用模式需要至少3个输入参数：sig_processed, fs, snr\n或者无参数运行以进入独立运行模式');
end

% 确保信号是行向量
if exist('sig_processed', 'var') && iscolumn(sig_processed)
    sig_processed = sig_processed';
end

% 验证信号长度
if length(sig_processed) < 100
    warning('警告：信号长度过短（<100样本），可能影响PSD估计精度');
end

% 验证采样频率
if fs <= 0
    error('错误：采样频率fs必须大于0');
end

% 解析可选参数
p = inputParser;
addParameter(p, 'plot', false, @islogical);
addParameter(p, 'welch_window', 100, @isnumeric);
addParameter(p, 'welch_overlap', 55, @isnumeric);
addParameter(p, 'welch_nfft', 8192, @isnumeric);
addParameter(p, 'ar_criterion', 'AIC', @(x) ischar(x) && (strcmpi(x, 'AIC') || strcmpi(x, 'FPE')));
parse(p, varargin{:});

plot_flag = p.Results.plot;
welch_window = p.Results.welch_window;
welch_overlap = p.Results.welch_overlap;
welch_nfft = p.Results.welch_nfft;
ar_criterion = p.Results.ar_criterion;

% 验证Welch算法参数
if welch_window <= 0 || welch_overlap < 0 || welch_overlap >= welch_window
    error('错误：Welch算法参数无效。要求：welch_window > 0, 0 <= welch_overlap < welch_window');
end

if welch_nfft < welch_window
    warning('警告：FFT点数（%d）小于窗长度（%d），建议welch_nfft >= welch_window', ...
        welch_nfft, welch_window);
end

% 显示处理信息（可选，用于调试）
if plot_flag
    fprintf('功率谱密度估计参数：\n');
    fprintf('  - 信号长度: %d 样本\n', length(sig_processed));
    fprintf('  - 采样频率: %.2f MHz\n', fs/1e6);
    fprintf('  - SNR: %.1f dB\n', snr);
    fprintf('  - Welch算法：窗长度=%d, 重叠=%d, FFT点数=%d\n', ...
        welch_window, welch_overlap, welch_nfft);
    fprintf('  - AR模型：阶数选择准则=%s\n', ar_criterion);
    fprintf('开始计算功率谱密度...\n');
end

%**************************************************************************
% Welch算法功率谱估计
%**************************************************************************
% 注意：对于实信号，频谱是共轭对称的，即 X(-f) = X*(f)
%      pwelch函数对实信号自动返回单边谱（0到fs/2），并且已经自动将
%      正频率部分的幅度乘以2（除了DC和Nyquist频率），所以输出已经是单边谱
%      单边谱幅度 = 2 × 双边谱幅度（对于0 < f < fs/2）
%      DC（f=0）和Nyquist频率（f=fs/2）的幅度保持不变
%
% 对于signalGenerate.m生成的信号：
%   - Rx_data或Tx_data都是实数基带信号（通过共轭对称映射实现）
%   - 可以直接使用pwelch进行功率谱估计
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

% 归一化处理（在转换为dB之前）
Pxx22 = Pxx2;
% 归一化到最小值（添加小值保护，避免除以0或极小值）
min_val = min(Pxx22);
if min_val <= 0 || min_val < max(Pxx22) * 1e-10
    % 如果最小值太小或为0，使用最大值归一化
    Pxx22 = Pxx22 / max(Pxx22);
else
    Pxx22 = Pxx22 / min_val;  % 归一化到最小值
end
Pxx22 = 10*log10(Pxx22 + eps);     % 转换为dB单位（添加eps避免log10(0)）
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
%
% 对于signalGenerate.m生成的信号：
%   - AR模型法可以提供更高的频率分辨率
%   - 特别适合短数据记录或需要精细频谱细节的场景
[Pxx1, f, p] = Burg(sig_processed, fs, ar_criterion);

% 输出AR功率谱（正频率部分：0到fs/2）
Pxx_ar = Pxx1;
f_ar = f;

% 显示完成信息（如果启用了绘图）
if plot_flag
    fprintf('功率谱密度计算完成！\n');
    fprintf('  - Welch算法：频率点数 = %d，频率范围 = %.3f MHz 到 %.3f MHz\n', ...
        length(f_welch), min(f_welch)/1e6, max(f_welch)/1e6);
    fprintf('  - AR模型法：频率点数 = %d，频率范围 = %.3f MHz 到 %.3f MHz，模型阶数 = %d\n', ...
        length(f_ar), min(f_ar)/1e6, max(f_ar)/1e6, p);
end

%**************************************************************************
% 绘制PSD图形（可选）
%**************************************************************************
% 如果启用绘图，将显示：
%   1. AR模型功率谱密度估计（单边谱）
%   2. Welch算法功率谱密度估计（单边谱）
%   3. OFDM信号频谱（单边谱，使用FFT直接计算）
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

%**************************************************************************
% 独立运行模式：显示总结信息
%**************************************************************************
if nargin == 0
    fprintf('\n========================================\n');
    fprintf('功率谱密度估计完成！\n');
    fprintf('========================================\n');
    fprintf('结果已保存在以下变量中：\n');
    fprintf('  - Pxx_welch: Welch算法PSD估计（dB，归一化）\n');
    fprintf('  - f_welch: Welch算法频率向量（Hz）\n');
    fprintf('  - Pxx_ar: AR模型PSD估计（dB，归一化）\n');
    fprintf('  - f_ar: AR模型频率向量（Hz）\n');
    fprintf('\nPSD估计参数：\n');
    fprintf('  - Welch算法：频率点数 = %d，频率范围 = %.3f MHz 到 %.3f MHz\n', ...
        length(f_welch), min(f_welch)/1e6, max(f_welch)/1e6);
    fprintf('  - AR模型法：频率点数 = %d，频率范围 = %.3f MHz 到 %.3f MHz，模型阶数 = %d\n', ...
        length(f_ar), min(f_ar)/1e6, max(f_ar)/1e6, p);
    fprintf('\n提示：可以使用以下命令进行带宽估计：\n');
    fprintf('  [B_estimated, results] = estimateBandwidthMethods(Pxx_welch, f_welch, ''method'', ''all'');\n');
    fprintf('========================================\n');
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
    
    % 计算0到p阶的误差（用于选择最优阶数）
    % 注意：computeARpara返回的E长度为p+1，E(1)对应0阶，E(2)对应1阶，...，E(p+1)对应p阶
    [~, E] = computeARpara(x, p);
    
    % 计算目标函数的最小值
    % kc是阶数索引：0到p，对应E(1)到E(p+1)
    kc = 0:p;  % 阶数从0到p
    switch criterion
        case 'FPE'
            % FPE准则：FPE(k) = E(k+1) * (N + k + 1) / (N - k - 1)
            % 其中k是阶数（0到p），E(k+1)是k阶对应的误差
            goalF = E(kc+1).*(N + (kc + 1))./(N - (kc + 1));
        case 'AIC'
            % AIC准则：AIC(k) = N*log(E(k+1)) + 2*k
            % 其中k是阶数（0到p），E(k+1)是k阶对应的误差
            goalF = N.*log(E(kc+1)) + 2.*kc;
    end
    [~, p_idx] = min(goalF); % p_idx是目标函数最小值的索引位置
    p = kc(p_idx);  % 对应的阶数值（0到p之间的某个值）
    
    % 使用p值重新计算AR模型参数
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

