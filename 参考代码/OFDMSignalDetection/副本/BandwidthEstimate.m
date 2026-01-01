function [B_welch, B_ar] = BandwidthEstimate(signal, fs, snr, method)
%**************************************************************************
%功能：使用Welch算法或AR模型估计信号带宽
%输入参数：
%   signal：输入信号（可以是复数或实数）
%   fs：采样频率（Hz）
%   snr：信噪比（dB），用于AR模型阈值选择
%   method：估计方法选择
%           'welch' - 仅使用Welch算法
%           'ar' - 仅使用AR模型
%           'both' - 同时使用两种方法（默认）
%输出参数：
%   B_welch：Welch算法估计的带宽值（Hz）
%   B_ar：AR模型估计的带宽值（Hz）
%**************************************************************************

% 参数检查
if nargin < 3
    error('至少需要3个输入参数：signal, fs, snr');
end
if nargin < 4
    method = 'both';  % 默认同时使用两种方法
end

% 确保信号是列向量
signal = signal(:);

% 如果输入是复数信号，转换为实数
if ~isreal(signal)
    signal = real(signal);
end

% 初始化输出
B_welch = [];
B_ar = [];

% ==================== Welch算法带宽估计 ====================
if strcmpi(method, 'welch') || strcmpi(method, 'both')
    % 使用Welch算法进行功率谱估计
    % 参数：hanning窗口长度100，重叠55，FFT点数4096*2
    [Pxx2, f1] = pwelch(signal, hanning(100), 55, 4096*2, fs);
    
    % 归一化处理
    Pxx22 = Pxx2;
    Pxx22 = Pxx22 / min(Pxx22);  % 归一化到最小值
    Pxx22 = 10 * log10(Pxx22);   % 转换为dB单位
    Pxx22 = Pxx22 - max(Pxx22);  % 归一化到最大值
    
    % 计算带宽：找到-3dB点
    L1 = ceil(length(Pxx22) / 2);
    P1 = Pxx22(1:L1, 1);
    P2 = Pxx22(L1:end, 1);
    
    % 找到最接近-3dB的点
    [as1, ~] = Proximate(-3, P1);
    band1 = f1(as1);
    [as2, ~] = Proximate(-3, P2);
    band2 = f1(as2 + L1 - 1);
    
    B_welch = abs(band1 - band2);
end

% ==================== AR模型带宽估计 ====================
if strcmpi(method, 'ar') || strcmpi(method, 'both')
    % 使用Burg算法（AR模型）进行功率谱估计
    % 使用AIC准则自动确定模型阶数
    [Pxx1, f, ~] = Burg(signal, fs, 'AIC');
    
    % 计算带宽：根据SNR选择不同的阈值
    L2 = ceil(length(Pxx1) / 2);
    P3 = Pxx1(1:L2, 1);
    P4 = Pxx1(L2:end, 1);
    
    % 根据SNR选择阈值
    if snr > 4
        threshold = -5;  % 高SNR使用-5dB阈值
    elseif snr > 0 && snr <= 4
        threshold = -4;  % 中等SNR使用-4dB阈值
    else
        threshold = -3;  % 低SNR使用-3dB阈值
    end
    
    % 找到最接近阈值的点
    [as3, ~] = Proximate(threshold, P3);
    band3 = f(as3);
    [as4, ~] = Proximate(threshold, P4);
    band4 = f(as4 + L2 - 1);
    
    B_ar = abs(band4 - band3);
end

end

