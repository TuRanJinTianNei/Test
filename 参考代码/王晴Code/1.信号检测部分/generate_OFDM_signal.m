function [TrData, ReData] = generate_OFDM_signal(SubCarryNN, SubCarryN, SymbN, ratio, SNR)
% 生成OFDM信号
% 
% 输入参数：
%   SubCarryNN  - 有效子载波数（默认：256）
%   SubCarryN   - 总子载波数/FFT长度（默认：1024）
%   SymbN       - OFDM符号个数（默认：20）
%   ratio       - 循环前缀比例（默认：1/4）
%   SNR         - 信噪比，单位dB（可选，默认不添加噪声）
%
% 输出参数：
%   TrData      - 发射信号（时域，复数）
%   ReData      - 接收信号（添加噪声后，如果SNR未指定则等于TrData）
%
% 示例：
%   [TrData, ReData] = generate_OFDM_signal(256, 1024, 20, 1/4, 10);
%   [TrData, ReData] = generate_OFDM_signal();  % 使用默认参数

% 设置默认参数
if nargin < 1 || isempty(SubCarryNN)
    SubCarryNN = 256;      % 有效子载波数
end
if nargin < 2 || isempty(SubCarryN)
    SubCarryN = 1024;      % 总子载波数
end
if nargin < 3 || isempty(SymbN)
    SymbN = 20;            % OFDM符号个数
end
if nargin < 4 || isempty(ratio)
    ratio = 1/4;           % 循环前缀比例
end

fftLen = SubCarryN;       % FFT长度
GuardLen = SubCarryN * ratio;    % 保护时隙的长度

%% 1. 生成随机比特流
SignalLen = SubCarryNN * SymbN * 2;           % 比特序列长度
Signal = round(rand(1, SignalLen));           % 生成随机二进制比特流

%% 2. 数据重组
ParaBitSig = reshape(Signal, SubCarryNN, SymbN * 2);

%% 3. QPSK调制
% 将数据分为I（同相）和Q（正交）两个通道
ich = zeros(SubCarryNN, SymbN);
qch = zeros(SubCarryNN, SymbN);
for j = 1:SymbN
    ich(:, j) = ParaBitSig(:, 2*j-1);  % 同相分量，奇数列
    qch(:, j) = ParaBitSig(:, 2*j);    % 正交分量，偶数列
end

% QPSK映射：0→-1, 1→1，然后归一化
kmod = 1 ./ sqrt(2);
ich0 = ich .* 2 - 1;      % 0→-1, 1→1
qch0 = qch .* 2 - 1;
ich1 = ich0 .* kmod;      % 归一化
qch1 = qch0 .* kmod;

% 形成复数信号
x = ich1 + qch1 .* sqrt(-1);

%% 4. 频域映射（将数据映射到子载波位置）
% 将有效子载波数据映射到FFT的对应位置
x = [x(1:SubCarryNN/2, 1:SymbN); ...
     zeros(SubCarryN - SubCarryNN, SymbN); ...
     x(SubCarryNN/2+1:SubCarryNN, 1:SymbN)];

%% 5. IFFT变换（频域转时域）
y = ifft(x, fftLen);      % 通过IFFT将频域数据转化为时域数据
ich2 = real(y);           % I信道取实部
qch2 = imag(y);           % Q信道取虚部

%% 6. 添加循环前缀（保护间隔）
ich3 = [ich2(fftLen - GuardLen + 1:fftLen, :); ich2];
qch3 = [qch2(fftLen - GuardLen + 1:fftLen, :); qch2];

%% 7. 并串变换
ich4 = reshape(ich3, 1, (fftLen + GuardLen) * SymbN);
qch4 = reshape(qch3, 1, (fftLen + GuardLen) * SymbN);

%% 8. 形成复数发射信号
TrData = ich4 + qch4 .* sqrt(-1);

%% 9. 添加噪声（如果指定了SNR）
if nargin >= 5 && ~isempty(SNR)
    ReData = awgn(TrData, SNR, 'measured');  % 添加高斯白噪声
else
    ReData = TrData;  % 不添加噪声
end

end

