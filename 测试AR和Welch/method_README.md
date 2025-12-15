# method.m - OFDM信号带宽估算方法说明文档

## 概述

`method.m` 实现了两种功率谱密度（PSD）估计方法来估算OFDM信号的带宽：
1. **Welch算法**（周期图平均法）
2. **AR模型法**（自回归模型，使用Burg算法）

本文档详细解释这两种方法的原理、实现过程和带宽计算方式。

---

## 信号处理流程

在带宽估算之前，信号经过以下处理步骤（在`PSD_OFDM_rayleigh`函数中）：

```matlab
% 1. 上变频：将基带信号上变频到载波频率fc
sig3_chnl = (sig1_chnl.*exp(1j*2*pi*fc/fs*(0:length(sig1_chnl)-1)));

% 2. 通过瑞利衰落信道（多径衰落）
sig3_chnl = real(MUL_RAYLEIGH(sig3_chnl,itau,power,itn,length(itau),length(sig3_chnl),1/fs,fmax,0));

% 3. 添加AWGN噪声
sig3_chnl = awgn(sig3_chnl,snr,'measured');
```

处理后的信号`sig3_chnl`是一个实数信号，然后分别使用两种方法进行功率谱估计。

---

## 方法一：Welch算法（周期图平均法）

### 原理

Welch算法是一种改进的周期图方法，通过以下步骤降低功率谱估计的方差：

1. **分段**：将信号分成多个重叠的段
2. **加窗**：对每段信号应用窗函数（减少频谱泄漏）
3. **FFT**：对每段进行FFT变换，得到周期图
4. **平均**：对所有段的周期图进行平均，得到最终的功率谱估计

### 实现代码

```matlab
[Pxx2, f1] = pwelch(sig3_chnl, hanning(100), 55, 4096*2, fs);
```

**参数说明**：
- `sig3_chnl`：输入信号（经过上变频、瑞利衰落和加噪声处理）
- `hanning(100)`：汉宁窗，窗长度为100个样本
- `55`：重叠样本数（overlap），约55%重叠
- `4096*2`：FFT点数（8192点），决定频率分辨率
- `fs`：采样频率（40 MHz）

**输出**：
- `Pxx2`：功率谱密度估计值（线性单位）
- `f1`：对应的频率向量（Hz）

### 功率谱归一化处理

```matlab
Pxx22 = Pxx2;
Pxx22 = Pxx22/min(Pxx22);        % 归一化到最小值
Pxx22 = 10*log10(Pxx22);          % 转换为dB单位
Pxx22 = Pxx22 - max(Pxx22);       % 归一化到最大值（峰值在0dB）
```

### 带宽计算

Welch算法使用**固定阈值法**计算带宽：

```matlab
% 将功率谱分成前后两半
L1 = ceil(length(Pxx22)/2);
P1 = Pxx22(1:L1,1);      % 前半部分（低频）
P2 = Pxx22(L1:end,1);    % 后半部分（高频）

% 在前半部分找到最接近-3dB的点（下边界）
[as1, as11] = Proximate(-3, P1);
band1 = f1(as1);

% 在后半部分找到最接近-3dB的点（上边界）
[as2, as22] = Proximate(-3, P2);
band2 = f1(as2+L1-1);

% 带宽 = 上边界 - 下边界
B_welch = abs(band1 - band2);
```

**关键点**：
- 使用**-3dB阈值**（功率下降一半的位置）
- 分别在前半部分和后半部分查找阈值点
- 带宽 = 上边界频率 - 下边界频率

### 优点
- 计算速度快
- 对噪声有较好的鲁棒性
- 实现简单，MATLAB内置函数

### 缺点
- 频率分辨率受窗长度限制
- 需要选择合适的窗函数和重叠率

---

## 方法二：AR模型法（自回归模型）

### 原理

AR（Auto-Regressive）模型将信号建模为：

\[
x(n) = -\sum_{k=1}^{p} a_k x(n-k) + e(n)
\]

其中：
- \(p\)：模型阶数
- \(a_k\)：AR模型系数
- \(e(n)\)：白噪声激励

通过估计AR模型参数，可以计算功率谱密度：

\[
P_{AR}(f) = \frac{\sigma^2}{|1 + \sum_{k=1}^{p} a_k e^{-j2\pi fk}|^2}
\]

其中\(\sigma^2\)是激励噪声的方差。

### 实现代码

#### 1. 模型阶数选择（AIC准则）

```matlab
[p, f, p] = Burg(sig3_chnl, fs, 'AIC');
```

**AIC（Akaike Information Criterion）准则**：
\[
AIC(p) = N \ln(E_p) + 2p
\]

其中：
- \(N\)：信号长度
- \(E_p\)：p阶AR模型的预测误差功率
- \(p\)：模型阶数

选择使AIC最小的阶数作为最优模型阶数。

#### 2. Burg算法估计AR参数

Burg算法通过最小化前向和后向预测误差的几何平均来估计反射系数：

```matlab
function [a, E] = computeARpara(x, p)
    % 初始化
    ef = x;  % 前向预测误差
    eb = x;  % 后向预测误差
    a = 1;   % AR模型系数（初始为1）
    
    for m = 1:p
        % 计算反射系数
        efm = ef(2:end);
        ebm = eb(1:end-1);
        k(m) = -2 * (ebm' * efm) / (efm' * efm + ebm' * ebm);
        
        % 更新预测误差
        ef = efm + k(m) * ebm;
        eb = ebm + conj(k(m)) * efm;
        
        % 更新AR系数
        a = [a; 0] + k(m) * [0; conj(flipud(a))];
        
        % 更新误差功率
        E(m+1) = (1 - abs(k(m))^2) * E(m);
    end
end
```

#### 3. 功率谱计算

```matlab
[h, f] = freqz(1, a, 20e5, Fs);              % 计算频率响应
psdviaBurg = E(end) * abs(h).^2 / Fs;       % 计算功率谱
psdviaBurg = psdviaBurg / abs(max(psdviaBurg));  % 归一化
psdviaBurg = 10*log10(abs(psdviaBurg));      % 转换为dB
```

**参数说明**：
- `freqz(1, a, 20e5, Fs)`：计算AR模型的频率响应
  - `1`：分子多项式（全通）
  - `a`：分母多项式（AR系数）
  - `20e5`：频率采样点数（200,000点，高分辨率）
  - `Fs`：采样频率

### 带宽计算

AR模型法使用**自适应阈值法**，根据SNR动态选择阈值：

```matlab
% 将功率谱分成前后两半
L2 = ceil(length(Pxx1)/2);
P3 = Pxx1(1:L2,1);      % 前半部分（低频）
P4 = Pxx1(L2:end,1);    % 后半部分（高频）

% 根据SNR选择阈值
if snr > 4
    threshold = -6;      % 高SNR：使用-6dB阈值
elseif snr > 0 && snr <= 4
    threshold = -5;      % 中等SNR：使用-5dB阈值
else
    threshold = -3;      % 低SNR：使用-3dB阈值
end

% 在前半部分找到最接近阈值的点（下边界）
[as3, as33] = Proximate(threshold, P3);
band3 = f(as3);

% 在后半部分找到最接近阈值的点（上边界）
[as4, as44] = Proximate(threshold, P4);
band4 = f(as4+L2-1);

% 带宽 = 上边界 - 下边界
B_ar = abs(band4 - band3);
```

**关键点**：
- **自适应阈值**：根据SNR动态调整（-6dB/-5dB/-3dB）
- SNR越高，使用更严格的阈值（-6dB），可以获得更精确的带宽估计
- SNR越低，使用更宽松的阈值（-3dB），避免噪声影响

### 优点
- 频率分辨率高（不受窗长度限制）
- 对短数据序列效果好
- 可以自适应选择模型阶数

### 缺点
- 计算复杂度较高
- 需要估计模型阶数
- 对模型阶数选择敏感

---

## 两种方法的对比

| 特性 | Welch算法 | AR模型法 |
|------|-----------|----------|
| **频率分辨率** | 受窗长度限制 | 高分辨率（可达200,000点） |
| **计算速度** | 快 | 较慢 |
| **阈值选择** | 固定-3dB | 自适应（-6dB/-5dB/-3dB） |
| **噪声鲁棒性** | 较好（通过平均降低方差） | 中等（依赖模型阶数） |
| **实现复杂度** | 低（MATLAB内置函数） | 高（需要实现Burg算法） |
| **适用场景** | 长信号、实时处理 | 短信号、高精度需求 |

---

## 辅助函数说明

### Proximate函数

用于查找数组中与目标值最接近的元素位置：

```matlab
function [as1, as11] = Proximate(b, aa)
    a = aa(:);                          % 转换为一维数组
    ab = (a(:) - b)';                   % 计算差值
    abc = abs(ab);                      % 取绝对值
    abc = sort(abc);                    % 排序
    [as1, as11] = find(abs((a(:) - b)) == abc(1,1));  % 找到最接近的位置
end
```

**用途**：在功率谱中找到最接近阈值（如-3dB、-5dB、-6dB）的频率点。

---

## 带宽检测率计算

`Bandwidth_rate_rayleigh`函数计算不同SNR下的带宽检测率：

```matlab
B_rate = 1 - abs((B_estimated - B_ideal)) / B_ideal
```

其中：
- `B_estimated`：估计的带宽值（Welch或AR方法）
- `B_ideal`：理想带宽值（8 MHz）

检测率越接近1（100%），表示估计越准确。

---

## 使用示例

```matlab
% 生成OFDM信号
[sig_ofdm] = PSD_generate(20);

% 设置参数
fc = 10e6;      % 载波频率
fs = 40e6;      % 采样频率
snr = 20;       % 信噪比（dB）

% 计算带宽（两种方法）
[B_welch, B_ar] = PSD_OFDM_rayleigh(sig_ofdm, fc, fs, snr, 1, ...
    itau, power, fmax, itn);

% 输出结果
fprintf('Welch算法估计带宽: %.2f MHz\n', B_welch/1e6);
fprintf('AR模型估计带宽: %.2f MHz\n', B_ar/1e6);
```

---

## 参考文献

1. **Welch算法**：
   - Welch, P. D. (1967). "The use of fast Fourier transform for the estimation of power spectra: A method based on time averaging over short, modified periodograms." *IEEE Transactions on Audio and Electroacoustics*, 15(2), 70-73.

2. **AR模型和Burg算法**：
   - Burg, J. P. (1967). "Maximum entropy spectral analysis." *37th Annual International SEG Meeting*, Oklahoma City.
   - Akaike, H. (1974). "A new look at the statistical model identification." *IEEE Transactions on Automatic Control*, 19(6), 716-723.

---

## 注意事项

1. **采样频率一致性**：确保用于功率谱估计的`fs`与信号生成时的采样频率一致。

2. **阈值选择**：
   - Welch算法使用固定-3dB阈值
   - AR模型法根据SNR自适应选择阈值

3. **模型阶数**：AR模型的阶数通过AIC准则自动选择，一般不超过信号长度的1/3。

4. **频率分辨率**：
   - Welch算法：\(\Delta f = \frac{fs}{N_{FFT}}\)
   - AR模型法：\(\Delta f = \frac{fs}{N_{freq}}\)（\(N_{freq} = 200,000\)）

5. **归一化处理**：两种方法都将功率谱归一化到最大值（峰值在0dB），便于比较和阈值检测。

---

## 更新日志

- **2017.03.01**：初始版本，实现Welch算法和AR模型法
- **2025.12.10**：添加详细注释和文档说明

