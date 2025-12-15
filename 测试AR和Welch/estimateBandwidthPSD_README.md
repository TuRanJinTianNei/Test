# estimateBandwidthPSD.m - 带宽估计工具函数说明

## 概述

`estimateBandwidthPSD.m` 是一个独立的工具函数，用于对已处理的OFDM信号进行功率谱密度（PSD）估计和带宽计算。该函数从 `method.m` 中的 `PSD_OFDM_rayleigh` 函数拆分而来，提供了更灵活和可重用的带宽估计功能。

## 功能特点

- **两种估计方法**：
  - Welch算法（周期图平均法）
  - AR模型法（自回归模型，使用Burg算法）

- **自适应阈值**：AR模型法根据SNR自动选择阈值（-6dB/-5dB/-3dB）

- **灵活的参数配置**：可自定义Welch算法和AR模型的参数

- **可选的可视化**：支持绘制PSD图形

## 函数签名

```matlab
[B_welch, B_ar, Pxx_welch, f_welch, Pxx_ar, f_ar] = estimateBandwidthPSD(...
    sig_processed, fs, snr, varargin)
```

## 输入参数

### 必需参数

| 参数 | 类型 | 说明 |
|------|------|------|
| `sig_processed` | 行向量 | 已处理的信号（实数信号）<br>通常已经过：上变频、Rayleigh衰落信道、AWGN噪声 |
| `fs` | 标量 | 采样频率（Hz） |
| `snr` | 标量 | 信噪比（dB），用于AR模型的自适应阈值选择 |

### 可选参数（名称-值对）

| 参数名 | 类型 | 默认值 | 说明 |
|--------|------|--------|------|
| `'plot'` | 逻辑值 | `false` | 是否绘制PSD图 |
| `'welch_window'` | 数值 | `100` | Welch算法窗长度（样本数） |
| `'welch_overlap'` | 数值 | `55` | Welch算法重叠样本数 |
| `'welch_nfft'` | 数值 | `8192` | Welch算法FFT点数 |
| `'ar_criterion'` | 字符串 | `'AIC'` | AR模型阶数选择准则（'AIC' 或 'FPE'） |

## 输出参数

| 参数 | 类型 | 说明 |
|------|------|------|
| `B_welch` | 标量 | Welch算法估计的带宽值（Hz） |
| `B_ar` | 标量 | AR模型法估计的带宽值（Hz） |
| `Pxx_welch` | 列向量 | Welch算法功率谱密度估计值（dB，归一化） |
| `f_welch` | 列向量 | Welch算法对应的频率向量（Hz） |
| `Pxx_ar` | 列向量 | AR模型功率谱密度估计值（dB，归一化） |
| `f_ar` | 列向量 | AR模型对应的频率向量（Hz） |

## 使用示例

### 示例1：基本用法

```matlab
% 假设 sig_processed 是已处理的信号
fs = 40e6;  % 采样频率 40 MHz
snr = 20;   % 信噪比 20 dB

% 估计带宽
[B_welch, B_ar] = estimateBandwidthPSD(sig_processed, fs, snr);

fprintf('Welch算法估计带宽: %.2f MHz\n', B_welch/1e6);
fprintf('AR模型估计带宽: %.2f MHz\n', B_ar/1e6);
```

### 示例2：绘制PSD图

```matlab
% 估计带宽并绘制PSD图
[B_welch, B_ar, Pxx_w, f_w, Pxx_a, f_a] = estimateBandwidthPSD(...
    sig_processed, fs, snr, 'plot', true);
```

### 示例3：自定义参数

```matlab
% 使用自定义的Welch算法参数
[B_welch, B_ar] = estimateBandwidthPSD(sig_processed, fs, snr, ...
    'welch_window', 256, ...
    'welch_overlap', 128, ...
    'welch_nfft', 16384, ...
    'ar_criterion', 'FPE');
```

### 示例4：在method.m中使用

```matlab
% 在PSD_OFDM_rayleigh函数中调用
sig3_chnl = (sig1_chnl.*exp(1j*2*pi*fc/fs*(0:length(sig1_chnl)-1)));
sig3_chnl = real(MUL_RAYLEIGH(sig3_chnl,itau,power,itn,length(itau),length(sig3_chnl),1/fs,fmax,0));
sig3_chnl = awgn(sig3_chnl,snr,'measured');

% 调用工具函数
[B_welch, B_ar] = estimateBandwidthPSD(sig3_chnl, fs, snr, 'plot', k==1);
```

## 带宽计算方法

### Welch算法

- **阈值**：固定使用 **-3dB** 阈值
- **方法**：在功率谱的前半部分和后半部分分别查找最接近-3dB的频率点
- **带宽**：上边界频率 - 下边界频率

### AR模型法

- **阈值**：根据SNR自适应选择
  - SNR > 4 dB：使用 **-6dB** 阈值
  - 0 < SNR ≤ 4 dB：使用 **-5dB** 阈值
  - SNR ≤ 0 dB：使用 **-3dB** 阈值
- **方法**：在功率谱的前半部分和后半部分分别查找最接近阈值的频率点
- **带宽**：上边界频率 - 下边界频率

## 内部函数

该工具函数包含以下内部辅助函数：

1. **Burg**：使用Burg算法实现AR模型功率谱估计
2. **computeARpara**：计算AR模型参数（Burg算法的子函数）
3. **Proximate**：查找数组中与目标值最接近的元素位置

## 依赖关系

- **MATLAB基础环境**
- **信号处理工具箱**（Signal Processing Toolbox）：
  - `pwelch`：Welch算法功率谱估计
  - `hanning`：汉宁窗函数
  - `freqz`：频率响应计算（用于AR模型）

## 注意事项

1. **输入信号**：函数期望输入已处理的实数信号。如果输入复数信号，需要先进行上变频和取实部处理。

2. **采样频率**：确保 `fs` 参数与信号的实际采样频率一致，否则频率轴和带宽计算会出错。

3. **SNR参数**：`snr` 参数主要用于AR模型的自适应阈值选择，不影响Welch算法的计算。

4. **频率分辨率**：
   - Welch算法：`Δf = fs / welch_nfft`
   - AR模型法：`Δf = fs / 200000`（固定200,000个频率点）

5. **归一化**：输出的功率谱密度都已归一化到最大值（峰值在0dB），便于比较和阈值检测。

## 与method.m的关系

- **原函数**：`method.m` 中的 `PSD_OFDM_rayleigh` 函数
- **拆分内容**：功率谱估计和带宽计算部分
- **保留内容**：信号处理部分（上变频、Rayleigh信道、加噪声）仍在 `PSD_OFDM_rayleigh` 中
- **兼容性**：`PSD_OFDM_rayleigh` 函数已更新，优先使用新工具函数，如果工具函数不存在则使用内部实现（向后兼容）

## 更新日志

- **2025.12.10**：创建工具函数，从 `method.m` 中拆分带宽估计功能

