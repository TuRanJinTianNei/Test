# estimateBandwidthPSD.m - 功率谱密度估计工具函数说明

## 概述

`estimateBandwidthPSD.m` 是一个独立的工具函数，用于对已处理的OFDM信号进行功率谱密度（PSD）估计。该函数从 `method.m` 中的 `PSD_OFDM_rayleigh` 函数拆分而来，提供了更灵活和可重用的功率谱密度估计功能。

**注意**：本函数仅进行功率谱密度估计，不包含带宽计算功能。如需带宽估计，请在使用本函数后自行实现带宽计算逻辑。

## 功能特点

- **两种估计方法**：
  - Welch算法（周期图平均法）
  - AR模型法（自回归模型，使用Burg算法）

- **灵活的参数配置**：可自定义Welch算法和AR模型的参数

- **可选的可视化**：支持绘制PSD图形（AR模型、Welch算法和原始信号频谱）

## 函数签名

```matlab
[Pxx_welch, f_welch, Pxx_ar, f_ar] = estimateBandwidthPSD(...
    sig_processed, fs, snr, varargin)
```

## 输入参数

### 必需参数

| 参数 | 类型 | 说明 |
|------|------|------|
| `sig_processed` | 行向量 | 已处理的信号（实数信号）<br>通常已经过：上变频、Rayleigh衰落信道、AWGN噪声 |
| `fs` | 标量 | 采样频率（Hz） |
| `snr` | 标量 | 信噪比（dB），用于绘图标题显示 |

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

% 估计功率谱密度
[Pxx_w, f_w, Pxx_a, f_a] = estimateBandwidthPSD(sig_processed, fs, snr);

% 可以进一步处理PSD数据，例如计算带宽
% [带宽计算代码...]
```

### 示例2：绘制PSD图

```matlab
% 估计功率谱密度并绘制PSD图
[Pxx_w, f_w, Pxx_a, f_a] = estimateBandwidthPSD(...
    sig_processed, fs, snr, 'plot', true);
```

### 示例3：自定义参数

```matlab
% 使用自定义的Welch算法参数
[Pxx_w, f_w, Pxx_a, f_a] = estimateBandwidthPSD(sig_processed, fs, snr, ...
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

% 调用工具函数进行功率谱估计
[Pxx_w, f_w, Pxx_a, f_a] = estimateBandwidthPSD(sig3_chnl, fs, snr, 'plot', k==1);

% 如果需要带宽估计，可以在此基础上实现
% [带宽计算代码...]
```

## 功率谱估计方法

### Welch算法

- **原理**：周期图平均法，通过分段、加窗、FFT和平均来降低功率谱估计的方差
- **特点**：
  - 计算速度快
  - 对噪声有较好的鲁棒性
  - 频率分辨率受窗长度和FFT点数限制
- **输出**：只返回正频率部分（0到fs/2），适用于实信号

### AR模型法

- **原理**：自回归模型，使用Burg算法估计AR模型参数，然后计算功率谱密度
- **特点**：
  - 频率分辨率高（不受窗长度限制）
  - 对短数据序列效果好
  - 可以自适应选择模型阶数（AIC或FPE准则）
- **输出**：返回0到Fs的完整频率范围，频率采样点数为200,000

## 内部函数

该工具函数包含以下内部辅助函数：

1. **Burg**：使用Burg算法实现AR模型功率谱估计
   - 支持AIC或FPE准则自动选择模型阶数
   - 使用`freqz`计算频率响应，频率采样点数为200,000
   
2. **computeARpara**：计算AR模型参数（Burg算法的子函数）
   - 通过最小化前向和后向预测误差的几何平均来估计反射系数
   - 迭代更新AR模型系数和误差功率

## 依赖关系

- **MATLAB基础环境**
- **信号处理工具箱**（Signal Processing Toolbox）：
  - `pwelch`：Welch算法功率谱估计
  - `hanning`：汉宁窗函数
  - `freqz`：频率响应计算（用于AR模型）

## 注意事项

1. **输入信号**：函数期望输入已处理的实数信号。如果输入复数信号，需要先进行上变频和取实部处理。

2. **采样频率**：确保 `fs` 参数与信号的实际采样频率一致，否则频率轴会出错。

3. **SNR参数**：`snr` 参数仅用于绘图标题显示，不影响功率谱估计的计算。

4. **频率分辨率**：
   - Welch算法：`Δf = fs / welch_nfft`（默认8192点，分辨率约4.88 kHz @ 40 MHz）
   - AR模型法：`Δf = fs / 200000`（固定200,000个频率点，分辨率约200 Hz @ 40 MHz）

5. **归一化**：输出的功率谱密度都已归一化到最大值（峰值在0dB），便于比较和分析。

6. **带宽计算**：本函数不包含带宽计算功能。如需带宽估计，可以在获得PSD后使用阈值法或其他方法自行实现。

## 与method.m的关系

- **原函数**：`method.m` 中的 `PSD_OFDM_rayleigh` 函数
- **拆分内容**：功率谱估计部分
- **保留内容**：信号处理部分（上变频、Rayleigh信道、加噪声）仍在 `PSD_OFDM_rayleigh` 中
- **兼容性**：`PSD_OFDM_rayleigh` 函数已更新，优先使用新工具函数，如果工具函数不存在则使用内部实现（向后兼容）

## 绘图功能

当设置 `'plot', true` 时，函数会生成三个图形：

1. **AR模型功率谱密度估计图**
   - 显示AR模型估计的功率谱密度
   - 标题包含SNR和AR模型阶数信息

2. **Welch算法功率谱密度估计图**
   - 显示Welch算法估计的功率谱密度
   - 标题包含SNR信息

3. **OFDM信号频谱图**
   - 显示原始信号的FFT频谱（归一化）
   - 只显示正频率部分（0到fs/2）
   - 使用信号长度作为FFT点数，提供高分辨率频谱

## 算法细节

### Welch算法实现

```matlab
% 功率谱估计
[Pxx2, f1] = pwelch(sig_processed, hanning(welch_window), welch_overlap, welch_nfft, fs);

% 归一化处理
Pxx22 = Pxx2;
Pxx22 = Pxx22 / min(Pxx22);      % 归一化到最小值
Pxx22 = 10*log10(Pxx22);          % 转换为dB
Pxx22 = Pxx22 - max(Pxx22);       % 归一化到最大值（峰值在0dB）
```

### AR模型实现

```matlab
% 使用Burg算法估计AR模型参数
[Pxx1, f, p] = Burg(sig_processed, fs, ar_criterion);

% Burg算法内部：
% 1. 根据AIC或FPE准则自动选择模型阶数
% 2. 使用freqz计算频率响应（200,000个频率点）
% 3. 计算功率谱密度并归一化
```

## 功率谱密度应用

获得功率谱密度后，可以用于：

1. **带宽估计**：使用阈值法（如-3dB、-6dB）在PSD上查找带宽边界
2. **频谱分析**：分析信号的频率分布特性
3. **信号检测**：检测信号的存在和频率范围
4. **性能评估**：比较不同估计方法的性能

## 更新日志

- **2025.12.10**：创建工具函数，从 `method.m` 中拆分功率谱估计功能
- **2025.12.10**：更新文档，准确描述峰值查找的带宽计算方法
- **2025.12.10**：移除带宽计算功能，仅保留功率谱密度估计
