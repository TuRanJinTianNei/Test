# 短波OFDM信号侦测与参数估计 

[![MATLAB](https://img.shields.io/badge/MATLAB-R2014b+-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

一个基于MATLAB的**短波OFDM信号侦测与参数估计**系统，用于OFDM信号的调制识别和关键参数盲估计。项目支持在**高斯信道（AWGN）**和**多径瑞利衰落信道**两种环境下进行仿真分析。

## 📋 目录

- [项目概述](#项目概述)
- [主要功能](#主要功能)
- [快速开始](#快速开始)
- [项目结构](#项目结构)
- [环境要求](#环境要求)
- [使用示例](#使用示例)
- [文档说明](#文档说明)
- [参考资源](#参考资源)
- [作者信息](#作者信息)

## 🎯 项目概述

本项目实现了完整的OFDM信号检测与参数估计流程，包括：

1. **OFDM调制识别**：区分OFDM信号与单载波调制信号（QPSK、QAM16、QAM64、FSK8等）
2. **OFDM参数估计**：盲估计OFDM信号的各项关键参数

项目采用模块化设计，每个功能模块都有独立的函数实现，便于使用和维护。

## ✨ 主要功能

### 1. 信号调制识别

- **基于高阶累积量的识别方法**
  - 利用OFDM信号的高斯特性（高阶累积量C42≈0）来区分OFDM信号和单载波信号
  - 在信噪比大于10dB时效果较好

- **基于小波变换的识别方法**
  - 利用小波变换的尺度参数和方差特征来识别信号类型
  - 对信号特征提取更稳定

### 2. 参数估计功能

#### 2.1 带宽估计

**算法流程**：

1. **功率谱密度（PSD）估计**
   - **Welch算法**：
     - 使用`pwelch`函数，窗函数为Hanning窗（窗长100，重叠55）
     - FFT点数为8192（4096×2）
     - 对PSD进行归一化处理：`Pxx = Pxx/min(Pxx)`，转换为dB单位，峰值归一化到0dB
   - **AR模型法（Burg算法）**：
     - 使用`Burg`函数估计AR模型参数
     - 模型阶数通过AIC准则自动选择
     - 计算AR模型的频率响应，得到PSD估计

2. **带宽提取**
   - **Welch算法**：
     - 将PSD分为前后两半（以中心频率为界）
     - 使用`Proximate`函数在前半部分找到最接近-3dB的频率点`band1`
     - 在后半部分找到最接近-3dB的频率点`band2`
     - 带宽 = `|band1 - band2|`
   - **AR模型法**：
     - 同样将PSD分为前后两半
     - 根据SNR自适应选择阈值：
       - SNR > 4dB：使用-5dB阈值（AWGN）或-6dB阈值（Rayleigh）
       - 0dB < SNR ≤ 4dB：使用-4dB阈值（AWGN）或-5dB阈值（Rayleigh）
       - SNR ≤ 0dB：使用-3dB阈值
     - 找到对应的频率点，计算带宽

3. **性能评估**
   - 通过蒙特卡洛仿真（默认10次）计算带宽检测准确率
   - 准确率 = `1 - |估计值 - 理想值| / 理想值`

**相关函数**：
- `PSD_OFDM.m`：AWGN信道下的带宽估计
- `PSD_OFDM_rayleigh.m`：Rayleigh信道下的带宽估计
- `Bandwidth_rate.m`：带宽检测准确率计算（AWGN）
- `Bandwidth_rate_rayleigh.m`：带宽检测准确率计算（Rayleigh）

---

#### 2.2 载频估计

**算法流程**：

1. **循环谱计算**
   - 使用平滑周期图法计算循环谱`S(α, f)`
   - 循环频率分辨率：`d_alpha = fs/N`（N为数据长度）
   - 循环频率范围：`-fs`到`fs`
   - 频率分辨率：`d_alpha × M`（M为平滑长度）

2. **循环自相关函数提取**
   - 提取循环频率α=0的二维切片：`S(0, f)`
   - 对切片进行插值处理（插值倍数100），提高频率分辨率

3. **载频估计**
   - 在频率域前半部分找到最大峰值位置`f1`
   - 在频率域后半部分找到最大峰值位置`f2`
   - 载频估计值：`fc = (f1 + f2) / 2`

**原理**：
  - 利用OFDM信号的循环平稳特性
- 在循环频率α=0处，循环谱的峰值位置对应载频

**相关函数**：
- `cyclic_spectrum.m`：循环谱计算和载频估计

---

#### 2.3 有效数据长度、符号总长度、循环前缀长度估计

**算法流程**：

1. **有效数据长度（Tu）估计**
   - **可变延时自相关计算**：
     - 计算信号的可变延时自相关函数`Rx(n, τ)`
     - 延时范围：`-P`到`P`（P为符号长度采样点数）
     - 对所有时间点n的自相关函数进行平均：`R_tao = mean(Rx(n, τ))`
   - **峰值检测**：
     - 对自相关函数前半部分进行排序
     - 找到第二高的峰值位置`Tu1`
     - 有效数据长度：`Tu = (length(R_tao)/2 - Tu1) × t`（t为采样间隔）

2. **符号总长度（Ts）估计**
   - **固定延时自相关计算**：
     - 固定延时为`length(R_tao)/2 - Tu1`（即有效数据长度对应的采样点数）
     - 计算固定延时自相关函数`RR(n, k)`，k的范围为`1:L_fft`（L_fft = 6×P）
   - **FFT变换**：
     - 对所有时间点n的自相关函数进行平均：`RRR_tao = mean(RR(n, k))`
     - 对`RRR_tao`进行IFFT变换：`RRR_fft = abs(IFFT(RRR_tao))`
   - **峰值间隔检测**：
     - 找到IFFT结果中的前三个最高峰值位置`T1, T2, T3`
     - 计算峰值间隔：`Ts1 = |T2-T1|`, `Ts2 = |T3-T2|`, `Ts3 = |T3-T1|`
     - 取最小间隔`TT = min([Ts1, Ts2, Ts3])`
     - 符号总长度：`Ts = L_fft / TT × t`

3. **循环前缀长度（Tg）估计**
   - 循环前缀长度：`Tg = Ts - Tu`

**原理**：
- OFDM信号的循环前缀导致自相关函数在特定延时处出现峰值
- 可变延时自相关的峰值位置对应有效数据长度
- 固定延时自相关的FFT峰值间隔对应符号总长度

**相关函数**：
- `effectivelength.m`：有效长度估计（AWGN）
- `effectivelength_rayleigh.m`：有效长度估计（Rayleigh）
- `auto_xcorr.m`：自相关函数计算

---

#### 2.4 子载波数目估计

**算法流程**：

1. **信号过采样**
   - 对输入信号进行3倍过采样（`fff = 3`）
   - 通过矩阵复制和重塑实现：`sig_2 = reshape(repmat(yd, 3, 1), 1, 3×length(yd))`

2. **自相关函数计算**
   - 计算自相关函数`Rxx(n, 1)`（固定延时为1）
   - 对每个时间点n，计算：`Rx_c(n, τ) = sig_2(n) × conj(sig_2(n+τ))`
   - 延时范围：`τ = 1:P`（P为符号长度采样点数）

3. **FFT变换和峰值检测**
   - 对`Rx_c(:, 1)`进行IFFT变换（FFT点数6000）：`RRR_c = IFFT(Rx_c(:, 1))`
   - 去除前后边界点，分析中间部分：`Ln = fft_num/fff - 100`
   - 在三个区间内分别找到最高峰值：
     - 区间1：`[Ln+1, 2×Ln]`，最高峰值位置`a1`
     - 区间2：`[2×Ln+1, 3×Ln]`，最高峰值位置`a3`
   - 计算峰值间隔：`interval = a3 + Ln - a1`
   - 采样倍数：`q = fft_num / interval`

4. **过采样有效数据长度估计**
   - 使用`overfind_num`函数估计过采样信号的有效数据长度`Tu_over`
   - 该函数通过可变延时自相关函数的峰值位置估计有效数据长度

5. **子载波数目计算**
   - 子载波数目：`carry_num = Tu_over / q`

**原理**：
- 过采样后的信号自相关函数FFT结果中，峰值间隔对应采样倍数
- 结合过采样有效数据长度，可以反推原始信号的子载波数目

**相关函数**：
- `carrier_number.m`：子载波数目估计（AWGN）
- `carrier_number_rayleigh.m`：子载波数目估计（Rayleigh）
- `overfind_num.m`：过采样信号的有效数据长度估计

## 🚀 快速开始

### 环境要求

- MATLAB R2014b 或更高版本
- Signal Processing Toolbox（信号处理工具箱）
- Wavelet Toolbox（小波工具箱，用于小波变换识别）

### 安装步骤

1. **克隆或下载项目**
   ```bash
   git clone <repository-url>
   cd OFDMSignalDetection
   ```

2. **添加MATLAB路径**
   ```matlab
   addpath(genpath('OFDM_Signal_Detection'))
   ```

3. **运行主程序**
   ```matlab
   cd OFDM_Signal_Detection
   OFDM_DEMO
   ```

### 基本使用

主程序 `OFDM_DEMO.m` 是项目的入口文件。当前版本主要演示**带宽估计**功能：

```matlab
% 运行主程序
OFDM_DEMO
```

程序会自动：
- 生成OFDM测试信号
- 计算AWGN信道下的带宽检测率
- 计算Rayleigh信道下的带宽检测率
- 生成功率谱密度（PSD）图形
- 绘制性能对比图

## 📁 项目结构

```
OFDMSignalDetection/
├── README.md                    # 项目说明文档（本文件）
├── OFDM_DEMO.md                # OFDM_DEMO.m 运行说明
├── 使用说明.md                 # 完整功能使用说明
└── OFDM_Signal_Detection/       # 主程序目录
    ├── OFDM_DEMO.m             # 主程序入口
    │
    ├── 信号生成模块/
    │   ├── ofdm.m              # 生成OFDM信号
    │   ├── QPSK.m              # 生成QPSK信号
    │   ├── QAM16.m             # 生成16QAM信号
    │   ├── QAM64.m             # 生成64QAM信号
    │   ├── FSK8.m              # 生成8FSK信号
    │   ├── ASK8.m              # 生成8ASK信号
    │   └── random_binary.m     # 生成随机二进制序列
    │
    ├── 识别算法模块/
    │   ├── cumulant.m          # 高阶累积量计算
    │   ├── sc_ofdm_wavelet.m   # 小波变换识别（AWGN）
    │   ├── sc_ofdm_wavelet_Ray.m  # 小波变换识别（Rayleigh）
    │   ├── OFDM_rate.m         # OFDM识别准确率（高阶累积量法）
    │   └── OFDM_rate_xiaobo.m  # OFDM识别准确率（改进方法）
    │
    ├── 参数估计模块/
    │   ├── PSD_OFDM.m          # 功率谱密度估计（AWGN）
    │   ├── PSD_OFDM_rayleigh.m # 功率谱密度估计（Rayleigh）
    │   ├── Bandwidth_rate.m    # 带宽估计准确率（AWGN）
    │   ├── Bandwidth_rate_rayleigh.m  # 带宽估计准确率（Rayleigh）
    │   ├── cyclic_spectrum.m   # 循环谱分析
    │   ├── effectivelength.m   # 有效长度估计（AWGN）
    │   ├── effectivelength_rayleigh.m  # 有效长度估计（Rayleigh）
    │   ├── carrier_number.m    # 子载波数目估计（AWGN）
    │   └── carrier_number_rayleigh.m   # 子载波数目估计（Rayleigh）
    │
    └── 辅助函数模块/
        ├── auto_xcorr.m        # 自相关函数计算
        ├── Burg.m              # Burg算法（AR模型参数估计）
        ├── MUL_RAYLEIGH.m      # 多径瑞利信道模拟
        └── Proximate.m         # 近似计算
```

## 💡 使用示例

### 示例1：OFDM信号识别（高阶累积量法）

```matlab
K = 1;  % K=1时显示图形
N = 20;
para = 128;
ratio = 1/4;
snr1 = -5:5:25;

[c42_ofdm, c42_qpsk, c42_qam16, c42_qam64, c42_fsk8] = ...
    cumulant(snr1, N, para, ratio, K);
```

### 示例2：带宽估计

```matlab
N = 20;
fc = 10e6;  % 载频 10MHz
fs = 40e6;  % 采样频率 40MHz
snr = -4:2:10;

% 生成OFDM信号
[sig_ofdm] = PSD_generate(N);

% 计算带宽检测率
[B_rate_welch, B_rate_ar] = Bandwidth_rate(sig_ofdm, fc, fs, snr);
```

### 示例3：有效数据长度估计

```matlab
N = 20;
para = 128;
ratio = 1/4;
fc = 10e6;
fs = 40e6;
K = 1;

% 生成OFDM信号并上变频
sig = ofdm(N, para, ratio);
sig_a = sig .* exp(1i*2*pi*fc/fs*(0:length(sig)-1));

% 估计有效数据长度、符号总长度、循环前缀长度
[Tu, Ts, Tg] = effectivelength(sig_a, fs, 20, N, K);
```

### 示例4：Rayleigh信道下的分析

```matlab
% 设置Rayleigh信道参数
itau = floor([0,1e-8,2e-8,5e-8,2e-7,5e-7].*fs);
power = [0,-1.0,-7.0,-10.0,-12.0,-17.0];
fmax = 20;
itn = [10000,20000,30000,40000,50000,60000];

% Rayleigh信道下的带宽估计
[B_rate_welch_rayleigh, B_rate_ar_rayleigh] = ...
    Bandwidth_rate_rayleigh(sig_ofdm, fc, fs, snr, itau, power, fmax, itn);
```

## 📚 文档说明

项目包含以下文档：

- **`README.md`**（本文件）：项目概述和快速开始指南
- **`OFDM_DEMO.md`**：`OFDM_DEMO.m` 主程序的详细运行说明，包括当前激活功能的完整解释
- **`使用说明.md`**：项目所有功能的完整使用说明，包括原理解释、参数说明和使用示例

**文档区别**：
- `OFDM_DEMO.md` 专门解释 `OFDM_DEMO.m` 中**当前激活的功能**
- `使用说明.md` 包含项目**所有功能的完整说明**，包括被注释掉的功能

## 📖 参考资源

项目相关的技术博客：

- [短波信道模型--多径瑞利信道原理详解及matlab实现](https://blog.csdn.net/WilsonSong1024/article/details/79449425)
- [用Burg法估计AR模型的参数原理详解及matlab实现](https://blog.csdn.net/WilsonSong1024/article/details/79449161)
- [OFDM信号循环谱原理详解及matlab实现](https://blog.csdn.net/WilsonSong1024/article/details/79449213)
- [OFDM信号的循环自相关函数原理详解及matlab实现](https://blog.csdn.net/WilsonSong1024/article/details/79449318)

## ⚙️ 参数说明

### 常用参数

- **N**：OFDM符号数量（默认：20）
- **para**：子载波数目（常用值：64, 128, 256, 512）
- **ratio**：循环前缀比例（常用值：1/8, 1/4, 3/16）
- **snr**：信噪比范围（dB），例如：`-4:2:10` 表示从-4dB到10dB，步长为2dB
- **fc**：载频（Hz），例如：`10e6` 表示10MHz
- **fs**：采样频率（Hz），通常为载频的2-4倍，例如：`40e6` 表示40MHz
- **K**：图形显示标志，`K=1` 时显示图形，`K=0` 时不显示

### Rayleigh信道参数

- **itau**：多径延时数组（样本数）
- **power**：多径功率数组（dB）
- **fmax**：最大多普勒频率（Hz）
- **itn**：瑞利信道记录次数数组

## ⚠️ 注意事项

1. **运行时间**：程序包含蒙特卡洛仿真，运行时间可能较长（约1-3分钟），请耐心等待
2. **内存占用**：某些函数（如循环谱分析）计算量较大，建议在配置较好的计算机上运行
3. **参数选择**：
   - 采样频率 `fs` 应至少为信号带宽的2倍（满足奈奎斯特定理）
   - 信噪比范围应根据实际应用场景设置
   - 子载波数目 `para` 建议使用2的幂次（64, 128, 256, 512）
4. **图形显示**：设置 `K=1` 时会显示图形，`K=0` 时不显示（适合批量处理）
5. **函数路径**：确保所有 `.m` 文件都在MATLAB路径中

## 🔧 故障排除

### 常见问题

1. **函数未找到错误**
   - 确保所有 `.m` 文件都在MATLAB路径中
   - 运行前执行：`addpath(genpath('OFDM_Signal_Detection'))`

2. **工具箱缺失错误**
   - 确保已安装 Signal Processing Toolbox
   - 检查：`ver` 命令查看已安装的工具箱

3. **内存不足**
   - 减少SNR数组的长度
   - 减少蒙特卡洛仿真的重复次数

4. **图形不显示**
   - 检查MATLAB图形显示设置
   - 确保 `K=1`（图形显示标志）

## 👤 作者信息

**作者**：Songzhiyong  
**创建日期**：2017.03.01  
**版本**：v1.0

## 📄 许可证

本项目采用 MIT 许可证。

---

**提示**：每个函数文件都包含详细的参数注释和功能说明，如有问题请查看函数文件中的注释或参考 `使用说明.md` 文档。
