# OFDM信号检测与参数估计项目

## 📋 项目概述

本项目是一个基于MATLAB的**OFDM（正交频分复用）信号检测与参数估计**仿真系统，主要用于：

- OFDM信号的生成、识别与分类
- 信号参数估计（带宽、子载波数量、循环前缀长度等）
- 不同信道环境下的性能评估（AWGN、Rayleigh衰落信道）
- 多种信号识别算法的实现与对比

本项目适用于通信系统研究、信号处理算法验证、以及毕业设计等场景。

## ✨ 主要功能

### 1. 信号生成模块
- **OFDM信号生成** (`ofdm.m`, `PSD_generate.m`)
  - 支持16QAM调制
  - 可配置子载波数量（默认128）
  - 可配置循环前缀比例（默认1/4）
  - 支持过采样和窗函数平滑过渡

- **其他调制信号生成**
  - QPSK (`QPSK.m`)
  - QAM16 (`QAM16.m`)
  - QAM64 (`QAM64.m`)
  - FSK8 (`FSK8.m`)
  - ASK8 (`ASK8.m`)

### 2. 信号识别模块
- **高阶累积量识别** (`cumulant.m`)
  - 使用四阶累积量（C40, C42）进行信号分类
  - 支持OFDM、QPSK、QAM16、QAM64、FSK8等信号识别

- **小波变换识别** (`sc_ofdm_wavelet.m`, `sc_ofdm_wavelet_Ray.m`)
  - 基于小波变换的信号识别方法
  - 支持AWGN和Rayleigh信道

- **OFDM信号识别率计算** (`OFDM_rate.m`, `OFDM_rate_xiaobo.m`, `OFDM_rate_xiaobo_Ray.m`)
  - 计算不同SNR下的识别准确率

### 3. 参数估计模块

#### 3.1 带宽估计
- **AWGN信道** (`Bandwidth_rate.m`)
  - AR模型法（Burg算法）
  - Welch算法
  - 计算不同SNR下的带宽检测准确率

- **Rayleigh信道** (`Bandwidth_rate_rayleigh.m`)
  - 支持多径衰落和多普勒效应
  - 对比AR模型和Welch算法性能

- **功率谱密度估计** (`PSD_OFDM.m`, `PSD_OFDM_rayleigh.m`)
  - AR模型功率谱估计
  - Welch算法功率谱估计
  - 直接FFT频谱分析

#### 3.2 子载波数量估计
- **AWGN信道** (`carrier_number.m`, `Rate_Carrynum.m`)
- **Rayleigh信道** (`carrier_number_rayleigh.m`, `Rate_Carrynum_rayleigh.m`)
- **基于带宽估计** (`solve_carrynum_rate.m`)

#### 3.3 循环前缀长度估计
- **有效数据长度估计** (`effectivelength.m`, `effectivelength_rayleigh.m`)
- **符号总长度估计**
- **循环前缀长度估计**
- **不同参数影响分析** (`Length_rate.m`, `Length_rate_difpara.m`, `Length_rate_difratio.m`)

### 4. 信道仿真模块
- **Rayleigh多径衰落信道** (`MUL_RAYLEIGH.m`)
  - 支持多径延时配置
  - 支持多普勒频移
  - 可配置各径功率衰减

### 5. 信号分析模块
- **循环谱估计** (`cyclic_spectrum.m`)
  - **载波频率估计**：基于循环谱方法估计OFDM信号的载波频率
    - 使用平滑周期图法计算循环谱
    - 提取α=0切片（功率谱）找到峰值位置
    - 适用于AWGN和Rayleigh信道
    - 详细实现说明请参考 [`载波频率估计说明.md`](载波频率估计说明.md)
  - 支持AWGN和Rayleigh信道

- **高阶累积量分析** (`ofdm_dif_para.m`)
  - 分析子载波数量对累积量估计的影响

- **自相关分析** (`auto_xcorr.m`)
  - 用于符号长度估计

## 🚀 快速开始

### 环境要求

- **MATLAB版本**：R2014b或更高版本
- **必需工具箱**：
  - Signal Processing Toolbox（信号处理工具箱）
  - Communications Toolbox（通信工具箱，用于QAM调制）

### 安装步骤

1. **下载项目**
   ```bash
   # 将所有文件下载到本地目录
   ```

2. **添加路径**
   在MATLAB命令窗口中执行：
   ```matlab
   cd '项目目录路径'
   addpath(genpath('.'))
   ```

### 运行示例

#### 示例1：运行主演示程序（test.m）

```matlab
% 直接运行主程序
test
```

该程序将：
- 生成OFDM测试信号
- 计算Rayleigh信道下的带宽检测率
- 生成功率谱密度估计图
- 绘制性能对比曲线

#### 示例2：运行简化演示程序（OFDM_DEMO.m）

```matlab
% 运行演示程序
OFDM_DEMO
```

#### 示例3：自定义参数运行

```matlab
% 设置参数
N = 20;              % OFDM符号数量
para = 128;          % 子载波数量
ratio = 1/4;         % 循环前缀比例
snr = -4:2:10;       % 信噪比范围
fc = 10e6;           % 载波频率
fs = 40e6;           % 采样频率

% 生成OFDM信号
sig_ofdm = PSD_generate(N);

% 计算带宽检测率（Rayleigh信道）
itau = floor([0,1e-8,2e-8,5e-8,2e-7,5e-7].*fs);
power = [0,-1.0,-7.0,-10.0,-12.0,-17.0];
fmax = 20;
itn = [10000,20000,30000,40000,50000,60000];

[B_rate_welch, B_rate_ar] = Bandwidth_rate_rayleigh(...
    sig_ofdm, fc, fs, snr, itau, power, fmax, itn);

% 绘制结果
figure;
plot(snr, B_rate_ar, 'r-o');
hold on;
plot(snr, B_rate_welch, 'b--s');
xlabel('SNR (dB)');
ylabel('检测准确率 (%)');
legend('AR模型法', 'Welch算法');
title('Rayleigh信道下的带宽检测率');
grid on;
```

## 📁 项目结构

```
OFDMSignalDetection/
├── README.md                    # 项目说明文档（本文件）
├── test.md                      # test.m详细功能说明
│
├── 主程序文件/
│   ├── test.m                   # 主测试程序（完整功能）
│   ├── OFDM_DEMO.m              # 简化演示程序
│
├── 信号生成模块/
│   ├── ofdm.m                   # OFDM信号生成
│   ├── PSD_generate.m           # 用于PSD测试的OFDM信号生成
│   ├── QPSK.m                   # QPSK信号生成
│   ├── QAM16.m                  # 16QAM信号生成
│   ├── QAM64.m                  # 64QAM信号生成
│   ├── FSK8.m                   # 8FSK信号生成
│   ├── ASK8.m                   # 8ASK信号生成
│   └── random_binary.m          # 随机二进制序列生成
│
├── 信号识别模块/
│   ├── cumulant.m               # 高阶累积量信号识别
│   ├── OFDM_rate.m              # OFDM信号识别率（AWGN）
│   ├── OFDM_rate_xiaobo.m       # OFDM信号识别率（小波方法，AWGN）
│   ├── OFDM_rate_xiaobo_Ray.m   # OFDM信号识别率（小波方法，Rayleigh）
│   ├── sc_ofdm_wavelet.m        # 小波变换识别（AWGN）
│   └── sc_ofdm_wavelet_Ray.m    # 小波变换识别（Rayleigh）
│
├── 带宽估计模块/
│   ├── Bandwidth_rate.m         # 带宽检测率（AWGN）
│   ├── Bandwidth_rate_rayleigh.m # 带宽检测率（Rayleigh）
│   ├── BandwidthEstimate.m      # 带宽估计通用函数
│   ├── PSD_OFDM.m               # 功率谱密度估计（AWGN）
│   ├── PSD_OFDM_rayleigh.m      # 功率谱密度估计（Rayleigh）
│   ├── Burg.m                   # Burg算法（AR模型）
│   └── computeARpara.m          # AR参数计算
│
├── 子载波估计模块/
│   ├── carrier_number.m         # 子载波数量估计（AWGN）
│   ├── carrier_number_rayleigh.m # 子载波数量估计（Rayleigh）
│   ├── Rate_Carrynum.m          # 子载波数量估计准确率（AWGN）
│   ├── Rate_Carrynum_rayleigh.m # 子载波数量估计准确率（Rayleigh）
│   └── solve_carrynum_rate.m    # 基于带宽的子载波数量估计
│
├── 循环前缀估计模块/
│   ├── effectivelength.m        # 有效长度估计（AWGN）
│   ├── effectivelength_rayleigh.m # 有效长度估计（Rayleigh）
│   ├── Length_rate.m            # 长度估计准确率（AWGN）
│   ├── Length_rate_rayleigh.m   # 长度估计准确率（Rayleigh）
│   ├── Length_rate_difpara.m    # 不同子载波数的影响
│   └── Length_rate_difratio.m   # 不同循环前缀比例的影响
│
├── 信道仿真模块/
│   └── MUL_RAYLEIGH.m           # 多径Rayleigh衰落信道
│
├── 信号分析模块/
│   ├── cyclic_spectrum.m        # 循环谱估计
│   ├── auto_xcorr.m             # 自相关分析
│   ├── ofdm_dif_para.m          # 不同参数对累积量的影响
│   ├── rt_C21.m                 # 累积量C21计算
│   ├── rt_C40.m                 # 累积量C40计算
│   ├── rt_C42.m                 # 累积量C42计算
│   └── rt_C421.m                # 累积量C421计算
│
└── 辅助函数/
    ├── Proximate.m              # 近似计算
    └── overfind_num.m           # 数值查找辅助函数
```

## ⚙️ 系统参数配置

### 默认参数设置

| 参数 | 默认值 | 说明 |
|------|--------|------|
| **N** | 20 | OFDM符号数量 |
| **para** | 128 | 子载波数量（需在脚本中设置） |
| **ratio** | 1/4 | 循环前缀比例（需在脚本中设置） |
| **M** | 16 | 调制阶数（16QAM） |
| **fc** | 10 MHz | 载波频率 |
| **fs** | 40 MHz | 采样频率 |
| **snr** | -4:2:10 dB | 信噪比范围 |

### Rayleigh信道参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| **多径数量** | 6径 | 多径信道路径数 |
| **多径延时** | [0, 10, 20, 50, 200, 500] ns | 各径延时 |
| **多径功率** | [0, -1, -7, -10, -12, -17] dB | 各径功率衰减 |
| **fmax** | 20 Hz | 最大多普勒频率 |

## 📊 输出结果说明

### 1. 控制台输出
程序运行时会显示：
- 仿真参数信息
- 各阶段执行进度
- 计算耗时统计
- 最终结果摘要

### 2. 图形输出
- **输入信号时域/频域图**：显示上变频并加噪声后的信号特征
- **功率谱密度估计图**：AR模型、Welch算法、FFT三种方法的PSD对比
- **性能对比曲线**：不同SNR下的检测准确率对比

### 3. 工作空间变量
- `sig_ofdm`：生成的OFDM信号
- `B_rate_welch_rayleigh`：Welch算法带宽检测准确率
- `B_rate_ar_rayleigh`：AR模型带宽检测准确率
- 其他中间结果变量

## 🔧 主要函数说明

### 信号生成函数

#### `ofdm(N, para, ratio)`
生成OFDM信号
- **输入**：
  - `N`：符号个数
  - `para`：子载波数量
  - `ratio`：循环前缀比例
- **输出**：归一化的OFDM时域信号

#### `PSD_generate(N)`
生成用于功率谱密度测试的OFDM信号（带窗函数平滑）
- **输入**：`N`：符号个数
- **输出**：处理后的OFDM信号

### 带宽估计函数

#### `Bandwidth_rate_rayleigh(sig, fc, fs, snr, itau, power, fmax, itn)`
计算Rayleigh信道下的带宽检测准确率
- **输入**：
  - `sig`：输入信号
  - `fc`：载波频率
  - `fs`：采样频率
  - `snr`：信噪比数组
  - `itau`：多径延时数组
  - `power`：多径功率数组
  - `fmax`：最大多普勒频率
  - `itn`：Rayleigh信道记录次数
- **输出**：
  - `B_rate_welch`：Welch算法检测准确率
  - `B_rate_ar`：AR模型检测准确率

#### `PSD_OFDM_rayleigh(sig, fc, fs, snr, k, itau, power, fmax, itn)`
生成Rayleigh信道下的功率谱密度估计图
- **输入**：同上，`k`为绘图标志
- **输出**：
  - `B_welch`：Welch算法估计的带宽
  - `B_ar`：AR模型估计的带宽

### 子载波估计函数

#### `carrier_number_rayleigh(yc, N, snr, K, itau, power, fmax, fs, itn)`
估计OFDM信号的子载波数量（Rayleigh信道）
- **输入**：
  - `yc`：接收信号
  - `N`：符号个数
  - `snr`：信噪比
  - `K`：绘图标志
  - 其他Rayleigh信道参数
- **输出**：估计的子载波数量

## 📖 使用示例

### 示例1：OFDM信号识别

```matlab
% 设置参数
N = 20;
para = 128;
ratio = 1/4;
snr = -5:5:25;
K = 1;  % 绘图标志

% 计算不同信号的高阶累积量
[c42_ofdm, c42_qpsk, c42_qam16, c42_qam64, c42_fsk8] = ...
    cumulant(snr, N, para, ratio, K);

% 计算OFDM信号识别率
rate_ofdm = OFDM_rate(snr, N, para, ratio);
```

### 示例2：带宽估计性能对比

```matlab
% 生成信号
sig_ofdm = PSD_generate(20);

% AWGN信道
[B_rate_welch_awgn, B_rate_ar_awgn] = ...
    Bandwidth_rate(sig_ofdm, 10e6, 40e6, -4:2:10);

% Rayleigh信道
itau = floor([0,1e-8,2e-8,5e-8,2e-7,5e-7].*40e6);
power = [0,-1.0,-7.0,-10.0,-12.0,-17.0];
fmax = 20;
itn = [10000,20000,30000,40000,50000,60000];

[B_rate_welch_ray, B_rate_ar_ray] = ...
    Bandwidth_rate_rayleigh(sig_ofdm, 10e6, 40e6, -4:2:10, ...
    itau, power, fmax, itn);

% 对比绘图
snr = -4:2:10;
figure;
subplot(1,2,1);
plot(snr, B_rate_ar_awgn, 'r-o', snr, B_rate_welch_awgn, 'b--s');
legend('AR (AWGN)', 'Welch (AWGN)');
title('AWGN信道');
xlabel('SNR (dB)'); ylabel('准确率 (%)');

subplot(1,2,2);
plot(snr, B_rate_ar_ray, 'r-o', snr, B_rate_welch_ray, 'b--s');
legend('AR (Rayleigh)', 'Welch (Rayleigh)');
title('Rayleigh信道');
xlabel('SNR (dB)'); ylabel('准确率 (%)');
```

### 示例3：循环前缀长度估计

```matlab
% 生成信号
sig = ofdm(20, 128, 1/4);
sig_up = sig.*exp(1i*2*pi*10e6/40e6*(0:length(sig)-1));

% AWGN信道估计
[Tu, Ts, Tg] = effectivelength(sig_up, 40e6, 20, 20, 1);

% 不同SNR下的准确率
snr = -4:2:10;
[Tu_rate, Ts_rate, Tg_rate] = Length_rate(sig_up, 40e6, snr, 20, 0);

% 绘图
figure;
plot(snr, Tu_rate, 'r-o', 'DisplayName', '有效数据长度');
hold on;
plot(snr, Ts_rate, 'b--s', 'DisplayName', '符号总长度');
plot(snr, Tg_rate, 'g-*', 'DisplayName', '循环前缀长度');
legend;
xlabel('SNR (dB)');
ylabel('估计准确率 (%)');
title('循环前缀长度估计性能');
```

## ⚠️ 注意事项

1. **参数设置**：部分参数（如`para`、`ratio`）在脚本中默认被注释，需要使用时请取消注释或手动赋值

2. **运行时间**：包含蒙特卡洛仿真的程序运行时间较长（约1-2分钟），请耐心等待

3. **内存占用**：处理较长信号时可能占用较多内存

4. **函数依赖**：确保所有依赖函数都在MATLAB路径中

5. **工具箱要求**：需要Signal Processing Toolbox和Communications Toolbox

6. **版本兼容性**：代码中使用了`qammod`函数，需要MATLAB R2014b或更高版本

## 🔍 技术细节

### 带宽估计方法

1. **AR模型法（Burg算法）**
   - 使用AIC准则自动选择模型阶数
   - 根据SNR选择不同的阈值点（-6dB、-5dB或-3dB）
   - 计算功率谱密度并提取带宽

2. **Welch算法**
   - 使用汉宁窗，重叠55个采样点，每段100个采样点
   - 计算功率谱密度
   - 使用-3dB点提取带宽

### 信道模型

1. **AWGN信道**：加性高斯白噪声
2. **Rayleigh信道**：
   - 多径传播环境
   - 6条多径，不同延时和功率
   - 考虑多普勒效应

### 信号生成流程

1. 生成随机16QAM调制数据
2. IFFT变换将频域数据转换为时域OFDM符号
3. 添加循环前缀（比例1/4）
4. 补零处理（用于测试）
5. 功率归一化

## 📝 版本信息

- **版本**：v1.0
- **创建日期**：2017.03.01
- **作者**：Songzhiyong
- **最后更新**：2025

## 📚 相关文档

- `test.md`：test.m程序的详细功能说明文档
- `载波频率估计说明.md`：载波频率估计算法的详细实现说明
- 各函数文件中的注释：包含详细的参数说明和功能描述

## 🤝 贡献

欢迎提出问题和改进建议！

## 📄 许可证

本项目仅供学习和研究使用。

---

**提示**：如有问题，请查看各函数文件中的注释说明，每个函数都包含详细的参数说明和功能描述。

