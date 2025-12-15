# OFDM项目文件关系说明

本文档详细说明了 `test3.m`、`signalGenerate.m`、`method.m`、`estimateBandwidth.m` 四个文件之间的关系和使用方法。

---

## 📋 文件概览

| 文件名 | 类型 | 主要功能 | 状态 |
|--------|------|----------|------|
| `test3.m` | 完整仿真脚本 | OFDM系统完整链路仿真（发送端→信道→接收端） | 主程序 |
| `signalGenerate.m` | 信号生成脚本 | 生成OFDM发送信号（Tx_data） | 独立模块 |
| `method.m` | 信号检测脚本 | OFDM信号识别和子载波检测 | 独立模块 |
| `estimateBandwidth.m` | 函数 | 带宽估计（Welch算法和AR模型法） | 工具函数 |

---

## 🔗 文件关系图

```
┌─────────────────────────────────────────────────────────────┐
│                    test3.m (完整仿真)                        │
│  ┌──────────────────────────────────────────────────────┐   │
│  │  发送端: 信号生成部分 (已拆分到 signalGenerate.m)   │   │
│  │  - 16QAM调制                                        │   │
│  │  - IFFT变换                                         │   │
│  │  - 循环前缀(CP)                                     │   │
│  │  - 升余弦加窗                                       │   │
│  └──────────────────────────────────────────────────────┘   │
│  ┌──────────────────────────────────────────────────────┐   │
│  │  信道: Jakes瑞利衰落 + AWGN                          │   │
│  └──────────────────────────────────────────────────────┘   │
│  ┌──────────────────────────────────────────────────────┐   │
│  │  接收端: 去CP、FFT、信道估计、均衡、解调             │   │
│  │  - BER-SNR曲线计算                                   │   │
│  │  - 可视化分析                                       │   │
│  └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
                            │
                            │ 信号生成部分拆分
                            ▼
┌─────────────────────────────────────────────────────────────┐
│              signalGenerate.m (信号生成模块)                │
│  - 生成 Tx_data (OFDM发送信号)                             │
│  - 可独立运行                                              │
│  - 基于 test3.m 的信号生成部分                             │
└─────────────────────────────────────────────────────────────┘
                            │
                            │ 生成的信号可用于
                            ▼
        ┌───────────────────┴───────────────────┐
        │                                       │
        ▼                                       ▼
┌──────────────────────┐          ┌──────────────────────┐
│ estimateBandwidth.m  │          │     method.m         │
│ (带宽估计函数)        │          │ (信号识别与检测)      │
│                      │          │                      │
│ - 自动调用           │          │ - OFDM信号识别       │
│   signalGenerate.m   │          │ - 子载波检测         │
│ - Welch算法          │          │ - 带宽检测率         │
│ - AR模型法           │          │ - Rayleigh信道分析   │
└──────────────────────┘          └──────────────────────┘
```

---

## 📝 详细说明

### 1. test3.m - 完整OFDM系统仿真

**功能：**
- 实现完整的OFDM系统链路仿真
- 发送端：16QAM调制、导频插入、IFFT、CP、加窗
- 信道：Jakes瑞利衰落 + AWGN加性高斯白噪声
- 接收端：去CP、FFT、信道估计、均衡、16QAM解调
- 计算BER-SNR曲线（10-30dB，步进2dB）
- 生成14-15个Figure进行可视化分析

**特点：**
- 完整的端到端仿真
- 包含详细的性能分析和可视化
- 支持导频辅助信道估计与均衡

**使用场景：**
- 需要完整OFDM系统性能评估
- 需要BER-SNR曲线分析
- 需要接收端处理流程验证

**输出：**
- `Tx_data`: 发送信号
- `Rx_data`: 接收信号
- BER-SNR曲线
- 多个可视化图表

---

### 2. signalGenerate.m - OFDM信号生成模块

**功能：**
- 生成完整的OFDM发送信号（基带复数信号）
- 支持16QAM调制、IFFT、循环前缀、升余弦加窗等完整流程
- 采用LTE风格的符号重叠相加机制

**特点：**
- **独立运行**：不依赖其他脚本
- **模块化设计**：从test3.m拆分出来，专注于信号生成
- **参数一致**：与test3.m使用相同的参数配置

**参数配置：**
```matlab
subcarrier_spacing = 15e3;      % 子载波间隔：15 kHz
fs = 15.36e6;                   % 采样频率：15.36 MHz
carrier_count = 300;            % 有效数据子载波数
IFFT_bin_length = 1024;        % IFFT点数
GI = 72;                        % 循环前缀长度
beta = 1/16;                    % 加窗滚降系数
```

**使用场景：**
- 只需要生成OFDM信号，不需要完整系统仿真
- 作为其他模块的信号源
- 快速测试信号生成功能

**输出：**
- `Tx_data`: OFDM发送信号（行向量）
- `fs`: 采样频率
- `carrier_count`: 子载波数
- `subcarrier_spacing`: 子载波间隔

**使用示例：**
```matlab
% 直接运行生成信号
signalGenerate

% 信号保存在变量 Tx_data 中
% 其他参数保存在工作空间
```

---

### 3. method.m - OFDM信号识别与检测

**功能：**
- OFDM信号识别和子载波检测
- 带宽检测率计算（Rayleigh信道）
- 功率谱密度（PSD）估计
- 支持多种检测算法

**特点：**
- 独立的信号检测模块
- 支持Rayleigh信道分析
- 使用不同的参数配置（与test3.m不同）

**参数配置：**
```matlab
N = 20;                         % OFDM符号数量
fc = 10e6;                      % 载波频率：10 MHz
fs = 40e6;                      % 采样频率：40 MHz
snr = -4:2:10;                  % 信噪比范围
fmax = 20;                      % 最大多普勒频率
```

**使用场景：**
- OFDM信号识别
- 子载波数量估计
- 带宽检测率分析
- Rayleigh信道下的性能评估

**输出：**
- 带宽检测率（Welch算法和AR模型法）
- PSD图形
- 检测性能曲线

**注意：**
- 使用不同的参数配置，与test3.m/signalGenerate.m的参数不同
- 主要用于信号检测和识别任务

---

### 4. estimateBandwidth.m - 带宽估计函数

**功能：**
- 使用Welch算法和AR模型法估计OFDM信号的带宽
- 支持AWGN信道下的带宽估计
- 自动调用signalGenerate.m生成信号（如果工作空间没有信号）

**特点：**
- **智能依赖**：如果工作空间没有Tx_data，会自动调用signalGenerate.m
- **灵活调用**：支持无参数调用和完整参数调用
- **详细输出**：输出估计结果、参数和仿真时间

**函数签名：**
```matlab
[B_welch, B_ar, B_ideal, results] = estimateBandwidth(...
    Tx_data, fs, carrier_count, subcarrier_spacing, ...
    'fc', 0, 'snr', 20, 'plot', true)
```

**使用场景：**
- 快速估计OFDM信号带宽
- 评估带宽估计算法性能
- 对比Welch算法和AR模型法

**使用示例：**

**方法1：无参数调用（自动生成信号）**
```matlab
% 直接运行，如果没有信号会自动调用signalGenerate.m
[B_welch, B_ar, B_ideal, results] = estimateBandwidth;
```

**方法2：使用已生成的信号**
```matlab
% 先运行signalGenerate.m
signalGenerate

% 然后调用estimateBandwidth（会使用已生成的信号）
[B_welch, B_ar, B_ideal, results] = estimateBandwidth;
```

**方法3：提供完整参数**
```matlab
[B_welch, B_ar, B_ideal, results] = estimateBandwidth(...
    Tx_data, fs, carrier_count, subcarrier_spacing, ...
    'fc', 0, 'snr', 20, 'plot', true);
```

**输出：**
- `B_welch`: Welch算法估计的带宽（Hz）
- `B_ar`: AR模型法估计的带宽（Hz）
- `B_ideal`: 理论带宽值（Hz）
- `results`: 结构体，包含详细的估计结果和误差信息

---

## 🔄 文件间依赖关系

### 依赖关系图

```
test3.m
  └─ (独立运行，包含完整功能)

signalGenerate.m
  └─ (独立运行，不依赖其他文件)

estimateBandwidth.m
  ├─ 可选依赖: signalGenerate.m (自动调用)
  └─ 需要: PSD_OFDM.m, Burg.m, Proximate.m 等函数

method.m
  └─ (独立运行，使用不同的参数配置)
```

### 详细依赖说明

1. **test3.m**
   - ✅ 完全独立，不依赖其他文件
   - 包含完整的信号生成、信道传输、接收处理功能

2. **signalGenerate.m**
   - ✅ 完全独立，不依赖其他文件
   - 从test3.m的信号生成部分拆分而来
   - 生成的信号可以被其他模块使用

3. **estimateBandwidth.m**
   - 🔗 **可选依赖**: signalGenerate.m（如果工作空间没有信号会自动调用）
   - 🔗 **需要**: PSD_OFDM.m及相关函数（Burg.m, Proximate.m等）
   - 可以接受test3.m或signalGenerate.m生成的信号

4. **method.m**
   - ✅ 完全独立，不依赖其他文件
   - 使用不同的参数配置，与其他文件不兼容

---

## 📊 参数对比

| 参数 | test3.m / signalGenerate.m | method.m |
|------|----------------------------|----------|
| 子载波间隔 | 15 kHz | - |
| 采样频率 | 15.36 MHz | 40 MHz |
| 载波频率 | 0 (基带) | 10 MHz |
| 子载波数 | 300 | 128 |
| IFFT点数 | 1024 | - |
| CP长度 | 72 | - |
| OFDM符号数 | 100 | 20 |
| SNR范围 | 15 dB (固定) | -4:2:10 dB |

**注意：** `method.m` 使用不同的参数配置，与其他文件生成的信号不兼容。

---

## 🎯 使用建议

### 场景1：完整OFDM系统仿真
```matlab
% 使用 test3.m
test3
% 输出：完整的系统性能分析、BER-SNR曲线、可视化图表
```

### 场景2：仅生成信号
```matlab
% 使用 signalGenerate.m
signalGenerate
% 输出：Tx_data 信号及相关参数
```

### 场景3：快速带宽估计
```matlab
% 使用 estimateBandwidth.m（自动生成信号）
estimateBandwidth
% 或
signalGenerate
estimateBandwidth
```

### 场景4：信号识别与检测
```matlab
% 使用 method.m
method
% 输出：带宽检测率、PSD图形等
```

---

## 📚 文件历史

- **test3.m**: 原始完整仿真程序
- **signalGenerate.m**: 2025.12.10 - 从test3.m拆分出来的信号生成模块
- **estimateBandwidth.m**: 2025.12.10 - 从signalGenerate.m拆分出来的带宽估计函数
- **method.m**: 2017.03.01 - 独立的信号识别和检测模块（使用不同参数）

---

## ⚠️ 注意事项

1. **参数兼容性**
   - `test3.m` 和 `signalGenerate.m` 使用相同的参数，生成的信号兼容
   - `method.m` 使用不同的参数配置，与其他文件不兼容
   - `estimateBandwidth.m` 可以处理 `test3.m` 或 `signalGenerate.m` 生成的信号

2. **工作空间变量**
   - `signalGenerate.m` 和 `test3.m` 会在工作空间创建 `Tx_data`、`fs`、`carrier_count` 等变量
   - `estimateBandwidth.m` 会尝试从工作空间获取这些变量
   - 如果工作空间没有信号，`estimateBandwidth.m` 会自动调用 `signalGenerate.m`

3. **函数依赖**
   - `estimateBandwidth.m` 需要 `PSD_OFDM.m`、`Burg.m`、`Proximate.m` 等函数
   - 确保这些函数在MATLAB路径中

4. **随机性**
   - 所有文件都使用 `rng('shuffle')` 生成随机数据
   - 每次运行会生成不同的信号（除非设置固定随机种子）

---

## 📞 联系信息

如有问题或建议，请参考各文件头部的注释信息。

---

**最后更新：** 2025.12.10

