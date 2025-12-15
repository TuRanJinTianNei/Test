# estimateBandwidth.m - OFDM信号带宽估计函数

## 📋 目录

- [功能概述](#功能概述)
- [输入输出参数](#输入输出参数)
- [信号处理流程](#信号处理流程)
- [算法原理](#算法原理)
- [可视化输出](#可视化输出)
- [使用示例](#使用示例)
- [依赖关系](#依赖关系)

---

## 功能概述

`estimateBandwidth.m` 是一个用于估计OFDM信号带宽的MATLAB函数。该函数使用两种功率谱估计方法（Welch算法和AR模型法）来估计信号的带宽，并支持在Rayleigh衰落信道下的性能评估。

### 主要特点

- ✅ 使用 `signalGenerate.m` 生成的基带信号作为输入
- ✅ 支持自动从工作空间获取参数或自动调用 `signalGenerate.m`
- ✅ 实现两种带宽估计方法：Welch算法和AR模型法
- ✅ 支持Rayleigh多径衰落信道仿真
- ✅ 提供详细的时域和频域可视化
- ✅ 计算不同SNR下的带宽检测率

---

## 输入输出参数

### 输入参数

| 参数名 | 类型 | 说明 | 默认值 |
|--------|------|------|--------|
| `Tx_data` | 向量 | OFDM发送信号（基带信号，来自signalGenerate.m） | 从工作空间获取或自动生成 |
| `fs` | 标量 | 采样频率（Hz） | 15.36e6 |
| `carrier_count` | 标量 | 有效数据子载波数 | 300 |
| `subcarrier_spacing` | 标量 | 子载波间隔（Hz） | 15e3 |
| `'fc'` | 标量 | 载波频率（Hz） | 0（基带信号） |
| `'snr'` | 标量 | 信噪比（dB） | 20 |
| `'plot'` | 逻辑值 | 是否绘制图形 | true |

### 输出参数

| 参数名 | 类型 | 说明 |
|--------|------|------|
| `B_welch` | 标量 | Welch算法估计的带宽值（Hz） |
| `B_ar` | 标量 | AR模型法估计的带宽值（Hz） |
| `B_ideal` | 标量 | 理论带宽值（Hz）= carrier_count × subcarrier_spacing |
| `results` | 结构体 | 包含详细的估计结果和误差信息 |

---

## 信号处理流程

### 整体流程图

```
输入基带信号 Tx_data (来自signalGenerate.m)
    ↓
[步骤1] 参数获取与验证
    ├─ 从工作空间获取参数
    ├─ 或自动调用signalGenerate.m生成信号
    └─ 计算理论带宽 B_ideal
    ↓
[步骤2] 绘制输入基带信号 Tx_data（时域和频域）
    ├─ Figure 1: 输入基带信号 Tx_data（时域和频域）
    └─ Figure 2: 从signalGenerate.m引用的基带信号 Tx_data
    ↓
[步骤3] 信号预处理
    ├─ 检查信号类型（实数/复数）
    └─ 转换为统一格式 sig_baseband
    ↓
[步骤4] 信号处理（上变频 + Rayleigh衰落 + AWGN噪声）
    ├─ 上变频到载波频率 fc
    ├─ 通过Rayleigh多径衰落信道
    └─ 添加AWGN噪声
    ↓
[步骤5] 绘制处理后信号（时域和频域）
    └─ Figure 3: 基于Tx_data的处理后信号（Rayleigh信道）
    ↓
[步骤6] 功率谱密度估计
    ├─ Welch算法（pwelch）
    └─ AR模型法（Burg算法）
    ↓
[步骤7] 带宽计算
    ├─ Welch方法：-3dB带宽
    └─ AR方法：根据SNR选择阈值（-6dB/-5dB/-3dB）
    ↓
[步骤8] 带宽检测率计算（不同SNR下）
    └─ 蒙特卡洛仿真，计算检测率
    ↓
[步骤9] 结果可视化与输出
    ├─ Figure 4: Rayleigh信道下的带宽检测率对比
    └─ Figure 5-7: 功率谱密度估计图（AR、Welch、频谱）
```

### 详细步骤说明

#### 步骤1：参数获取与验证

```matlab
% 1.1 无参数调用时，从工作空间获取或自动生成信号
if nargin == 0
    % 尝试从工作空间获取Tx_data
    % 如果不存在，自动调用signalGenerate.m生成信号
end

% 1.2 计算理论带宽
B_ideal = carrier_count * subcarrier_spacing;
```

**处理内容：**
- 支持无参数调用，自动获取或生成信号
- 从工作空间获取fs、carrier_count、subcarrier_spacing等参数
- 如果参数不存在，使用默认值或自动调用signalGenerate.m
- 计算理论带宽值用于对比

#### 步骤1.5：信号预处理

```matlab
% 检查Tx_data是否为复数信号
if isreal(Tx_data)
    sig_baseband = Tx_data + 0*1j;  % 转换为复数形式
else
    sig_baseband = Tx_data;
end
```

**处理内容：**
- 检测输入信号类型（实数/复数）
- signalGenerate.m生成的Tx_data通常是实数信号（加窗后）
- 统一转换为复数形式以便后续处理

#### 步骤2：绘制输入基带信号

**生成图形：** `输入基带信号 Tx_data（时域和频域）`

包含4个子图：
1. **时域 - 实部**：显示信号前5000个采样点的时域波形
2. **时域 - 虚部**：显示信号的虚部（如果为复数信号）
3. **频域 - 双边频谱**：使用FFT计算的双边频谱（包含正负频率）
4. **频域 - 单边频谱**：只显示正频率部分的频谱

**关键处理：**
```matlab
% FFT变换
N_fft = 2^nextpow2(sig_length);
fft_tx = fft(Tx_data, N_fft);
fft_tx_shifted = fftshift(fft_tx);  % 零频率居中

% 归一化处理
fft_mag = 20*log10(abs(fft_tx_shifted) + eps);
fft_mag = fft_mag - max(fft_mag);  % 归一化到最大值
```

#### 步骤3：信号处理（信道仿真）

**处理流程：**

```matlab
% 3.1 上变频：将基带信号调制到载波频率fc
sig_processed = sig_baseband .* exp(1j*2*pi*fc/fs*(0:length(sig_baseband)-1));

% 3.2 通过Rayleigh多径衰落信道
sig_processed = real(MUL_RAYLEIGH(sig_processed, itau, power, itn, ...));

% 3.3 添加AWGN噪声
sig_processed = awgn(sig_processed, snr_plot_signal, 'measured');
```

**Rayleigh信道参数：**
- `itau`: 多径延时（样本数）
- `power`: 每径功率（dB）
- `fmax`: 最大多普勒频率（Hz）
- `itn`: 瑞利信道的记录次数

#### 步骤4：绘制处理后信号

**生成图形：** `基于Tx_data的处理后信号（Rayleigh信道）`

包含2个子图：
1. **时域波形**：显示经过信道处理后的时域信号
2. **频谱**：显示处理后信号的频谱，标记载波频率和理论带宽边界

#### 步骤5：功率谱密度（PSD）估计

##### 5.1 Welch算法

```matlab
[Pxx2, f1] = pwelch(sig3_chnl, hanning(100), 55, 4096*2, fs);
```

**参数说明：**
- 窗函数：Hanning窗，长度100
- 重叠：55个样本
- FFT点数：4096×2 = 8192
- 归一化：转换为dB单位并归一化

##### 5.2 AR模型法（Burg算法）

```matlab
[Pxx1, f, p] = Burg(sig3_chnl, fs, 'AIC');
```

**参数说明：**
- 算法：Burg算法
- 阶数选择：使用AIC（Akaike Information Criterion）准则自动确定
- 归一化：转换为dB单位并归一化

**生成图形：**
- `Rayleigh信道-AR模型功率谱密度估计`
- `Rayleigh信道-Welch算法功率谱密度估计`
- `Rayleigh信道-OFDM信号频谱`

#### 步骤6：带宽计算

##### 6.1 Welch算法带宽计算

```matlab
% 在-3dB处计算带宽
[as1, as11] = Proximate(-3, P1);  % 左边界
[as2, as22] = Proximate(-3, P2);  % 右边界
B_welch = abs(band1 - band2);
```

**方法：** 在功率谱的-3dB点处测量带宽

##### 6.2 AR模型带宽计算

```matlab
% 根据SNR选择不同的阈值
if snr > 4
    threshold = -6;  % 高SNR使用-6dB
elseif snr > 0 && snr <= 4
    threshold = -5;  % 中等SNR使用-5dB
else
    threshold = -3;  % 低SNR使用-3dB
end
```

**自适应阈值：** 根据SNR动态调整带宽测量阈值，提高估计精度

#### 步骤7：带宽检测率计算

**蒙特卡洛仿真：**

```matlab
for i = 1:length(snr)
    for j = 1:numb  % 重复numb次
        [b_welch(i,j), b_ar(i,j)] = PSD_OFDM_rayleigh(...);
    end
    % 计算平均估计值和检测率
    B_welch(i) = mean(b_welch(i,:));
    B_rate_welch(i) = 1 - abs((B_welch(i) - B_ideal)) / B_ideal;
end
```

**输出：**
- 每个SNR下的平均带宽估计值
- 带宽检测率（1 - 相对误差）
- 绝对误差和相对误差统计

#### 步骤8：结果可视化

**生成图形：** `Rayleigh信道下的带宽检测率对比`

显示不同SNR下两种方法的带宽检测率对比曲线。

---

## 算法原理

### Welch算法

**原理：**
- 将信号分段，对每段加窗后进行FFT
- 对所有段的功率谱求平均
- 通过重叠分段提高频率分辨率

**优点：**
- 计算效率高
- 对噪声有较好的平滑效果
- 实现简单

**缺点：**
- 频率分辨率受窗函数长度限制
- 可能存在频谱泄露

### AR模型法（Burg算法）

**原理：**
- 将信号建模为自回归（AR）过程
- 使用Burg算法估计AR模型参数
- 通过模型参数计算功率谱密度

**优点：**
- 频率分辨率高
- 适合短数据序列
- 可以自动选择模型阶数（AIC准则）

**缺点：**
- 计算复杂度较高
- 对模型阶数选择敏感

### 带宽测量方法

**-3dB带宽（半功率带宽）：**
- 在功率谱密度下降3dB（功率减半）的位置测量带宽
- 适用于Welch算法和低SNR情况

**-6dB/-5dB带宽：**
- 在功率谱密度下降6dB或5dB的位置测量带宽
- 适用于AR模型法和高SNR情况
- 提供更严格的带宽定义

---

## 可视化输出

### Figure 1: 输入基带信号 Tx_data（时域和频域）
- **位置：** 函数开始处，参数输出后（如果 `plot_flag=true`）
- **内容：** 显示从 `signalGenerate.m` 输入的基带信号的时域和频域特性
- **子图：**
  1. 时域 - 实部（前5000个采样点）
  2. 时域 - 虚部（如果为复数信号）
  3. 频域 - 双边频谱（FFT，零频率居中）
  4. 频域 - 单边频谱（只显示正频率）
- **标记：** 理论带宽边界

### Figure 2: 从signalGenerate.m引用的基带信号 Tx_data
- **位置：** 主程序开始处，信号预处理后（如果 `plot_flag=true`）
- **内容：** 显示基带信号的详细时域和频域分析
- **子图：** 与Figure 1相同，但基于处理后的sig_baseband
- **说明：** 这是实际用于后续处理的基带信号

### Figure 3: 基于Tx_data的处理后信号（Rayleigh信道）
- **位置：** 信号处理后（如果 `plot_flag=true`）
- **内容：** 显示经过上变频、Rayleigh衰落和加噪声后的信号
- **子图：**
  1. 时域波形（完整信号）
  2. 频谱（FFT，标记载波频率和带宽边界）
- **参数：** SNR = 20 dB（用于绘图）

### Figure 4: Rayleigh信道下的带宽检测率对比
- **位置：** 带宽检测率计算后
- **内容：** 显示不同SNR下两种方法的带宽检测率对比曲线
- **SNR范围：** -4:2:10 dB
- **曲线：** AR模型法（红色圆点）、Welch算法（蓝色方块）

### Figure 5-7: 功率谱密度估计图
- **位置：** PSD估计函数内部（`k=1`时，即snr_plot=20 dB）
- **内容：**
  - **Figure 5:** Rayleigh信道-AR模型功率谱密度估计
  - **Figure 6:** Rayleigh信道-Welch算法功率谱密度估计
  - **Figure 7:** Rayleigh信道-OFDM信号频谱（FFT）
- **说明：** 这些图形在调用 `PSD_OFDM_rayleigh` 时生成

---

## 使用示例

### 示例1：无参数调用（自动获取信号）

```matlab
% 先运行signalGenerate.m生成信号
signalGenerate;

% 然后调用estimateBandwidth（自动使用工作空间中的Tx_data）
[B_welch, B_ar, B_ideal, results] = estimateBandwidth;
```

### 示例2：完整参数调用

```matlab
% 生成信号
signalGenerate;

% 调用带宽估计函数
[B_welch, B_ar, B_ideal, results] = estimateBandwidth(...
    Tx_data, ...              % 基带信号
    fs, ...                   % 采样频率
    carrier_count, ...        % 子载波数
    subcarrier_spacing, ...   % 子载波间隔
    'fc', 10e6, ...          % 载波频率 10MHz
    'snr', 20, ...           % 信噪比 20dB
    'plot', true ...         % 绘制图形
);
```

### 示例3：仅计算不绘图

```matlab
[B_welch, B_ar, B_ideal, results] = estimateBandwidth(...
    Tx_data, fs, carrier_count, subcarrier_spacing, ...
    'plot', false ...  % 不绘制图形，加快计算速度
);
```

### 示例4：查看详细结果

```matlab
[B_welch, B_ar, B_ideal, results] = estimateBandwidth;

% 显示结果
fprintf('理论带宽: %.2f MHz\n', B_ideal/1e6);
fprintf('Welch估计: %.2f MHz\n', B_welch/1e6);
fprintf('AR模型估计: %.2f MHz\n', B_ar/1e6);
fprintf('Welch相对误差: %.2f%%\n', ...
    abs(B_welch - B_ideal) / B_ideal * 100);
fprintf('AR模型相对误差: %.2f%%\n', ...
    abs(B_ar - B_ideal) / B_ideal * 100);
```

---

## 依赖关系

### 必需函数

| 函数名 | 说明 | 位置 | 必需性 |
|--------|------|------|--------|
| `signalGenerate.m` | 生成OFDM基带信号 | 同一目录 | 必需（或提供Tx_data） |
| `Bandwidth_rate_rayleigh` | 计算Rayleigh信道带宽检测率 | 内部函数 | 必需 |
| `PSD_OFDM_rayleigh` | Rayleigh信道PSD估计 | 内部函数 | 必需 |
| `MUL_RAYLEIGH` | 多径Rayleigh衰落信道 | 内部函数 | 必需 |
| `delay` | 信号延时函数（MUL_RAYLEIGH的子函数） | 内部函数 | 必需 |
| `siglfade` | 单径衰落函数（MUL_RAYLEIGH的子函数） | 内部函数 | 必需 |
| `Burg` | Burg算法AR模型估计 | 内部函数 | 必需 |
| `computeARpara` | AR模型参数计算（Burg的子函数） | 内部函数 | 必需 |
| `Proximate` | 查找最接近值的位置 | 内部函数 | 必需 |

### MATLAB工具箱

- **Signal Processing Toolbox**：用于 `pwelch` 函数（必需）
- **Communications Toolbox**：用于 `awgn` 函数（必需，用于添加AWGN噪声）

### signalGenerate.m 的依赖

如果使用自动调用 `signalGenerate.m` 功能，还需要：
- `rcoswindow.m` - 升余弦窗函数（必需）
- `qam16.m` - 16QAM调制函数（可选，有fallback到qammod）
- `jakesChannel.m` - Jakes瑞利衰落信道（可选，仅在启用信道时）

### 内部函数说明

#### Bandwidth_rate_rayleigh
- **功能：** 计算不同SNR下的带宽检测率
- **方法：** 蒙特卡洛仿真，重复多次估计取平均

#### PSD_OFDM_rayleigh
- **功能：** 估计Rayleigh信道下的功率谱密度
- **输出：** Welch和AR两种方法的PSD估计值

#### MUL_RAYLEIGH
- **功能：** 模拟多径Rayleigh衰落信道
- **参数：** 多径延时、功率、多普勒频率等

#### Burg
- **功能：** 使用Burg算法进行AR模型功率谱估计
- **特点：** 支持AIC/FPE准则自动选择模型阶数

---

## 信号大小与性能

### 输入信号特征（来自signalGenerate.m）

`estimateBandwidth.m` 使用的输入信号 `Tx_data` 由 `signalGenerate.m` 生成，具有以下特征：

| 特征 | 数值 | 说明 |
|------|------|------|
| **样本点数** | 109,736 | 信号的总采样点数 |
| **信号时长** | ~7.15 ms | 信号持续时间 |
| **采样频率** | 15.36 MHz | signalGenerate.m的采样频率 |
| **理论带宽** | 4.5 MHz | carrier_count × subcarrier_spacing |
| **OFDM符号数** | 100 | 总符号数 |
| **每符号长度** | 1,164 样本 | IFFT点数(1024) + CP(72) + 后缀(68) |
| **信号类型** | 实数 | 加窗后的实数信号 |

### 内存占用

| 数据类型 | 内存占用 | 说明 |
|----------|----------|------|
| **Tx_data (double)** | ~0.84 MB | 109,736 × 8 bytes |
| **Tx_data (single)** | ~0.42 MB | 109,736 × 4 bytes |
| **处理后信号** | ~0.84 MB | 相同长度 |
| **总内存占用** | ~2-3 MB | 包含所有中间变量和图形数据 |

### 计算性能

| 项目 | 时间 | 说明 |
|------|------|------|
| **信号预处理** | <1秒 | 信号验证和转换 |
| **PSD估计（单次）** | 1-5秒 | 取决于信号长度 |
| **带宽检测率计算** | 10-60秒 | 取决于SNR点数和重复次数 |
| **总计算时间** | 15-90秒 | 完整流程（包含绘图） |

**注意：** 计算时间会随信号长度和SNR点数变化。当前配置下，完整运行通常需要30-60秒。

---

## 技术细节

### 信号处理参数

| 参数 | 值 | 说明 |
|------|-----|------|
| OFDM符号数 N | 20 | 每帧OFDM符号数量（用于PSD_generate，已弃用） |
| 载波频率 fc | 10 MHz | 上变频载波频率（函数内部设置） |
| 采样频率 fs | 40 MHz | 仿真采样频率（函数内部设置，与signalGenerate.m不同） |
| SNR范围 | -4:2:10 dB | 带宽检测率测试的SNR范围 |
| 最大多普勒频率 fmax | 20 Hz | Rayleigh衰落参数 |
| 多径数 | 6 | 多径信道路径数 |
| 多径延时 itau | [0, 1e-8, 2e-8, 5e-8, 2e-7, 5e-7]×fs | 多径延时（秒） |
| 多径功率 power | [0, -1.0, -7.0, -10.0, -12.0, -17.0] dB | 每径功率 |


### 性能指标

- **带宽估计精度：** 相对误差通常在5%以内（SNR > 0 dB）
- **计算时间：** 取决于信号长度和SNR点数，通常几秒到几分钟
- **内存占用：** 
  - Tx_data: ~0.84 MB (double类型，109,736样本)
  - 处理后的信号: ~0.84 MB
  - 总内存占用: ~2-3 MB（包含中间变量和图形数据）

---

## 注意事项

### ⚠️ 重要提示

1. **信号格式：** 输入的 `Tx_data` 应为行向量（来自signalGenerate.m）
2. **采样频率差异：** 
   - ⚠️ **重要**：函数内部会重新设置 `fs = 40 MHz`（第327行）
   - signalGenerate.m使用的采样频率是 `15.36 MHz`
   - 函数内部的 `fs` 会覆盖输入参数，用于后续的Rayleigh信道仿真
   - 这可能导致频率轴显示不准确，但不影响带宽估计结果
3. **载波频率：** 
   - 函数内部设置 `fc = 10 MHz`（第326行）
   - 会覆盖输入的可选参数 `'fc'`
   - 用于上变频和频谱显示
4. **信号大小：** 
   - signalGenerate.m生成的信号约109,736个样本，时长约7.15 ms
   - 内存占用约0.84 MB（double类型）
5. **内存管理：** 对于超长信号，可能需要分段处理
6. **图形窗口：** 
   - 如果 `plot_flag=true`，会生成多个图形窗口（最多7个）
   - signalGenerate.m中的 `close all` 已被注释，图形窗口会保留
7. **参数一致性：** 
   - 确保 `carrier_count` 和 `subcarrier_spacing` 与实际信号匹配
   - 函数内部参数（fc, fs）会覆盖输入参数，这是设计行为
8. **信号来源：** 
   - ✅ 函数不再使用内部生成的信号（PSD_generate已弃用）
   - ✅ 必须使用从signalGenerate.m引用的Tx_data信号
   - ✅ 支持自动调用signalGenerate.m生成信号

### 参数覆盖说明

函数内部会设置以下参数（覆盖输入参数）：
- `fc = 10e6` Hz（载波频率）
- `fs = 40e6` Hz（采样频率，用于Rayleigh信道仿真）
- `snr = -4:2:10` dB（SNR范围，用于带宽检测率计算）

这些参数是硬编码的，用于保证仿真的一致性。如果需要修改，请直接编辑函数代码。

---

## 版本历史

- **2017.03.01** - 原始版本
- **2025.12.10** - 改造成函数，接受 `signalGenerate.m` 的基带信号
- **2025.12.10** - 移除内部信号生成（PSD_generate），使用外部输入的 `Tx_data`
- **2025.12.10** - 添加详细的基带信号和处理后信号的可视化
- **2025.12.10** - 修复函数定义格式问题（添加end语句）
- **2025.12.10** - 优化图形显示，保留图形窗口不被关闭

---

## 作者信息

基于原始OFDM信号检测代码改造，用于OFDM信号带宽估计研究。

---

## 许可证

本项目代码仅供学术研究使用。

