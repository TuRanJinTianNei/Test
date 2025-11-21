### test2.m 使用说明（README）

`test2.m` 是一个基于 LTE 思路实现的 OFDM 发射-信道-接收端完整仿真脚本，包含：
- 16QAM 比特映射与子载波映射（采用频域共轭对称，保证 IFFT 输出为实信号）
- IFFT 生成时域符号、添加循环前缀（区分长/短 CP）与升余弦加窗
- 帧/串行信号构造（在 CP 范围内进行重叠相加，降低旁瓣）
- AWGN 信道、接收端去 CP → FFT → 解调判决、误码率计算
- 多幅图形用于观测频域/时域/频谱、星座图、发收比特等

---

### 环境要求
- MATLAB R2016a 或更高版本（建议更高版本）
- 同目录下依赖函数：
  - `qam16.m`（16QAM 映射）
  - `demoduqam16.m`（16QAM 判决解调）
  - `rcoswindow.m`（升余弦窗）

---

### 运行方法
1. 打开 MATLAB，将当前工作目录切换到本脚本所在路径。
2. 直接运行：
   ```matlab
   test2
   ```

---

### 关键可调参数（在脚本顶部附近）
- `IFFT_bin_length`：IFFT 点数（LTE 10MHz 配置为 2048）
- 子载波配置：正频率数据子载波 `k=1:252`，负频率为其共轭，`carrier_count=252`
- CP 配置：
  - 长 CP 采样点 `CP_long_samples = 160`
  - 短 CP 采样点 `CP_short_samples = 144`
  - `CP_lengths` 根据符号顺序分配（前 2 符号长 CP，后 8 符号短 CP）
- 加窗配置：
  - 长 CP 过渡 `window_transition_long = 40`
  - 短 CP 过渡 `window_transition_short = 36`
- 符号数与比特数：
  - 总符号 `total_symbols = 100`
  - 每符号比特 `bits_per_symbol_total = carrier_count * 4`（16QAM）
- 信道 SNR：
  - `targetSNRdB = 15`（用于详细分析/单点 BER 展示）

提示：以上参数均可直接在脚本中修改后再次运行生效。

---

### 随机性与复现
- 本脚本每次运行都会随机生成新的发射比特与噪声：
  ```matlab
  rng('default'); rng('shuffle');
  baseband_out = randi([0 1], 1, baseband_out_length);
  ```
- 如需复现实验（固定随机序列），可将上述两行替换为固定种子，例如：
  ```matlab
  rng(20251116,'twister');  % 任意固定整数种子
  baseband_out = randi([0 1], 1, baseband_out_length);
  ```

---

### 帧结构与 CP 设计（为何有长/短 CP）
- CP 作用：循环前缀等效于在符号前复制尾部一段，用于对抗多径时延扩展并保持子载波正交，避免符号间干扰（ISI）和子载波间干扰（ICI）。
- LTE 正常 CP（Normal CP）中每个 0.5 ms 时隙含 7 个 OFDM 符号。为了在固定时隙长度内容纳 7 个符号并兼顾多径容限：
  - 第一个符号使用“较长 CP”（Long CP），以提供更大的时延容限，也便于时隙/子帧的定时对齐。
  - 其余符号使用“较短 CP”（Short CP），从而在时隙时长不变的前提下放入足够多的符号数。
- 在 10 MHz 配置下，脚本采用：`CP_long_samples = 160`，`CP_short_samples = 144`。这与 LTE 采样率 30.72 MHz 下的典型取值一致。
- 本脚本为“快速仿真”（`total_symbols = 10`），未完整搭建 1 ms 子帧的全部结构，但保留了“长/短 CP”的分配策略：前 2 个符号为长 CP，后续符号为短 CP（模拟两个时隙各自首符号较长 CP 的特性）。
- 加窗与重叠：脚本在 CP 范围内进行升余弦窗的平滑过渡与“重叠相加”，保证过渡仅发生在 CP 内，不影响有效符号主体数据：
  - 长 CP 对应的过渡长度 `window_transition_long = 40`
  - 短 CP 对应的过渡长度 `window_transition_short = 36`

---

### 时间结构层次：帧 - 子帧 - 时隙 - 符号 的关系
本脚本参考 LTE 正常 CP（Normal CP）的时间组织方式。四级结构与典型时长如下：

- 帧（Frame）
  - 无线帧时长为 10 ms。
  - 每帧包含 10 个子帧。
- 子帧（Subframe）
  - 子帧时长为 1 ms。
  - 每子帧包含 2 个时隙。
- 时隙（Slot）
  - 时隙时长为 0.5 ms。
  - 正常 CP 时，每个时隙包含 7 个 OFDM 符号；扩展 CP 时为 6 个。
  - 正常 CP 情况下，每个时隙的第 1 个符号使用较长 CP，其余符号使用较短 CP。
- 符号（OFDM Symbol）
  - 一个 OFDM 符号由“循环前缀（CP）+ 有效符号体（N 个采样）(+ 可选后缀/用于加窗过渡)”构成。
  - 正常 CP：在 30.72 MHz 采样率（10 MHz 系统）下，典型长 CP 约 160 点，短 CP 约 144 点；有效符号体为 IFFT 点数（本脚本为 2048）。

用包含关系表示：
- 1 帧 = 10 子帧 = 20 时隙
- 1 子帧 = 2 时隙
- 1 时隙（正常 CP）= 7 符号（第 1 个符号长 CP，其余短 CP）

与脚本配置的对应关系（快速仿真）
- 为了演示链路关键环节，脚本未逐层构造“10 ms 帧 → 1 ms 子帧 → 0.5 ms 时隙”的完整层级，而是采用 `total_symbols = 100` 的连续符号序列。
- 仍然保留了“每个时隙首符号长 CP、其余短 CP”的关键特征：在 100 个符号中，前 2 个符号使用长 CP，其余使用短 CP，可理解为“模拟两个连续时隙的首符号较长 CP”的特性在序列开始处得到体现。
- 若需要严格的帧/子帧/时隙层级组织，可将脚本中与层级相关的参数（如 `num_subframes、slots_per_subframe、symbols_per_slot`）用于构建多重循环，按层级拼接串行序列；当前脚本为了演示频谱/时域/星座与 BER 评估而进行了简化。

---

### 参数详解（test2.m 中的重要变量）
以下参数在脚本前半部分定义，决定系统配置、映射方式与可视化窗口。实际数值可在脚本中直接修改。

- 基本时间/频率与调制
  - `subcarrier_spacing = 15e3`：子载波间隔（Hz）。LTE 基本间隔为 15 kHz。
  - `fs = 30.72e6`：采样率（Hz）。对应 10 MHz 系统配置的典型采样率。
  - `symbol_duration = 1/subcarrier_spacing`：有效符号时长，约 66.67 μs。
  - `modulation_order = 16`、`bits_per_symbol = 4`：16QAM，每子载波承载 4 比特。

- 频域尺寸与子载波索引
  - `IFFT_bin_length = 2048`：IFFT 点数（N）。决定时域符号有效采样数。
  - 标准 k 索引范围：`-299..300`（DC 为 0）。为构造“实信号”时域输出，本脚本使用共轭对称：
    - 正频率数据子载波：`k_positive = 1:252`（252 个）
    - 负频率镜像载波：`k_negative = -251:-1`（251 个，用作共轭部分）
  - MATLAB IFFT bin 索引映射（从 1 开始）：
    - `map_k_to_bin = @(k) (k==0).*1 + (k>0).*(k+1) + (k<0).*(N+k+1)`
    - `carriers_positive = arrayfun(map_k_to_bin, k_positive)` → 正频率数据所在 bin（2..253）
    - `carriers_negative = arrayfun(map_k_to_bin, k_negative)` → 负频率镜像所在 bin（1798..2048）
  - 有效子载波计数：
    - `carrier_count = length(k_positive) = 252`（用于承载数据的正频率数量）
    - 含镜像的总有效子载波数为 252+251=503（不含 DC）
  - 空子载波（说明）：1 个 DC + 保护带（低频 48、高频 48）共 97 个。

- 帧/时隙/符号与快速仿真规模
  - `num_subframes = 1`、`slots_per_subframe = 2`、`symbols_per_slot = 7`：与 LTE 时间结构相关的参考值。
  - `total_symbols = 100`：采用 100 个连续符号，不严格展开完整帧-子帧-时隙层次，但在序列起始保留“长 CP”特性。

- 循环前缀（CP）与符号采样数
  - `CP_long_samples = 160`：长 CP（每个时隙第 1 个符号）。
  - `CP_short_samples = 144`：短 CP（其余符号）。
  - `CP_lengths`：长度数组，前 2 符号设为长 CP，其余为短 CP（模拟两个时隙的首符号较长 CP）。
  - 单符号总采样点（不含过渡后缀）：
    - 长 CP 符号：`symbol_samples_long = 2048 + 160 = 2208`
    - 短 CP 符号：`symbol_samples_short = 2048 + 144 = 2192`
  - `total_samples`：100 个符号的总采样数（用于规划内存或评估量级）。

- 加窗（升余弦窗）与过渡长度
  - 窗过渡长度与 CP 配合，保证“重叠相加”仅发生在 CP 范围内：
    - `window_transition_long = 40`（长 CP）
    - `window_transition_short = 36`（短 CP）
  - `window_transitions`：按符号顺序的过渡长度数组。
  - 窗函数生成：`rcos_win_full = rcoswindow(beta_local, N+CP_len)`，仅取前 `N+CP_len` 点用于 `CP+主体`。

- 比特/符号流水与矩阵维度
  - 每符号比特数：`bits_per_symbol_total = carrier_count * bits_per_symbol = 252 * 4 = 1008`
  - 总比特：`total_bits = total_symbols * bits_per_symbol_total = 100 * 1008 = 100800`
  - 发送随机比特：
    ```matlab
    rng('default'); rng('shuffle');
    baseband_out = randi([0 1], 1, baseband_out_length);
    ```
  - 16QAM 映射与重塑：
    - `complex_carrier_matrix = qam16(baseband_out);`
    - `complex_carrier_matrix = reshape(..., carrier_count, total_symbols)'` → 行为符号、列为载波

- 频域映射与实信号约束
  - 正频率：将 `complex_carrier_matrix(sym_idx, 1:carrier_count)` 放置到 `carriers_positive`。
  - 负频率：放置正频率的“反序共轭”到 `carriers_negative`，保证 IFFT 输出近似实值（数值误差下再取实部）。
  - DC（bin=1）设为 0。

- IFFT 与串行信号构造
  - `signal_after_IFFT = ifft(IFFT_modulation, IFFT_bin_length, 2);`
  - 对每符号构造结构 `[CP | 主体(N) | 后缀(window_trans)]`，后缀取自符号前部样本，用于与下一个符号 CP 的重叠。
  - 串行拼接时，“重叠相加”长度 `overlap_len = min(prev_window_trans, window_trans, CP_len)`，并仅限制在 CP 内。
  - 对比序列：`Tx_data_withoutwindow` 为未加窗、逐符号直接串接；`Tx_data` 为加窗+重叠相加的 LTE 风格串行。

- 可视化窗口与频谱
  - `zoom_len = 5000`：局部时域窗口长度（Figure 7/8/12/13 使用同一范围，便于对比）。
  - 使用 `fftshift(fft(...))` 生成双边谱，并归一化到自身最大值后转 dB。
  - 归一化频率轴：`(-N/2:(N/2-1))/N`，对应物理频率的 `±fs/2`。

- AWGN 信道与噪声参数
  - `targetSNRdB = 15`：分析用单点 SNR。
  - 噪声方差：`noise_sigma = var(Tx_data) / 10^(SNR/10)`；标准差 `sqrt(noise_sigma)`。
  - 接收信号：`Rx_data = Tx_data + noise`。
  - 局部窗口 SNR/MSE 打印：基于窗口内 `tx_seg` 与 `rx_seg` 的功率与误差。

- 接收端 FFT 与解调
  - 去除 CP/后缀，仅保留主体 `N` 点：构建 `Rx_data_matrix`，逐符号取 `CP_len+1 : CP_len+N`。
  - FFT 后取正频率数据 `carriers_positive`，再由极坐标转直角坐标得到 `Rx_complex_carrier_matrix`。
  - 判决：`demoduqam16`，并与 `baseband_out` 比较得到 BER。

- 图形与指标输出
  - Figures 1–13：详见“输出图形（主要）”一节。
  - 命令行打印包括：`total_symbols/total_bits/BER/Tx power/Noise variance/N/Active carriers/Null carriers/CP/window transition/solution method` 等。

---

### 输出图形（主要）
- Figure 1：IFFT 输入各频点幅度（展示子载波分配与空子载波）
- Figure 2：IFFT 输入相位（验证共轭对称）
- Figure 3：单符号时域（未加 CP/后缀）
- Figure 4：单符号时域（加 CP 与后缀）
- Figure 5：单符号时域（加窗后）
- Figure 6：发送端串行信号对比（未加窗逐符号串接 vs. LTE 风格重叠相加）
- Figure 7：局部窗口时域对比（未加窗 vs. 加窗）
- Figure 8：局部窗口频域对比（未加窗 vs. 加窗，双边谱）
- Figure 9：接收端 16QAM 星座图
- Figure 10：发/收比特流前 100 位对比
- Figure 11：单点 BER 显示（SNR = targetSNRdB）
- Figure 12：发送/接收信号局部窗口时域对比（并打印窗口内 SNR/MSE）
- Figure 13：与 Figure 12 同一窗口的双边谱对比

---

### 常见问题（FAQ）
- 运行报 “RNG 与旧生成器冲突”
  - 原因：在 MATLAB 会话中曾调用过 `rand('state',...)` / `rand('twister',...)`。
  - 解决：在运行脚本前调用 `rng('default');`，或重启 MATLAB 再运行。

- 星座图发散或 BER 异常高
  - 检查 `targetSNRdB` 是否过低；检查载波映射与 IFFT 长度是否匹配；确保依赖函数存在且无改动。

---

### 性能建议
- 将 `total_symbols` 适当增大以获得更稳定的 BER 统计，但运行时间会增加。
- 需要更快速度时，可关闭部分绘图或减少 `total_symbols`。

---

### 文件清单（与本脚本相关）
- `test2.m`：主仿真脚本（本文件说明对应）
- `qam16.m`：16QAM 调制
- `demoduqam16.m`：16QAM 解调
- `rcoswindow.m`：升余弦窗生成

---

如需进一步说明或扩展（如多径信道、均衡、编码/交织），可在现有流程基础上添加相应模块。祝你仿真顺利！


