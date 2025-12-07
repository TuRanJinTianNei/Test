%% 哈尔小波（Haar Wavelet）分解与重构仿真
% 详细演示尺度函数 φ 和小波函数 ψ 的作用
% 
% 尺度函数 φ: 用于捕捉信号的低频部分（近似信息/平均值）
% 小波函数 ψ: 用于捕捉信号的高频部分（细节信息/差异）

clear all;
close all;
clc;

fprintf('========== 哈尔小波（Haar Wavelet）分解与重构仿真 ==========\n\n');

%% ==================== 第一部分：定义哈尔尺度函数和小波函数 ====================
fprintf('第一部分：定义哈尔尺度函数 φ 和小波函数 ψ\n');

% 哈尔尺度函数 φ(t)
% φ(t) = 1, 0 ≤ t < 1
%      = 0, otherwise
phi_haar = @(t) (t >= 0 & t < 1) * 1;

% 哈尔小波函数 ψ(t)
% ψ(t) = 1,  0 ≤ t < 0.5
%      = -1, 0.5 ≤ t < 1
%      = 0,  otherwise
psi_haar = @(t) (t >= 0 & t < 0.5) * 1 + (t >= 0.5 & t < 1) * (-1);

% 可视化尺度函数和小波函数
figure('Position', [100, 100, 1200, 600], 'Name', '哈尔尺度函数和小波函数');

% 子图1：尺度函数 φ(t)
subplot(1, 2, 1);
t_phi = -0.5:0.001:1.5;
y_phi = arrayfun(phi_haar, t_phi);
plot(t_phi, y_phi, 'b-', 'LineWidth', 2);
xlabel('t', 'FontSize', 12);
ylabel('φ(t)', 'FontSize', 12);
title('哈尔尺度函数 φ(t)', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
xlim([-0.5, 1.5]);
ylim([-0.5, 1.5]);
text(0.5, 0.8, 'φ(t) = 1, 0 ≤ t < 1', 'FontSize', 11, 'HorizontalAlignment', 'center');

% 子图2：小波函数 ψ(t)
subplot(1, 2, 2);
t_psi = -0.5:0.001:1.5;
y_psi = arrayfun(psi_haar, t_psi);
plot(t_psi, y_psi, 'r-', 'LineWidth', 2);
xlabel('t', 'FontSize', 12);
ylabel('ψ(t)', 'FontSize', 12);
title('哈尔小波函数 ψ(t)', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
xlim([-0.5, 1.5]);
ylim([-1.5, 1.5]);
hold on;
plot([0.5, 0.5], [-1.5, 1.5], 'k--', 'LineWidth', 1);
text(0.25, 0.8, 'ψ(t) = 1', 'FontSize', 11, 'HorizontalAlignment', 'center');
text(0.75, -0.8, 'ψ(t) = -1', 'FontSize', 11, 'HorizontalAlignment', 'center');
hold off;

sgtitle('哈尔小波基础函数', 'FontSize', 16, 'FontWeight', 'bold');

fprintf('  ✓ 尺度函数 φ(t): 在 [0, 1) 上为 1，其他地方为 0\n');
fprintf('  ✓ 小波函数 ψ(t): 在 [0, 0.5) 上为 1，在 [0.5, 1) 上为 -1\n\n');

%% ==================== 第二部分：生成OFDM测试信号 ====================
fprintf('第二部分：生成OFDM测试信号\n');

% OFDM信号参数
N_subcarriers = 64;        % 子载波数量
N_cp = 16;                 % 循环前缀长度
M = 4;                     % QPSK调制
num_symbols = 4;           % OFDM符号数量（用于演示）
fs = 1e6;                  % 采样频率 1MHz

% 生成OFDM信号
% 1. 生成随机数据并进行QPSK调制
data = randi([0 M-1], N_subcarriers, num_symbols);
modulated_data = qammod(data, M, 'gray');

% 2. IFFT变换（将频域信号转换为时域）
ofdm_symbols = ifft(modulated_data, N_subcarriers);

% 3. 添加循环前缀
ofdm_symbols_cp = [ofdm_symbols(end-N_cp+1:end, :); ofdm_symbols];

% 4. 串行化，生成完整的OFDM信号
ofdm_signal = reshape(ofdm_symbols_cp, [], 1);

% 5. 归一化信号幅度（便于可视化）
ofdm_signal = ofdm_signal / max(abs(ofdm_signal)) * 2;

% 6. 创建时间轴
N_signal = length(ofdm_signal);
t_signal = (0:N_signal-1) / fs;  % 时间轴（秒）

% 使用OFDM信号作为测试信号
f_signal_normalized = real(ofdm_signal);  % 使用实部（OFDM信号通常是复数）

fprintf('  OFDM参数:\n');
fprintf('    子载波数量: %d\n', N_subcarriers);
fprintf('    循环前缀长度: %d\n', N_cp);
fprintf('    调制方式: QPSK\n');
fprintf('    OFDM符号数: %d\n', num_symbols);
fprintf('  信号长度: %d 个采样点\n', N_signal);
fprintf('  信号时长: %.4f 秒\n', max(t_signal));
fprintf('  信号范围: [%.2f, %.2f]\n', min(f_signal_normalized), max(f_signal_normalized));
fprintf('  信号特点: OFDM多载波信号，包含多个子载波频率成分\n\n');

%% ==================== 第三部分：哈尔小波分解 ====================
fprintf('第三部分：哈尔小波分解\n');
fprintf('  使用MATLAB内置的Haar小波进行分解\n\n');

% 使用MATLAB的Haar小波进行分解
wavelet_name = 'haar';
decomposition_level = 3;  % 分解3层

% 执行小波分解
[c, l] = wavedec(f_signal_normalized, decomposition_level, wavelet_name);

% 提取近似系数和细节系数
approximation = appcoef(c, l, wavelet_name);  % 近似系数（尺度函数系数）
details = cell(decomposition_level, 1);
for i = 1:decomposition_level
    details{i} = detcoef(c, l, i);  % 细节系数（小波函数系数）
end

fprintf('  分解层数: %d\n', decomposition_level);
fprintf('  原始信号长度: %d 个采样点\n', N_signal);
fprintf('  近似系数长度: %d (尺度函数 φ 的系数)\n', length(approximation));
for i = 1:decomposition_level
    fprintf('  第%d层细节系数长度: %d (小波函数 ψ 的系数)\n', i, length(details{i}));
end
fprintf('\n');

%% ==================== 系数单位和长度说明 ====================
fprintf('========== 系数单位和长度说明 ==========\n');
fprintf('\n1. 系数的单位：\n');
fprintf('   • 近似系数 c_{j,k} 的单位：与原始信号 f(t) 的单位相同\n');
fprintf('     对于OFDM信号，单位是：幅度（无量纲，或V/归一化单位）\n');
fprintf('   • 细节系数 d_{j,k} 的单位：与原始信号 f(t) 的单位相同\n');
fprintf('     对于OFDM信号，单位是：幅度（无量纲，或V/归一化单位）\n');
fprintf('   • 原因：系数是通过内积计算得到的\n');
fprintf('     c_{j,k} = <f(t), φ_{j,k}(t)> = ∫ f(t) · φ_{j,k}(t) dt\n');
fprintf('     d_{j,k} = <f(t), ψ_{j,k}(t)> = ∫ f(t) · ψ_{j,k}(t) dt\n');
fprintf('     由于 φ 和 ψ 是无量纲的，所以系数与 f(t) 的单位相同\n\n');

fprintf('2. 为什么长度不一样？\n');
fprintf('   原因：小波分解过程中的降采样（Downsampling）\n\n');
fprintf('   分解过程：\n');
fprintf('   第1层分解：\n');
fprintf('     • 原始信号（长度 N = %d）\n', N_signal);
fprintf('     ↓ [低通滤波 + 降采样（÷2）]\n');
fprintf('     • 第1层近似（长度 N/2 = %d）\n', floor(N_signal/2));
fprintf('     • 第1层细节（长度 N/2 = %d）\n', floor(N_signal/2));
fprintf('   \n');
fprintf('   第2层分解（对第1层近似继续分解）：\n');
fprintf('     • 第1层近似（长度 N/2 = %d）\n', floor(N_signal/2));
fprintf('     ↓ [低通滤波 + 降采样（÷2）]\n');
fprintf('     • 第2层近似（长度 N/4 = %d）\n', floor(N_signal/4));
fprintf('     • 第2层细节（长度 N/4 = %d）\n', floor(N_signal/4));
fprintf('   \n');
fprintf('   第3层分解（对第2层近似继续分解）：\n');
fprintf('     • 第2层近似（长度 N/4 = %d）\n', floor(N_signal/4));
fprintf('     ↓ [低通滤波 + 降采样（÷2）]\n');
fprintf('     • 第3层近似（长度 N/8 = %d）≈ 近似系数长度 %d\n', floor(N_signal/8), length(approximation));
fprintf('     • 第3层细节（长度 N/8 = %d）≈ 第3层细节系数长度 %d\n', floor(N_signal/8), length(details{3}));
fprintf('   \n');
fprintf('   长度关系（二进制规律）：\n');
fprintf('     • 原始信号：N = %d\n', N_signal);
fprintf('     • 第1层细节 D1：N/2 = %d\n', length(details{1}));
fprintf('     • 第2层细节 D2：N/4 = %d\n', length(details{2}));
fprintf('     • 第3层细节 D3：N/8 = %d\n', length(details{3}));
fprintf('     • 近似系数 A3：N/8 = %d\n', length(approximation));
fprintf('   \n');
fprintf('   验证：D1 + D2 + D3 + A3 = %d + %d + %d + %d = %d ≈ N = %d ✓\n', ...
    length(details{1}), length(details{2}), length(details{3}), length(approximation), ...
    length(details{1}) + length(details{2}) + length(details{3}) + length(approximation), N_signal);
fprintf('   \n');
fprintf('   为什么需要降采样？\n');
fprintf('     • 每层分解后，频率范围减半（奈奎斯特定理）\n');
fprintf('     • 降采样可以去除冗余信息，保持信息量不变\n');
fprintf('     • 降采样后，采样率减半，需要的采样点数也减半\n');
fprintf('     • 这是多分辨率分析（MRA）的核心特性\n\n');
fprintf('   时间覆盖范围：\n');
signal_duration = max(t_signal);  % 计算信号时长
fprintf('     • 重要：所有系数都覆盖完整的时间范围（0 到 %.4f 秒）\n', signal_duration);
fprintf('     • 长度不同是因为采样密度不同（采样率不同）\n');
fprintf('     • D1：最高采样率，时间分辨率最高\n');
fprintf('     • A3：最低采样率，时间分辨率最低\n\n');
fprintf('==========================================\n\n');

%% ==================== 详细解释：分解过程 ====================
fprintf('========== 详细解释：分解过程（用具体例子）==========\n');
fprintf('\n假设原始信号有 8 个采样点：[x1, x2, x3, x4, x5, x6, x7, x8]\n\n');

fprintf('【第1层分解】\n');
fprintf('输入：原始信号，8个点\n');
fprintf('操作：\n');
fprintf('  1. 低通滤波：提取低频成分（平滑/平均）\n');
fprintf('  2. 高通滤波：提取高频成分（差异/变化）\n');
fprintf('  3. 降采样：每2个点保留1个点（采样点数÷2）\n');
fprintf('输出：\n');
fprintf('  • 近似A1：4个点 [a1, a2, a3, a4]  ← 低频，平滑后的信号\n');
fprintf('  • 细节D1：4个点 [d1, d2, d3, d4]  ← 高频，变化的信息\n');
fprintf('说明：A1和D1各4个点，加起来还是8个点（信息量守恒）\n\n');

fprintf('【第2层分解】\n');
fprintf('输入：第1层的近似A1（4个点）[a1, a2, a3, a4]\n');
fprintf('注意：只对近似继续分解，细节D1不再分解！\n');
fprintf('操作：\n');
fprintf('  1. 对A1进行低通滤波：提取更低频成分\n');
fprintf('  2. 对A1进行高通滤波：提取中频成分\n');
fprintf('  3. 降采样：每2个点保留1个点（采样点数÷2）\n');
fprintf('输出：\n');
fprintf('  • 近似A2：2个点 [a1, a2]  ← 更低频，更平滑\n');
fprintf('  • 细节D2：2个点 [d1, d2]  ← 中频，中等变化\n');
fprintf('说明：A2和D2各2个点，加起来是4个点（等于A1的长度）\n\n');

fprintf('【第3层分解】\n');
fprintf('输入：第2层的近似A2（2个点）[a1, a2]\n');
fprintf('操作：\n');
fprintf('  1. 对A2进行低通滤波：提取最低频成分\n');
fprintf('  2. 对A2进行高通滤波：提取低频变化\n');
fprintf('  3. 降采样：每2个点保留1个点（采样点数÷2）\n');
fprintf('输出：\n');
fprintf('  • 近似A3：1个点 [a1]  ← 最低频，整体平均值\n');
fprintf('  • 细节D3：1个点 [d1]  ← 低频变化\n');
fprintf('说明：A3和D3各1个点，加起来是2个点（等于A2的长度）\n\n');

fprintf('【最终结果】\n');
fprintf('原始信号（8点）被分解为：\n');
fprintf('  • D1（4点）：最高频细节\n');
fprintf('  • D2（2点）：中频细节\n');
fprintf('  • D3（1点）：低频细节\n');
fprintf('  • A3（1点）：最低频近似（整体平均值）\n');
fprintf('总计：4 + 2 + 1 + 1 = 8点 ✓（信息量守恒）\n\n');

fprintf('【关键理解】\n');
fprintf('1. 为什么长度减半？\n');
fprintf('   → 因为降采样：每2个点只保留1个点\n');
fprintf('   → 频率范围减半，需要的采样点数也减半\n\n');
fprintf('2. 为什么只对近似继续分解？\n');
fprintf('   → 近似包含低频信息，可以进一步分解\n');
fprintf('   → 细节包含高频信息，已经是最细的细节了\n\n');
fprintf('3. 为什么信息量守恒？\n');
fprintf('   → 每层分解：输入点数 = 近似点数 + 细节点数\n');
fprintf('   → 所以：D1 + D2 + D3 + A3 = 原始信号长度\n\n');
fprintf('4. 时间覆盖范围？\n');
fprintf('   → 所有系数都覆盖完整时间范围\n');
fprintf('   → 长度不同 = 采样密度不同（时间分辨率不同）\n');
fprintf('   → D1：最密集，时间分辨率最高\n');
fprintf('   → A3：最稀疏，时间分辨率最低\n\n');
fprintf('==========================================\n\n');

%% ==================== 可视化：系数长度差异说明 ====================
figure('Position', [100, 100, 1400, 700], 'Name', '系数长度差异说明');

% 子图1：长度对比条形图
subplot(2, 2, 1);
lengths = [N_signal, length(details{1}), length(details{2}), length(details{3}), length(approximation)];
names = {'原始信号', 'D1', 'D2', 'D3', 'A3'};
bar(lengths, 'FaceColor', [0.3 0.6 0.9]);
set(gca, 'XTickLabel', names);
ylabel('长度（采样点数）', 'FontSize', 11);
title('各层系数长度对比', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
% 添加数值标签
for i = 1:length(lengths)
    text(i, lengths(i), sprintf('%d', lengths(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10);
end

% 子图2：长度关系（对数刻度）
subplot(2, 2, 2);
semilogy([1, 2, 3, 4, 5], lengths, 'o-', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'r');
set(gca, 'XTick', 1:5);
set(gca, 'XTickLabel', names);
ylabel('长度（对数刻度）', 'FontSize', 11);
title('长度关系：每层减半（二进制）', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

% 子图3：降采样过程示意图
subplot(2, 2, 3);
axis off;
text_str = {
    '降采样过程示意图：';
    '';
    sprintf('原始信号: %d 点', N_signal);
    '    ↓ [低通滤波 + 降采样 ÷2]';
    sprintf('第1层: 近似 %d 点 + 细节 %d 点', floor(N_signal/2), floor(N_signal/2));
    '    ↓ [对近似继续分解，降采样 ÷2]';
    sprintf('第2层: 近似 %d 点 + 细节 %d 点', floor(N_signal/4), floor(N_signal/4));
    '    ↓ [对近似继续分解，降采样 ÷2]';
    sprintf('第3层: 近似 %d 点 + 细节 %d 点', length(approximation), length(details{3}));
    '';
    '关键点：';
    '  • 每层分解后，采样点数减半';
    '  • 频率范围也减半（奈奎斯特定理）';
    '  • 但时间覆盖范围不变';
    '  • 总信息量守恒：N = D1 + D2 + D3 + A3';
};
text(0.1, 0.5, text_str, 'FontSize', 11, 'VerticalAlignment', 'middle', ...
    'FontName', 'FixedWidth');

% 子图4：时间覆盖范围说明
subplot(2, 2, 4);
axis off;
text_str2 = {
    '时间覆盖范围说明：';
    '';
    '所有系数都覆盖完整时间范围：';
    sprintf('  0 到 %.4f 秒（%.2f ms）', signal_duration, signal_duration*1000);
    '';
    '采样密度对比：';
    sprintf('  • D1: %d 点 / %.4f秒 = %.1f 点/秒', ...
        length(details{1}), signal_duration, length(details{1})/signal_duration);
    sprintf('  • D2: %d 点 / %.4f秒 = %.1f 点/秒', ...
        length(details{2}), signal_duration, length(details{2})/signal_duration);
    sprintf('  • D3: %d 点 / %.4f秒 = %.1f 点/秒', ...
        length(details{3}), signal_duration, length(details{3})/signal_duration);
    sprintf('  • A3: %d 点 / %.4f秒 = %.1f 点/秒', ...
        length(approximation), signal_duration, length(approximation)/signal_duration);
    '';
    '结论：';
    '  • 长度不同 = 采样密度不同';
    '  • 时间范围相同 = 都覆盖完整信号';
    '  • D1时间分辨率最高，A3时间分辨率最低';
};
text(0.1, 0.5, text_str2, 'FontSize', 10, 'VerticalAlignment', 'middle', ...
    'FontName', 'FixedWidth');

sgtitle('系数长度差异原因：降采样（Downsampling）', 'FontSize', 14, 'FontWeight', 'bold');

%% ==================== 第四部分：可视化分解结果 ====================
fprintf('第四部分：可视化分解结果\n\n');

figure('Position', [150, 150, 1400, 1000], 'Name', '哈尔小波分解结果');

% 计算各层的时间轴（映射到原始信号的时间范围）
signal_duration = max(t_signal);
t_approx = linspace(0, signal_duration, length(approximation));
t_details = cell(decomposition_level, 1);
for i = 1:decomposition_level
    t_details{i} = linspace(0, signal_duration, length(details{i}));
end

% 子图1：原始信号
subplot(decomposition_level + 2, 1, 1);
plot(t_signal * 1000, f_signal_normalized, 'b-', 'LineWidth', 2);  % 转换为毫秒显示
xlabel('时间 (ms)', 'FontSize', 11);
ylabel('幅度', 'FontSize', 11);
title('原始OFDM信号 f(t)', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
xlim([0, signal_duration * 1000]);
set(gca, 'XTickLabel', []);

% 子图2：近似系数（尺度函数 φ 的系数）
subplot(decomposition_level + 2, 1, 2);
stem(t_approx * 1000, approximation, 'r-', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('时间 (ms)', 'FontSize', 11);
ylabel('近似系数 c_{j,k}', 'FontSize', 11);
title(sprintf('近似系数（尺度函数 φ 的系数）- 长度=%d', length(approximation)), ...
    'FontSize', 12, 'FontWeight', 'bold');
grid on;
xlim([0, signal_duration * 1000]);
set(gca, 'XTickLabel', []);
text(signal_duration * 500, max(approximation)*0.8, '这些系数表示OFDM信号在不同尺度下的平均值（平滑趋势）', ...
    'FontSize', 10, 'HorizontalAlignment', 'center', 'BackgroundColor', 'yellow');

% 子图3-5：各层细节系数（小波函数 ψ 的系数）
colors = {'g', 'm', 'c'};
for i = 1:decomposition_level
    subplot(decomposition_level + 2, 1, i + 2);
    stem(t_details{i} * 1000, details{i}, 'Color', colors{i}, 'LineWidth', 1.5, 'MarkerSize', 6);
    xlabel('时间 (ms)', 'FontSize', 11);
    ylabel(sprintf('细节系数 d_{%d,k}', i), 'FontSize', 11);
    title(sprintf('第%d层细节系数（小波函数 ψ 的系数）- 长度=%d', i, length(details{i})), ...
        'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    xlim([0, signal_duration * 1000]);
    if i < decomposition_level
        set(gca, 'XTickLabel', []);
    end
    text(signal_duration * 500, max(abs(details{i}))*0.7, '这些系数表示OFDM信号在该尺度下的差异（变化细节）', ...
        'FontSize', 10, 'HorizontalAlignment', 'center', 'BackgroundColor', 'yellow');
end

sgtitle('哈尔小波分解：尺度函数 φ 与小波函数 ψ 的作用', 'FontSize', 16, 'FontWeight', 'bold');

%% ==================== 第五部分：内积计算演示 ====================
fprintf('第五部分：内积计算演示\n');
fprintf('  演示如何通过内积计算尺度系数和小波系数\n\n');

% 选择信号的一个区间进行演示（选择信号的前一部分，归一化到[0, 1)区间）
% 选择一个OFDM符号的长度作为演示区间
demo_length = N_subcarriers + N_cp;  % 一个OFDM符号的长度
if demo_length > length(f_signal_normalized)
    demo_length = floor(length(f_signal_normalized) / 2);
end
t_demo_original = t_signal(1:demo_length);
f_demo_original = f_signal_normalized(1:demo_length);

% 归一化到[0, 1)区间用于演示
t_demo_normalized = (t_demo_original - t_demo_original(1)) / (t_demo_original(end) - t_demo_original(1));
t_demo = 0:0.001:0.999;
f_demo = interp1(t_demo_normalized, f_demo_original, t_demo, 'linear', 'extrap');

% 计算与尺度函数 φ_0 的内积（整个区间的平均值）
c_00 = trapz(t_demo, f_demo .* arrayfun(phi_haar, t_demo));
fprintf('  尺度系数 c_{0,0} = <f, φ_0> = %.4f\n', c_00);
fprintf('    物理意义：OFDM信号在该区间上的平均值（平滑趋势）\n');

% 计算与小波函数 ψ_0 的内积（前后半段的差异）
d_00 = trapz(t_demo, f_demo .* arrayfun(psi_haar, t_demo));
fprintf('  小波系数 d_{0,0} = <f, ψ_0> = %.4f\n', d_00);
fprintf('    物理意义：OFDM信号前半段平均值 - 后半段平均值（变化细节）\n\n');

% 可视化内积计算过程
figure('Position', [200, 200, 1400, 800], 'Name', '内积计算演示');

% 子图1：信号和尺度函数
subplot(2, 2, 1);
plot(t_demo, f_demo, 'b-', 'LineWidth', 2, 'DisplayName', 'f(t)');
hold on;
plot(t_demo, arrayfun(phi_haar, t_demo), 'r-', 'LineWidth', 2, 'DisplayName', 'φ_0(t)');
fill([t_demo, fliplr(t_demo)], [f_demo .* arrayfun(phi_haar, t_demo), zeros(size(t_demo))], ...
    'g', 'FaceAlpha', 0.3, 'DisplayName', 'f(t)·φ_0(t)');
xlabel('t', 'FontSize', 11);
ylabel('幅度', 'FontSize', 11);
title(sprintf('内积计算：c_{0,0} = <f, φ_0> = %.4f', c_00), 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best');
grid on;
hold off;

% 子图2：信号和小波函数
subplot(2, 2, 2);
plot(t_demo, f_demo, 'b-', 'LineWidth', 2, 'DisplayName', 'f(t)');
hold on;
plot(t_demo, arrayfun(psi_haar, t_demo), 'r-', 'LineWidth', 2, 'DisplayName', 'ψ_0(t)');
fill([t_demo, fliplr(t_demo)], [f_demo .* arrayfun(psi_haar, t_demo), zeros(size(t_demo))], ...
    'g', 'FaceAlpha', 0.3, 'DisplayName', 'f(t)·ψ_0(t)');
plot([0.5, 0.5], [-1, 4], 'k--', 'LineWidth', 1);
xlabel('t', 'FontSize', 11);
ylabel('幅度', 'FontSize', 11);
title(sprintf('内积计算：d_{0,0} = <f, ψ_0> = %.4f', d_00), 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best');
grid on;
hold off;

% 子图3：尺度函数的作用（求平均）
subplot(2, 2, 3);
f_mean = mean(f_demo) * ones(size(t_demo));
plot(t_demo, f_demo, 'b-', 'LineWidth', 1.5, 'DisplayName', '原始信号');
hold on;
plot(t_demo, f_mean, 'r--', 'LineWidth', 2, 'DisplayName', sprintf('平均值 = %.4f', mean(f_demo)));
xlabel('t', 'FontSize', 11);
ylabel('幅度', 'FontSize', 11);
title('尺度函数 φ 的作用：求平均值（近似）', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best');
grid on;
hold off;

% 子图4：小波函数的作用（求差异）
subplot(2, 2, 4);
f_first_half = f_demo(1:500);
f_second_half = f_demo(500:end);
mean_first = mean(f_first_half);
mean_second = mean(f_second_half);
plot(t_demo(1:500), f_first_half, 'b-', 'LineWidth', 1.5, 'DisplayName', '前半段');
hold on;
plot(t_demo(500:end), f_second_half, 'g-', 'LineWidth', 1.5, 'DisplayName', '后半段');
plot([0, 0.5], [mean_first, mean_first], 'r--', 'LineWidth', 2, 'DisplayName', sprintf('前半段平均 = %.4f', mean_first));
plot([0.5, 1], [mean_second, mean_second], 'm--', 'LineWidth', 2, 'DisplayName', sprintf('后半段平均 = %.4f', mean_second));
plot([0.5, 0.5], [-1, 4], 'k--', 'LineWidth', 1);
xlabel('t', 'FontSize', 11);
ylabel('幅度', 'FontSize', 11);
title(sprintf('小波函数 ψ 的作用：求差异 = %.4f - %.4f = %.4f', ...
    mean_first, mean_second, mean_first - mean_second), 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best');
grid on;
hold off;

sgtitle('内积计算演示：OFDM信号的尺度函数 φ 与小波函数 ψ 的作用', 'FontSize', 16, 'FontWeight', 'bold');

%% ==================== 第六部分：信号重构 ====================
fprintf('第六部分：信号重构\n');
fprintf('  使用近似系数和细节系数重构原始信号\n\n');

% 使用所有系数重构
reconstructed_full = waverec(c, l, wavelet_name);

% 只使用近似系数重构（去除所有细节）
c_approx_only = c;
% 将细节系数置零
approx_length = length(approximation);
c_approx_only(approx_length+1:end) = 0;
reconstructed_approx_only = waverec(c_approx_only, l, wavelet_name);

% 计算重构误差
reconstruction_error = f_signal_normalized(1:length(reconstructed_full)) - reconstructed_full;

fprintf('  完整重构误差（最大）: %.6f\n', max(abs(reconstruction_error)));
fprintf('  完整重构误差（均方）: %.6f\n', mean(reconstruction_error.^2));
fprintf('  信噪比: %.2f dB\n\n', 10*log10(var(f_signal_normalized(1:length(reconstructed_full)))/var(reconstruction_error)));

% 可视化重构结果
figure('Position', [250, 250, 1400, 800], 'Name', '信号重构');

% 子图1：原始信号 vs 完整重构
subplot(2, 2, 1);
plot(t_signal * 1000, f_signal_normalized, 'b-', 'LineWidth', 2, 'DisplayName', '原始OFDM信号');
hold on;
plot(t_signal(1:length(reconstructed_full)) * 1000, reconstructed_full, 'r--', 'LineWidth', 2, 'DisplayName', '完整重构');
xlabel('时间 (ms)', 'FontSize', 11);
ylabel('幅度', 'FontSize', 11);
title('原始OFDM信号 vs 完整重构（近似+细节）', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best');
grid on;
hold off;

% 子图2：原始信号 vs 仅近似重构
subplot(2, 2, 2);
plot(t_signal * 1000, f_signal_normalized, 'b-', 'LineWidth', 2, 'DisplayName', '原始OFDM信号');
hold on;
plot(t_signal(1:length(reconstructed_approx_only)) * 1000, reconstructed_approx_only, 'g--', 'LineWidth', 2, 'DisplayName', '仅近似重构（无细节）');
xlabel('时间 (ms)', 'FontSize', 11);
ylabel('幅度', 'FontSize', 11);
title('原始OFDM信号 vs 仅近似重构（只有尺度函数系数）', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best');
grid on;
hold off;

% 子图3：重构误差
subplot(2, 2, 3);
plot(t_signal(1:length(reconstruction_error)) * 1000, reconstruction_error, 'r-', 'LineWidth', 1.5);
xlabel('时间 (ms)', 'FontSize', 11);
ylabel('误差', 'FontSize', 11);
title('重构误差（原始 - 完整重构）', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

% 子图4：系数能量分布
subplot(2, 2, 4);
energy_approx = sum(approximation.^2);
energy_details = zeros(decomposition_level, 1);
for i = 1:decomposition_level
    energy_details(i) = sum(details{i}.^2);
end
bar([energy_approx; energy_details], 'FaceColor', [0.3 0.6 0.9]);
set(gca, 'XTickLabel', {'近似系数', '细节1', '细节2', '细节3'});
ylabel('能量', 'FontSize', 11);
title('系数能量分布', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

sgtitle('OFDM信号重构：使用尺度函数系数和小波函数系数', 'FontSize', 16, 'FontWeight', 'bold');

%% ==================== 第七部分：总结 ====================
fprintf('========== 总结 ==========\n');
fprintf('对于OFDM信号，哈尔小波分解的作用：\n\n');
fprintf('尺度函数 φ 的作用：\n');
fprintf('  • 计算内积 <f, φ_{j,k}> 得到近似系数 c_{j,k}\n');
fprintf('  • 近似系数表示OFDM信号在不同尺度下的平均值（平滑趋势）\n');
fprintf('  • 捕获OFDM信号的低频、平滑信息（整体包络）\n\n');

fprintf('小波函数 ψ 的作用：\n');
fprintf('  • 计算内积 <f, ψ_{j,k}> 得到细节系数 d_{j,k}\n');
fprintf('  • 细节系数表示OFDM信号在不同尺度下的差异（变化细节）\n');
fprintf('  • 捕获OFDM信号的高频、突变信息（子载波变化、符号边界等）\n\n');

fprintf('哈尔分解公式：\n');
fprintf('  f(t) = Σ c_{j,k} · φ_{j,k}(t) + Σ d_{j,k} · ψ_{j,k}(t)\n');
fprintf('        = 近似部分（平滑趋势）+ 细节部分（变化细节）\n\n');

fprintf('应用：\n');
fprintf('  • 可以用于OFDM信号的去噪（去除细节系数中的噪声）\n');
fprintf('  • 可以用于OFDM信号的压缩（保留主要近似系数）\n');
fprintf('  • 可以用于OFDM信号的边缘检测（细节系数大的地方）\n\n');

fprintf('========== 仿真完成 ==========\n');

