%===============================================================================
% signalGenerate_Test.m - 循环谱分析测试脚本
% 
% 功能说明:
%   1. 加载 signalGenerate.m 生成的 RX_signal.mat 接收信号
%   2. 对接收信号进行循环谱分析
%   3. 提取 alpha = 0 的切片（即功率谱密度）
%
% 循环谱定义:
%   S_x^α(f) = lim E[X_T(f + α/2) * X_T^*(f - α/2)]
%   其中 X_T(f) 是信号的傅里叶变换
%
% 当 alpha = 0 时:
%   S_x^0(f) = |X(f)|^2  (功率谱密度)
%
% 创建日期: 2025.12.23
%===============================================================================

clc;
clear all;
close all;

fprintf('========================================\n');
fprintf('循环谱分析测试程序\n');
fprintf('========================================\n\n');

%===============================================================================
% 【步骤1: 加载接收信号】
%===============================================================================
fprintf('步骤1: 加载接收信号...\n');

% 检查文件是否存在（先检查当前目录，再检查上一级目录）
if exist('RX_signal.mat', 'file')
    mat_file = 'RX_signal.mat';
elseif exist('../RX_signal.mat', 'file')
    mat_file = '../RX_signal.mat';
else
    error('错误：找不到文件 RX_signal.mat，请先运行 signalGenerate.m 生成信号');
end

% 加载信号数据
load(mat_file);

% 检查必要的变量是否存在
if ~exist('Rx_data', 'var')
    error('错误：RX_signal.mat 中未找到 Rx_data 变量');
end

if ~exist('fs', 'var')
    error('错误：RX_signal.mat 中未找到 fs 变量');
end

fprintf('  - 接收信号长度: %d 样本\n', length(Rx_data));
fprintf('  - 采样频率: %.2f MHz\n', fs/1e6);
fprintf('  - 信号时长: %.4f ms\n', length(Rx_data)/fs*1e3);
fprintf('  完成！\n\n');

%===============================================================================
% 【步骤2: 循环谱分析参数设置】
%===============================================================================
fprintf('步骤2: 设置循环谱分析参数...\n');

% 循环谱分析参数
Nfft = 2^12;                    % FFT点数（频率方向，控制频率分辨率）
Nfft_alpha = 2^10;              % FFT点数（循环频率方向，控制alpha分辨率）
overlap_ratio = 0.5;            % 重叠比例
window_length = min(2^11, floor(length(Rx_data)/8));  % 窗长度

% 确保窗长度不超过信号长度
if window_length > length(Rx_data)
    window_length = floor(length(Rx_data) / 2);
end

% 计算重叠样本数
noverlap = floor(window_length * overlap_ratio);

% 计算分段数
nsegments = max(1, floor((length(Rx_data) - noverlap) / (window_length - noverlap)));

fprintf('  - FFT点数（频率方向）: %d\n', Nfft);
fprintf('  - FFT点数（循环频率方向）: %d\n', Nfft_alpha);
fprintf('  - 窗长度: %d 样本\n', window_length);
fprintf('  - 重叠样本数: %d\n', noverlap);
fprintf('  - 分段数: %d\n', nsegments);
fprintf('  完成！\n\n');

%===============================================================================
% 【步骤3: 计算循环谱（频域方法）】
%===============================================================================
fprintf('步骤3: 计算循环谱（使用频域方法）...\n');

% 初始化循环谱矩阵
% 循环谱 S_x^alpha(f) 的维度：[频率 f × 循环频率 alpha]
cyclic_spectrum = zeros(Nfft, Nfft_alpha);

% 使用汉宁窗
window = hanning(window_length);
window = window(:)';  % 确保是行向量

% 频率轴（零频率居中）
f_axis = (-Nfft/2:(Nfft/2-1)) * fs / Nfft;

% 循环频率轴（零频率居中）
alpha_axis = (-Nfft_alpha/2:(Nfft_alpha/2-1)) * fs / Nfft_alpha;

% 对信号进行分段处理并计算循环谱
fprintf('  正在处理 %d 个分段...\n', nsegments);
for seg = 1:nsegments
    % 计算当前段的起始和结束位置
    start_idx = (seg - 1) * (window_length - noverlap) + 1;
    end_idx = start_idx + window_length - 1;
    
    if end_idx > length(Rx_data)
        break;
    end
    
    % 提取当前段并加窗
    x_segment = Rx_data(start_idx:end_idx) .* window;
    
    % 计算FFT（零频率居中）
    X_f = fft(x_segment, Nfft);
    X_f_shifted = fftshift(X_f);
    
    % 计算循环谱：S_x^alpha(f) = X(f + alpha/2) * X^*(f - alpha/2)
    for alpha_idx = 1:Nfft_alpha
        alpha = alpha_axis(alpha_idx);
        
        for f_idx = 1:Nfft
            f = f_axis(f_idx);
            
            % 计算 f + alpha/2 和 f - alpha/2
            f1 = f + alpha/2;
            f2 = f - alpha/2;
            
            % 转换为索引（考虑fftshift，零频率在中间）
            % f_axis的范围是 [-fs/2, fs/2)，索引1对应-fs/2，索引Nfft/2+1对应0
            idx1 = round((f1 + fs/2) * Nfft / fs) + 1;
            idx2 = round((f2 + fs/2) * Nfft / fs) + 1;
            
            % 边界检查
            if idx1 >= 1 && idx1 <= Nfft && idx2 >= 1 && idx2 <= Nfft
                cyclic_spectrum(f_idx, alpha_idx) = cyclic_spectrum(f_idx, alpha_idx) + ...
                    X_f_shifted(idx1) * conj(X_f_shifted(idx2));
            end
        end
    end
    
    % 显示进度
    if mod(seg, max(1, floor(nsegments/10))) == 0 || seg == nsegments
        fprintf('    进度: %d/%d (%.1f%%)\n', seg, nsegments, 100*seg/nsegments);
    end
end

% 平均化（除以分段数）
cyclic_spectrum = cyclic_spectrum / nsegments;

fprintf('  完成！\n\n');

%===============================================================================
% 【步骤4: 提取 alpha = 0 的切片】
%===============================================================================
fprintf('步骤4: 提取 alpha = 0 的切片...\n');

% 找到 alpha = 0 对应的索引
[~, alpha_zero_idx] = min(abs(alpha_axis));  % 找到最接近0的索引

% 提取 alpha = 0 的切片（功率谱密度）
PSD_alpha_zero = cyclic_spectrum(:, alpha_zero_idx);

% 转换为dB
PSD_alpha_zero_dB = 10*log10(abs(PSD_alpha_zero) + eps);

fprintf('  - Alpha = 0 对应的循环频率索引: %d\n', alpha_zero_idx);
fprintf('  - 对应的循环频率: %.6f Hz\n', alpha_axis(alpha_zero_idx));
fprintf('  完成！\n\n');

%===============================================================================
% 【步骤5: 使用Welch方法计算PSD作为对比】
%===============================================================================
fprintf('步骤5: 使用Welch方法计算PSD作为对比...\n');

[Pxx_welch, f_welch] = pwelch(Rx_data, window_length, noverlap, Nfft, fs);
Pxx_welch_dB = 10*log10(Pxx_welch + eps);

% 将Welch PSD的频率轴转换为零频率居中
f_welch_shifted = fftshift(f_welch - fs/2);
Pxx_welch_dB_shifted = fftshift(Pxx_welch_dB);

fprintf('  完成！\n\n');

%===============================================================================
% 【步骤6: 可视化结果】
%===============================================================================
fprintf('步骤6: 绘制循环谱分析结果...\n');

% 图1: 循环谱分析结果
figure('Name', '循环谱分析结果', 'Position', [100, 100, 1400, 900]);

% 子图1: 循环谱的2D等高线图
subplot(2, 3, 1);
cyclic_spectrum_dB = 10*log10(abs(cyclic_spectrum) + eps);
imagesc(alpha_axis/1e3, f_axis/1e6, cyclic_spectrum_dB);
colorbar;
xlabel('循环频率 α (kHz)', 'FontSize', 11);
ylabel('频率 f (MHz)', 'FontSize', 11);
title('循环谱 |S_x^α(f)| (dB)', 'FontSize', 12, 'FontWeight', 'bold');
axis xy;
colormap('jet');
caxis([max(cyclic_spectrum_dB(:))-60, max(cyclic_spectrum_dB(:))]);

% 子图2: Alpha = 0 的切片（功率谱密度）
subplot(2, 3, 2);
plot(f_axis/1e6, PSD_alpha_zero_dB, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('频率 f (MHz)', 'FontSize', 11);
ylabel('功率谱密度 (dB)', 'FontSize', 11);
title('Alpha = 0 切片（功率谱密度）', 'FontSize', 12, 'FontWeight', 'bold');
xlim([-fs/2/1e6, fs/2/1e6]);

% 子图3: 循环谱在alpha方向的切片（固定频率f=0）
subplot(2, 3, 3);
f_zero_idx = round(Nfft/2) + 1;  % f=0对应的索引
cyclic_slice_f0 = cyclic_spectrum(f_zero_idx, :);
plot(alpha_axis/1e3, 10*log10(abs(cyclic_slice_f0) + eps), 'r-', 'LineWidth', 1.5);
grid on;
xlabel('循环频率 α (kHz)', 'FontSize', 11);
ylabel('幅度 (dB)', 'FontSize', 11);
title('频率 f = 0 处的循环谱切片', 'FontSize', 12, 'FontWeight', 'bold');
xlim([alpha_axis(1)/1e3, alpha_axis(end)/1e3]);

% 子图4: Alpha = 0 切片与Welch方法PSD对比
subplot(2, 3, 4);
plot(f_axis/1e6, PSD_alpha_zero_dB, 'b-', 'LineWidth', 1.5, 'DisplayName', '循环谱 α=0');
hold on;
plot(f_welch_shifted/1e6, Pxx_welch_dB_shifted, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Welch PSD');
grid on;
xlabel('频率 f (MHz)', 'FontSize', 11);
ylabel('功率谱密度 (dB)', 'FontSize', 11);
title('Alpha = 0 切片 vs Welch PSD', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
xlim([-fs/2/1e6, fs/2/1e6]);
hold off;

% 子图5: 循环谱的3D视图（仅显示主要部分）
subplot(2, 3, 5);
% 选择alpha范围的中心部分进行3D显示
alpha_center_idx = round(Nfft_alpha/2);
alpha_range = max(1, alpha_center_idx-50):min(Nfft_alpha, alpha_center_idx+50);
surf(alpha_axis(alpha_range)/1e3, f_axis/1e6, cyclic_spectrum_dB(:, alpha_range));
xlabel('循环频率 α (kHz)', 'FontSize', 11);
ylabel('频率 f (MHz)', 'FontSize', 11);
zlabel('幅度 (dB)', 'FontSize', 11);
title('循环谱3D视图（Alpha中心区域）', 'FontSize', 12, 'FontWeight', 'bold');
colormap('jet');
shading interp;
view(45, 30);

% 子图6: Alpha = 0 切片（线性尺度）
subplot(2, 3, 6);
plot(f_axis/1e6, abs(PSD_alpha_zero), 'g-', 'LineWidth', 1.5);
grid on;
xlabel('频率 f (MHz)', 'FontSize', 11);
ylabel('功率谱密度（线性尺度）', 'FontSize', 11);
title('Alpha = 0 切片（线性尺度）', 'FontSize', 12, 'FontWeight', 'bold');
xlim([-fs/2/1e6, fs/2/1e6]);

sgtitle('循环谱分析结果（Alpha = 0 切片提取）', 'FontSize', 14, 'FontWeight', 'bold');

fprintf('  图形绘制完成！\n\n');

%===============================================================================
% 【步骤7: 保存结果】
%===============================================================================
fprintf('步骤7: 保存分析结果...\n');

% 保存结果到MAT文件
save('cyclic_spectrum_results.mat', 'cyclic_spectrum', 'PSD_alpha_zero', ...
    'f_axis', 'alpha_axis', 'alpha_zero_idx', 'fs', 'Nfft', 'Nfft_alpha', ...
    'window_length', 'noverlap', 'nsegments', 'Pxx_welch', 'f_welch', ...
    'Pxx_welch_dB', 'f_welch_shifted', 'Pxx_welch_dB_shifted', '-v7.3');

fprintf('  - 结果已保存到: cyclic_spectrum_results.mat\n');
fprintf('  - 包含变量: cyclic_spectrum, PSD_alpha_zero, f_axis, alpha_axis, ...\n');
fprintf('  完成！\n\n');

%===============================================================================
% 【结果总结】
%===============================================================================
fprintf('========================================\n');
fprintf('循环谱分析完成！\n');
fprintf('========================================\n');
fprintf('分析结果：\n');
fprintf('  - 循环谱矩阵大小: [%d × %d] (频率 × 循环频率)\n', size(cyclic_spectrum));
fprintf('  - Alpha = 0 切片长度: %d\n', length(PSD_alpha_zero));
fprintf('  - 频率分辨率: %.4f kHz\n', fs/Nfft/1e3);
fprintf('  - 循环频率分辨率: %.4f Hz\n', fs/Nfft_alpha);
fprintf('  - Alpha = 0 对应的循环频率: %.6f Hz\n', alpha_axis(alpha_zero_idx));
fprintf('  - 频率范围: %.3f MHz 到 %.3f MHz\n', f_axis(1)/1e6, f_axis(end)/1e6);
fprintf('  - 循环频率范围: %.3f kHz 到 %.3f kHz\n', alpha_axis(1)/1e3, alpha_axis(end)/1e3);
fprintf('========================================\n\n');

fprintf('所有步骤完成！\n');
fprintf('提示：Alpha = 0 的切片即为信号的功率谱密度（PSD）\n');
