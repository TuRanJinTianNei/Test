% ============================================================================
% 文件名: test6.m
% 功能: 二维高斯分布到瑞利分布的推导步骤 + 可视化示例
% 描述:
%   1. 展示二维高斯分布（两个独立的高斯随机变量X和Y）
%   2. 推导瑞利分布：R = sqrt(X^2 + Y^2) 的幅度分布
%   3. 推导相位分布：θ = atan2(Y, X) 的均匀分布
%   4. 生成多个Figure展示：
%      - Figure 1: 二维高斯分布的散点图（X-Y平面）
%      - Figure 2: 二维高斯分布的联合概率密度函数（3D曲面）
%      - Figure 3: 二维高斯分布的联合概率密度函数（等高线图）
%      - Figure 4: 幅度R的直方图与理论瑞利分布PDF对比
%      - Figure 5: 相位θ的直方图与理论均匀分布PDF对比
%      - Figure 6: 幅度R的累积分布函数（CDF）对比
%      - Figure 7: 极坐标表示（幅度-相位关系）
%      - Figure 8: 瑞利分布参数对PDF的影响
% 原理: 
%   设 X ~ N(0, σ²), Y ~ N(0, σ²)，且X和Y独立
%   定义 R = sqrt(X² + Y²) 为幅度，θ = atan2(Y, X) 为相位
%   则 R 服从瑞利分布：f_R(r) = (r/σ²) * exp(-r²/(2σ²)), r ≥ 0
%   且 θ 服从均匀分布：f_θ(θ) = 1/(2π), θ ∈ [0, 2π)
%   瑞利分布参数：σ（尺度参数）
% ============================================================================

clc; clear; close all;

%===============================================================================
% 【系统参数配置】
%===============================================================================
N_samples = 100000;              % 样本数量（足够大以保证统计准确性）
sigma = 1.0;                     % 高斯分布的标准差（瑞利分布的尺度参数）
% 注意：瑞利分布的PDF为 f_R(r) = (r/σ²) * exp(-r²/(2σ²))
%      其中σ是尺度参数，对应二维高斯分布的标准差

%===============================================================================
% 【步骤1: 生成二维高斯分布样本】
%===============================================================================
% 生成两个独立的标准高斯随机变量
% X ~ N(0, σ²), Y ~ N(0, σ²)，且X和Y独立
randn('state', 0);               % 设置随机种子以保证可重复性
X = sigma * randn(1, N_samples); % X分量：均值为0，标准差为σ
Y = sigma * randn(1, N_samples); % Y分量：均值为0，标准差为σ

%===============================================================================
% 【步骤2: 计算幅度和相位】
%===============================================================================
% 幅度：R = sqrt(X² + Y²)
R = sqrt(X.^2 + Y.^2);

% 相位：θ = atan2(Y, X)，范围 [-π, π]
theta = atan2(Y, X);
% 将相位转换到 [0, 2π] 范围以便与理论分布对比
theta_0_2pi = mod(theta, 2*pi);

%===============================================================================
% 【可视化分析】
%===============================================================================

%------------------------------------------------------------------------------
% Figure 1: 二维高斯分布的散点图（X-Y平面）
%------------------------------------------------------------------------------
figure(1);
scatter(X(1:5000), Y(1:5000), 5, 'filled', 'MarkerFaceAlpha', 0.3);
grid on; hold on;
xlabel('X (实部)');
ylabel('Y (虚部)');
title(sprintf('二维高斯分布散点图 (σ=%.2f, 显示前5000个样本)', sigma));
axis equal;
xlim([-4*sigma, 4*sigma]);
ylim([-4*sigma, 4*sigma]);

% 添加理论等高线（2D高斯分布的等高线是同心圆）
theta_circle = linspace(0, 2*pi, 100);
for r_mult = 1:3
    r_contour = r_mult * sigma;
    x_contour = r_contour * cos(theta_circle);
    y_contour = r_contour * sin(theta_circle);
    plot(x_contour, y_contour, 'r--', 'LineWidth', 1.5);
end
legend({'样本点', '理论等高线 (1σ, 2σ, 3σ)'}, 'Location', 'best');
hold off;

%------------------------------------------------------------------------------
% Figure 2: 二维高斯分布的联合概率密度函数（3D曲面）
%------------------------------------------------------------------------------
figure(2);
x_grid = linspace(-4*sigma, 4*sigma, 100);
y_grid = linspace(-4*sigma, 4*sigma, 100);
[X_grid, Y_grid] = meshgrid(x_grid, y_grid);

% 理论联合PDF：f(x,y) = (1/(2πσ²)) * exp(-(x²+y²)/(2σ²))
Z_theory = (1/(2*pi*sigma^2)) * exp(-(X_grid.^2 + Y_grid.^2)/(2*sigma^2));

surf(X_grid, Y_grid, Z_theory);
xlabel('X (实部)');
ylabel('Y (虚部)');
zlabel('联合概率密度 f(x,y)');
title(sprintf('二维高斯分布联合PDF (σ=%.2f)', sigma));
colorbar;
shading interp;
view(45, 30);

%------------------------------------------------------------------------------
% Figure 3: 二维高斯分布的联合概率密度函数（等高线图）
%------------------------------------------------------------------------------
figure(3);
contour(X_grid, Y_grid, Z_theory, 20);
xlabel('X (实部)');
ylabel('Y (虚部)');
title(sprintf('二维高斯分布联合PDF等高线图 (σ=%.2f)', sigma));
colorbar;
grid on;
axis equal;

%------------------------------------------------------------------------------
% Figure 4: 幅度R的直方图与理论瑞利分布PDF对比
%------------------------------------------------------------------------------
figure(4);
% 直方图
r_max = 5*sigma;
r_bins = linspace(0, r_max, 100);
hist_counts = histcounts(R, r_bins);
r_centers = (r_bins(1:end-1) + r_bins(2:end)) / 2;
hist_pdf = hist_counts / (N_samples * (r_bins(2) - r_bins(1))); % 归一化为PDF

% 理论瑞利分布PDF：f_R(r) = (r/σ²) * exp(-r²/(2σ²)), r ≥ 0
r_theory = linspace(0, r_max, 1000);
rayleigh_pdf_theory = (r_theory / sigma^2) .* exp(-r_theory.^2 / (2*sigma^2));

% 绘制对比
bar(r_centers, hist_pdf, 'FaceColor', [0.7 0.7 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
hold on;
plot(r_theory, rayleigh_pdf_theory, 'r-', 'LineWidth', 2);
grid on;
xlabel('幅度 R = \sqrt{X^2 + Y^2}');
ylabel('概率密度 f_R(r)');
title(sprintf('幅度R的分布：仿真直方图 vs 理论瑞利分布 (σ=%.2f)', sigma));
legend({'仿真直方图', '理论瑞利分布'}, 'Location', 'best');
xlim([0, r_max]);

% 计算统计量
R_mean_sim = mean(R);
R_var_sim = var(R);
R_mean_theory = sigma * sqrt(pi/2);  % 理论均值
R_var_theory = (4-pi)/2 * sigma^2;   % 理论方差
fprintf('\n==== 幅度R的统计量对比 ====\n');
fprintf('仿真均值: %.6f, 理论均值: %.6f (误差: %.2f%%)\n', ...
    R_mean_sim, R_mean_theory, 100*abs(R_mean_sim-R_mean_theory)/R_mean_theory);
fprintf('仿真方差: %.6f, 理论方差: %.6f (误差: %.2f%%)\n', ...
    R_var_sim, R_var_theory, 100*abs(R_var_sim-R_var_theory)/R_var_theory);
fprintf('================================\n');

%------------------------------------------------------------------------------
% Figure 5: 相位θ的直方图与理论均匀分布PDF对比
%------------------------------------------------------------------------------
figure(5);
% 直方图
theta_bins = linspace(0, 2*pi, 50);
hist_counts_theta = histcounts(theta_0_2pi, theta_bins);
theta_centers = (theta_bins(1:end-1) + theta_bins(2:end)) / 2;
hist_pdf_theta = hist_counts_theta / (N_samples * (theta_bins(2) - theta_bins(1))); % 归一化为PDF

% 理论均匀分布PDF：f_θ(θ) = 1/(2π), θ ∈ [0, 2π)
uniform_pdf_theory = ones(size(theta_centers)) / (2*pi);

% 绘制对比
bar(theta_centers, hist_pdf_theta, 'FaceColor', [0.9 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
hold on;
plot(theta_centers, uniform_pdf_theory, 'b-', 'LineWidth', 2);
grid on;
xlabel('相位 θ = atan2(Y, X) (弧度)');
ylabel('概率密度 f_θ(θ)');
title('相位θ的分布：仿真直方图 vs 理论均匀分布');
legend({'仿真直方图', '理论均匀分布 (1/(2π))'}, 'Location', 'best');
xlim([0, 2*pi]);
xticks(0:pi/2:2*pi);
xticklabels({'0', 'π/2', 'π', '3π/2', '2π'});

%------------------------------------------------------------------------------
% Figure 6: 幅度R的累积分布函数（CDF）对比
%------------------------------------------------------------------------------
figure(6);
% 仿真CDF（经验CDF）
r_sorted = sort(R);
cdf_sim = (1:N_samples) / N_samples;

% 理论瑞利分布CDF：F_R(r) = 1 - exp(-r²/(2σ²)), r ≥ 0
r_cdf = linspace(0, r_max, 1000);
cdf_theory = 1 - exp(-r_cdf.^2 / (2*sigma^2));

plot(r_sorted, cdf_sim, 'b-', 'LineWidth', 1.5, 'DisplayName', '仿真CDF');
hold on;
plot(r_cdf, cdf_theory, 'r--', 'LineWidth', 2, 'DisplayName', '理论瑞利CDF');
grid on;
xlabel('幅度 R');
ylabel('累积分布函数 F_R(r)');
title(sprintf('幅度R的累积分布函数对比 (σ=%.2f)', sigma));
legend('Location', 'best');
xlim([0, r_max]);

%------------------------------------------------------------------------------
% Figure 7: 极坐标表示（幅度-相位关系）
%------------------------------------------------------------------------------
figure(7);
% 子图1：极坐标散点图
subplot(2, 2, 1);
polarscatter(theta_0_2pi(1:5000), R(1:5000), 5, 'filled', 'MarkerFaceAlpha', 0.3);
title('极坐标表示：幅度-相位散点图（前5000个样本）');

% 子图2：幅度分布（极坐标下的径向分布）
subplot(2, 2, 2);
polarhistogram(R(1:5000), 30, 'FaceColor', [0.7 0.7 0.9], 'FaceAlpha', 0.6);
title('幅度R的极坐标直方图');

% 子图3：相位分布（极坐标下的角度分布）
subplot(2, 2, 3);
polarhistogram(theta_0_2pi(1:5000), 30, 'FaceColor', [0.9 0.7 0.7], 'FaceAlpha', 0.6);
title('相位θ的极坐标直方图');

% 子图4：幅度-相位联合分布（2D直方图）
subplot(2, 2, 4);
histogram2(theta_0_2pi, R, 30, 'DisplayStyle', 'tile', 'ShowEmptyBins', 'off');
xlabel('相位 θ (弧度)');
ylabel('幅度 R');
title('幅度-相位联合分布（2D直方图）');
colorbar;
xticks(0:pi/2:2*pi);
xticklabels({'0', 'π/2', 'π', '3π/2', '2π'});

%------------------------------------------------------------------------------
% Figure 8: 瑞利分布参数σ对PDF的影响
%------------------------------------------------------------------------------
figure(8);
sigma_values = [0.5, 1.0, 1.5, 2.0];  % 不同的σ值
colors_8 = lines(length(sigma_values));
r_plot = linspace(0, 6, 1000);

hold on;
for i = 1:length(sigma_values)
    sig = sigma_values(i);
    rayleigh_pdf = (r_plot / sig^2) .* exp(-r_plot.^2 / (2*sig^2));
    plot(r_plot, rayleigh_pdf, 'Color', colors_8(i,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('σ = %.1f', sig));
end
grid on;
xlabel('幅度 r');
ylabel('概率密度 f_R(r)');
title('瑞利分布PDF随参数σ的变化');
legend('Location', 'best');
xlim([0, 6]);

% 标注峰值位置（瑞利分布的众数在 r = σ 处）
for i = 1:length(sigma_values)
    sig = sigma_values(i);
    peak_value = (sig / sig^2) * exp(-sig^2 / (2*sig^2));
    plot(sig, peak_value, 'o', 'Color', colors_8(i,:), 'MarkerSize', 8, 'MarkerFaceColor', colors_8(i,:));
end

%===============================================================================
% 【命令行输出摘要】
%===============================================================================
fprintf('\n==== 二维高斯分布到瑞利分布推导总结 ====\n');
fprintf('样本数量: %d\n', N_samples);
fprintf('高斯分布标准差 σ: %.2f\n', sigma);
fprintf('\n【推导步骤】\n');
fprintf('1. 生成两个独立的高斯随机变量：X ~ N(0,σ²), Y ~ N(0,σ²)\n');
fprintf('2. 计算幅度：R = sqrt(X² + Y²)\n');
fprintf('3. 计算相位：θ = atan2(Y, X)\n');
fprintf('4. 结果：R 服从瑞利分布，θ 服从均匀分布\n');
fprintf('\n【理论结果】\n');
fprintf('幅度R的PDF: f_R(r) = (r/σ²) * exp(-r²/(2σ²)), r ≥ 0\n');
fprintf('幅度R的均值: E[R] = σ * sqrt(π/2) ≈ %.6f\n', sigma*sqrt(pi/2));
fprintf('幅度R的方差: Var[R] = (4-π)/2 * σ² ≈ %.6f\n', (4-pi)/2*sigma^2);
fprintf('相位θ的PDF: f_θ(θ) = 1/(2π), θ ∈ [0, 2π)\n');
fprintf('相位θ的均值: E[θ] = π\n');
fprintf('相位θ的方差: Var[θ] = π²/3\n');
fprintf('\n【应用场景】\n');
fprintf('瑞利分布常用于描述：\n');
fprintf('  - 无线通信中的瑞利衰落信道（多径传播）\n');
fprintf('  - 复高斯随机过程的幅度分布\n');
fprintf('  - 二维随机游走的距离分布\n');
fprintf('==========================================\n\n');

