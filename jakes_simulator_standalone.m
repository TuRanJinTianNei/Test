%% ============================================================================
% Jakes模型Rayleigh衰落信道仿真器 - 独立运行版本
% ============================================================================
% 功能：实现经典的Jakes模型来模拟Rayleigh衰落信道
% 特点：完全独立运行，包含所有必要的函数和绘图功能
% 
% 作者：基于原始项目代码整合
% 日期：2024
% ============================================================================

function jakes_simulator_standalone()
     %% ========================================================================
     % 1. 仿真参数设置
     % ========================================================================
     clear;
     close all;
     clc;
     
     % 仿真时间参数
     t_sim = 15;          % 仿真持续时间（秒）
     delta_t = 1e-3;      % 时间步长（秒）
     N_step = t_sim/delta_t;  % 时间步数
     time = linspace(0, t_sim, N_step)';  % 时间向量
     
     % Jakes模型参数
     M = 50;              % 振荡器数量参数（N = 4*M+2，共202个振荡器）
     fd = 50;             % 多普勒频率（Hz）
     
     fprintf('========================================\n');
     fprintf('Jakes模型Rayleigh衰落信道仿真\n');
     fprintf('========================================\n');
     fprintf('仿真参数：\n');
     fprintf('  M = %d (N = %d个振荡器)\n', M, 4*M+2);
     fprintf('  多普勒频率 fd = %.1f Hz\n', fd);
     fprintf('  仿真时长 = %.1f 秒\n', t_sim);
     fprintf('  时间步数 = %d\n', N_step);
     fprintf('========================================\n\n');
     
     %% ========================================================================
     % 2. 生成Jakes模型信道
     % ========================================================================
     fprintf('正在生成Jakes模型信道...\n');
     [g, pdf, R] = gen_jakes_channel(M, fd, time);
     fprintf('信道生成完成！\n\n');
     
     %% ========================================================================
     % 3. 生成理想参考模型（用于对比）
     % ========================================================================
     fprintf('正在生成理想参考模型...\n');
     [pdf_id, R_id] = gen_ideal_reference(fd, delta_t, N_step);
     fprintf('参考模型生成完成！\n\n');
     
     %% ========================================================================
     % 4. 绘制所有验证图表
     % ========================================================================
     fprintf('正在绘制验证图表...\n');
     plot_jakes_results(g, pdf, R, pdf_id, R_id, time, fd, delta_t, t_sim, M);
     fprintf('所有图表绘制完成！\n');
     fprintf('========================================\n');
 end
 
 %% ============================================================================
 % Jakes模型信道生成函数
 % ============================================================================
 function [g, pdf, R] = gen_jakes_channel(M, fd, time)
     % 输入：
     %   M    - 振荡器数量参数（N = 4*M+2）
     %   fd   - 多普勒频率（Hz）
     %   time - 时间向量
     % 输出：
     %   g    - 复衰落过程
     %   pdf  - 概率密度函数（cell数组：pdf{1}幅度，pdf{2}相位）
     %   R    - 相关函数（cell数组，包含6种相关函数）
     
     N = 4*M + 2;  % 总振荡器数量
     
     % 初始振荡器（n=0）
     beta_0 = pi/4;
     a_0 = sqrt(2)*cos(beta_0);
     b_0 = sqrt(2)*sin(beta_0);
     w_0 = 2*pi*fd;
     phi = 0;
     
     % 初始化同相和正交分量
     gr = a_0*cos(w_0*time + phi);
     gi = b_0*cos(w_0*time + phi);
     
     % 叠加M个振荡器
     for n = 1:M
         beta_n = n*pi/M;
         a_n = 2*cos(beta_n);
         b_n = 2*sin(beta_n);
         w_n = w_0*cos((2*pi*n)/N);
         
         gr = gr + a_n*cos(w_n*time + phi);
         gi = gi + b_n*cos(w_n*time + phi);
     end
     
     % 归一化并生成复衰落过程
     g = 2/sqrt(N)*(gr + 1i*gi);
     
     % 计算概率密度函数
     F = abs(g);  % 衰落幅度
     [pdf_abs, ~] = ksdensity(F, linspace(0, 4, 100));  % 幅度PDF
     TH = angle(g);  % 相位
     [pdf_ph, ~] = ksdensity(TH, linspace(-pi, pi, 100));  % 相位PDF
     
     pdf = cell(1, 2);
     pdf{1} = pdf_abs;
     pdf{2} = pdf_ph;
     
     % 计算相关函数
     R = compute_correlations(g, gr, gi);
 end
 
 %% ============================================================================
 % 理想参考模型生成函数
 % ============================================================================
 function [pdf_id, R_id] = gen_ideal_reference(fd, delta_t, N_step)
     % 生成理想Rayleigh衰落信道的理论参考值
     % 用于与Jakes模型仿真结果对比
     
     % 计算归一化时间延迟
     lags = (0:N_step-1)*(fd*delta_t);
     lags = lags';
     
     % 理论相关函数（基于Bessel函数J₀）
     R_id = cell(1, 6);
     R_id{1} = besselj(0, 2*pi*lags);  % 同相分量自相关
     R_id{2} = besselj(0, 2*pi*lags);  % 正交分量自相关
     R_id{3} = zeros(N_step, 1);       % 互相关（应该为0）
     R_id{4} = zeros(N_step, 1);       % 互相关（应该为0）
     R_id{5} = 2*besselj(0, 2*pi*lags);  % 复衰落自相关
     R_id{6} = 4 + 4*besselj(0, 2*pi*lags).^2;  % 功率自相关
     
     % 理论概率密度函数
     r = linspace(0, 4, 500);
     pdf_abs_id = r.*exp(-r.^2/2);  % 理论Rayleigh分布（σ²=1）
     pdf_ph_id = (1/(2*pi))*ones(1, 500);  % 理论均匀分布
     
     pdf_id = cell(1, 2);
     pdf_id{1} = pdf_abs_id;
     pdf_id{2} = pdf_ph_id;
 end
 
 %% ============================================================================
 % 相关函数计算函数
 % ============================================================================
 function R = compute_correlations(g, gr, gi)
     % 计算6种相关函数
     % R{1} = Rgcgc (同相分量自相关)
     % R{2} = Rgsgs (正交分量自相关)
     % R{3} = Rgcgs (互相关)
     % R{4} = Rgsgc (互相关)
     % R{5} = Rgg (复衰落自相关)
     % R{6} = Rg2g2 (功率自相关)
     
     N_step = length(g);
     R = cell(6, 1);
     
     Rgcgc = xcorr(gr, gr, 'coeff');
     Rgcgc = Rgcgc(N_step:end);
     
     Rgsgs = xcorr(gi, gi, 'coeff');
     Rgsgs = Rgsgs(N_step:end);
     
     Rgcgs = xcorr(gr, gi, 'coeff');
     Rgcgs = Rgcgs(N_step:end);
     
     Rgsgc = xcorr(gi, gr, 'coeff');
     Rgsgc = Rgsgc(N_step:end);
     
     Rgg = 2*xcorr(g, g, 'coeff');
     Rgg = Rgg(N_step:end);
     
     Rg2g2 = 8*xcorr(abs(g).^2, abs(g).^2, 'coeff');
     Rg2g2 = Rg2g2(N_step:end);
     
     R{1} = Rgcgc;
     R{2} = Rgsgs;
     R{3} = Rgcgs;
     R{4} = Rgsgc;
     R{5} = Rgg;
     R{6} = Rg2g2;
 end
 
 %% ============================================================================
 % 绘图函数
 % ============================================================================
 function plot_jakes_results(g, pdf, R, pdf_id, R_id, time, fd, delta_t, t_sim, M)
     % 绘制所有Jakes模型的验证图表
     
     ln_wdt = 2;      % 线宽
     f_size = 20;     % 字体大小
     N_step = length(g);
     
     % 归一化时间延迟
     lags = (0:N_step-1)*(fd*delta_t);
     lags = lags';
     
     %% Figure 1: 衰落包络概率密度函数（PDF）
     figure('Name', 'Jakes - Fading Envelope', 'NumberTitle', 'on');
     set(gcf, 'Position', [200, 100, 800, 600]);
     grid on; hold on;
     text(0.02, 0.98, 'Figure 1', 'Units', 'normalized', 'FontSize', 18, ...
          'FontWeight', 'bold', 'VerticalAlignment', 'top', ...
          'BackgroundColor', [1 1 1], 'EdgeColor', [0 0 0], 'Margin', 3);
     
     x = linspace(0, 4, 500);
     x_pdf = linspace(0, 4, 100);
     pdf_interp = interp1(x_pdf, pdf{1}, x, 'linear', 'extrap');
     
     plot(x, pdf_id{1}, 'k', 'LineWidth', ln_wdt);
     plot(x, pdf_interp, 'r--', 'LineWidth', ln_wdt);
     xlabel('$r$', 'Interpreter', 'latex', 'FontSize', f_size);
     ylabel('$f_{|g|}(r)$', 'Interpreter', 'latex', 'FontSize', f_size);
     xlim([0, 4]);
     legend({'Reference (Ideal)', sprintf('Jakes ($M=%d$)', M)}, ...
            'Interpreter', 'latex', 'FontSize', f_size);
     title('衰落包络概率密度函数', 'FontSize', f_size);
     
     %% Figure 2: 相位概率密度函数（PDF）
     figure('Name', 'Jakes - Phase PDF', 'NumberTitle', 'on');
     set(gcf, 'Position', [200, 100, 800, 600]);
     grid on; hold on;
     text(0.02, 0.98, 'Figure 2', 'Units', 'normalized', 'FontSize', 18, ...
          'FontWeight', 'bold', 'VerticalAlignment', 'top', ...
          'BackgroundColor', [1 1 1], 'EdgeColor', [0 0 0], 'Margin', 3);
     
     x_ph_id = linspace(-1, 1, 500);
     x_ph = linspace(-1, 1, 100);
     
     plot(x_ph_id, pdf_id{2}, 'k', 'LineWidth', ln_wdt);
     plot(x_ph, pdf{2}, 'r--', 'LineWidth', ln_wdt);
     xlabel('$\theta$ ($\times \pi$)', 'Interpreter', 'latex', 'FontSize', f_size);
     ylabel('$f_{\theta}(\theta)$', 'Interpreter', 'latex', 'FontSize', f_size);
     xlim([-1, 1]);
     legend({'Reference (Ideal)', sprintf('Jakes ($M=%d$)', M)}, ...
            'Interpreter', 'latex', 'FontSize', f_size);
     title('相位概率密度函数', 'FontSize', f_size);
     
     %% Figure 3: 复衰落过程自相关函数
     figure('Name', 'Jakes - Second Order Statistic', 'NumberTitle', 'on');
     set(gcf, 'Position', [200, 100, 800, 600]);
     grid on; hold on;
     text(0.02, 0.98, 'Figure 3', 'Units', 'normalized', 'FontSize', 18, ...
          'FontWeight', 'bold', 'VerticalAlignment', 'top', ...
          'BackgroundColor', [1 1 1], 'EdgeColor', [0 0 0], 'Margin', 3);
     
     plot(lags, R_id{5}, 'k', 'LineWidth', ln_wdt);
     plot(lags, R{5}, 'r--', 'LineWidth', ln_wdt);
     xlabel('Normalized Time: $f_d \tau$', 'Interpreter', 'latex', 'FontSize', f_size);
     ylabel('$R_{g,g}(\tau)$', 'Interpreter', 'latex', 'FontSize', f_size);
     xlim([0, min(15, t_sim*fd)]);
     legend({'Reference (Ideal)', sprintf('Jakes ($M=%d$)', M)}, ...
            'Interpreter', 'latex', 'FontSize', f_size);
     title('复衰落过程自相关函数', 'FontSize', f_size);
     
     %% Figure 4: 同相分量自相关函数
     figure('Name', 'Jakes - Second Order Statistic', 'NumberTitle', 'on');
     set(gcf, 'Position', [200, 100, 800, 600]);
     grid on; hold on;
     text(0.02, 0.98, 'Figure 4', 'Units', 'normalized', 'FontSize', 18, ...
          'FontWeight', 'bold', 'VerticalAlignment', 'top', ...
          'BackgroundColor', [1 1 1], 'EdgeColor', [0 0 0], 'Margin', 3);
     
     plot(lags, R_id{1}, 'k', 'LineWidth', ln_wdt);
     plot(lags, R{1}, 'r--', 'LineWidth', ln_wdt);
     xlabel('Normalized Time: $f_d \tau$', 'Interpreter', 'latex', 'FontSize', f_size);
     ylabel('$R_{g_c,g_c}(\tau)$', 'Interpreter', 'latex', 'FontSize', f_size);
     xlim([0, min(15, t_sim*fd)]);
     legend({'Reference (Ideal)', sprintf('Jakes ($M=%d$)', M)}, ...
            'Interpreter', 'latex', 'FontSize', f_size);
     title('同相分量自相关函数', 'FontSize', f_size);
     
     %% Figure 5: 同相与正交分量互相关函数
     figure('Name', 'Jakes - Second Order Statistic', 'NumberTitle', 'on');
     set(gcf, 'Position', [200, 100, 800, 600]);
     grid on; hold on;
     text(0.02, 0.98, 'Figure 5', 'Units', 'normalized', 'FontSize', 18, ...
          'FontWeight', 'bold', 'VerticalAlignment', 'top', ...
          'BackgroundColor', [1 1 1], 'EdgeColor', [0 0 0], 'Margin', 3);
     
     plot(lags, R_id{3}, 'k', 'LineWidth', ln_wdt);
     plot(lags, R{3}, 'r--', 'LineWidth', ln_wdt);
     xlabel('Normalized Time: $f_d \tau$', 'Interpreter', 'latex', 'FontSize', f_size);
     ylabel('$R_{g_c,g_s}(\tau)$', 'Interpreter', 'latex', 'FontSize', f_size);
     ylim([-0.5, 0.5]);
     xlim([0, min(15, t_sim*fd)]);
     legend({'Reference (Ideal)', sprintf('Jakes ($M=%d$)', M)}, ...
            'Interpreter', 'latex', 'FontSize', f_size);
     title('同相与正交分量互相关函数', 'FontSize', f_size);
     
     %% Figure 6: 时域波形
     figure('Name', 'Jakes - Time Domain Waveform', 'NumberTitle', 'on');
     set(gcf, 'Position', [200, 100, 1000, 700]);
     
     % 幅度时域波形
     subplot(2, 1, 1);
     hold on;
     text(0.02, 0.98, 'Figure 6', 'Units', 'normalized', 'FontSize', 18, ...
          'FontWeight', 'bold', 'VerticalAlignment', 'top', ...
          'BackgroundColor', [1 1 1], 'EdgeColor', [0 0 0], 'Margin', 3);
     plot(time(1:min(5000, length(time))), abs(g(1:min(5000, length(g)))), ...
          'b-', 'LineWidth', 1);
     grid on;
     xlabel('时间 t (s)', 'FontSize', f_size);
     ylabel('幅度 |g(t)|', 'FontSize', f_size);
     title('Jakes信道时域波形 - 幅度', 'FontSize', f_size);
     xlim([0, min(5, t_sim)]);
     
     % 相位时域波形
     subplot(2, 1, 2);
     hold on;
     text(0.02, 0.98, 'Figure 6-cont', 'Units', 'normalized', 'FontSize', 18, ...
          'FontWeight', 'bold', 'VerticalAlignment', 'top', ...
          'BackgroundColor', [1 1 1], 'EdgeColor', [0 0 0], 'Margin', 3);
     plot(time(1:min(5000, length(time))), angle(g(1:min(5000, length(g))))/pi, ...
          'r-', 'LineWidth', 1);
     grid on;
     xlabel('时间 t (s)', 'FontSize', f_size);
     ylabel('相位 \theta(t)/\pi', 'FontSize', f_size);
     title('Jakes信道时域波形 - 相位', 'FontSize', f_size);
     ylim([-1, 1]);
     xlim([0, min(5, t_sim)]);
     
     %% Figure 7: 幅度直方图
     figure('Name', 'Jakes - Amplitude Histogram', 'NumberTitle', 'on');
     set(gcf, 'Position', [200, 100, 800, 600]);
     grid on; hold on;
     text(0.02, 0.98, 'Figure 7', 'Units', 'normalized', 'FontSize', 18, ...
          'FontWeight', 'bold', 'VerticalAlignment', 'top', ...
          'BackgroundColor', [1 1 1], 'EdgeColor', [0 0 0], 'Margin', 3);
     
     abs_g = abs(g);
     [counts, centers] = hist(abs_g, 50);
     bin_width = centers(2) - centers(1);
     pdf_hist = counts / (sum(counts) * bin_width);
     
     % 理论Rayleigh分布
     r = linspace(0, max(abs_g)*1.2, 500);
     sigma_r = sqrt(mean(abs_g.^2)/2);
     pdf_rayleigh = (r / sigma_r^2) .* exp(-r.^2 / (2*sigma_r^2));
     pdf_rayleigh(r < 0) = 0;
     
     bar(centers, pdf_hist, 'FaceColor', [0.7 0.7 0.9], 'EdgeColor', 'none', 'BarWidth', 1);
     plot(r, pdf_rayleigh, 'r-', 'LineWidth', ln_wdt);
     xlabel('幅度 r', 'FontSize', f_size);
     ylabel('概率密度 f_{|g|}(r)', 'FontSize', f_size);
     legend({'仿真直方图', '理论Rayleigh分布'}, 'FontSize', f_size-2);
     title('幅度分布验证 - 瑞利分布', 'FontSize', f_size);
     xlim([0, max(abs_g)*1.1]);
     
     %% Figure 8: 多普勒频谱（核心验证）
     figure('Name', 'Jakes - Doppler Spectrum', 'NumberTitle', 'on');
     set(gcf, 'Position', [200, 100, 800, 600]);
     grid on; hold on;
     text(0.02, 0.98, 'Figure 8', 'Units', 'normalized', 'FontSize', 18, ...
          'FontWeight', 'bold', 'VerticalAlignment', 'top', ...
          'BackgroundColor', [1 1 1], 'EdgeColor', [0 0 0], 'Margin', 3);
     
     % 计算功率谱密度
     N_fft = 2^nextpow2(length(g));
     G_fft = fft(g, N_fft);
     PSD = abs(G_fft).^2;
     PSD = fftshift(PSD);
     
     % 频率轴
     f = (-N_fft/2:N_fft/2-1) / (N_fft * delta_t);
     f_norm = f / fd;
     
     % 理论U型频谱
     f_theory = linspace(-0.99, 0.99, 1000);
     PSD_theory = 1 ./ (pi * sqrt(1 - f_theory.^2));
     PSD_theory = PSD_theory / max(PSD_theory);
     
     % 归一化仿真频谱
     PSD_norm = PSD / max(PSD);
     
    % 绘制频谱（只显示多普勒频率范围内）
    idx_plot = abs(f_norm) <= 1;
    plot(f_norm(idx_plot), 10*log10(PSD_norm(idx_plot) + eps), 'b-', 'LineWidth', ln_wdt);
    plot(f_theory, 10*log10(PSD_theory + eps), 'r--', 'LineWidth', ln_wdt);
    xlabel('归一化频率 f/f_d', 'FontSize', f_size);
    ylabel('归一化功率谱密度 (dB)', 'FontSize', f_size);
    xlim([-1, 1]);
    ylim([-40, 5]);
    legend({'仿真结果', '理论U型频谱'}, 'FontSize', f_size-2);
    title('多普勒功率谱密度 - U型频谱验证', 'FontSize', f_size);
     
     %% Figure 9: 自相关函数J₀验证（详细视图）
     figure('Name', 'Jakes - Autocorrelation J0 Verification', 'NumberTitle', 'on');
     set(gcf, 'Position', [200, 100, 800, 600]);
     grid on; hold on;
     text(0.02, 0.98, 'Figure 9', 'Units', 'normalized', 'FontSize', 18, ...
          'FontWeight', 'bold', 'VerticalAlignment', 'top', ...
          'BackgroundColor', [1 1 1], 'EdgeColor', [0 0 0], 'Margin', 3);
     
    idx_detail = lags <= 5;
    plot(lags(idx_detail), R_id{5}(idx_detail), 'k', 'LineWidth', ln_wdt);
    plot(lags(idx_detail), R{5}(idx_detail), 'r--', 'LineWidth', ln_wdt);
    xlabel('归一化时间 f_d \tau', 'FontSize', f_size);
    ylabel('自相关函数 R_{g,g}(\tau)', 'FontSize', f_size);
    xlim([0, 5]);
    legend({'理论 J_0(2\pi f_d \tau)', '仿真结果'}, ...
           'FontSize', f_size-2);
    title('自相关函数J₀验证（详细视图）', 'FontSize', f_size);
     
     %% Figure 10: 实部/虚部分布
     figure('Name', 'Jakes - Real/Imaginary Distribution', 'NumberTitle', 'on');
     set(gcf, 'Position', [200, 100, 1200, 600]);
     
     % 实部分布
     subplot(1, 2, 1);
     hold on;
     text(0.02, 0.98, 'Figure 10', 'Units', 'normalized', 'FontSize', 18, ...
          'FontWeight', 'bold', 'VerticalAlignment', 'top', ...
          'BackgroundColor', [1 1 1], 'EdgeColor', [0 0 0], 'Margin', 3);
     
     gc = real(g);
     [counts_gc, centers_gc] = hist(gc, 50);
     bin_width_gc = centers_gc(2) - centers_gc(1);
     pdf_hist_gc = counts_gc / (sum(counts_gc) * bin_width_gc);
     
     % 理论高斯分布
     x_gc = linspace(min(gc), max(gc), 500);
     mu_gc = mean(gc);
     sigma_gc = std(gc);
     pdf_gauss_gc = (1/(sigma_gc*sqrt(2*pi))) * exp(-(x_gc-mu_gc).^2/(2*sigma_gc^2));
     
     bar(centers_gc, pdf_hist_gc, 'FaceColor', [0.7 0.9 0.7], 'EdgeColor', 'none', 'BarWidth', 1);
     plot(x_gc, pdf_gauss_gc, 'r-', 'LineWidth', ln_wdt);
     xlabel('实部 $g_c$', 'Interpreter', 'latex', 'FontSize', f_size);
     ylabel('概率密度 $f_{g_c}(x)$', 'Interpreter', 'latex', 'FontSize', f_size);
     legend({'仿真直方图', '理论高斯分布'}, 'FontSize', f_size-2);
     title('实部分布验证', 'FontSize', f_size);
     grid on;
     
     % 虚部分布
     subplot(1, 2, 2);
     hold on;
     text(0.02, 0.98, 'Figure 10-cont', 'Units', 'normalized', 'FontSize', 18, ...
          'FontWeight', 'bold', 'VerticalAlignment', 'top', ...
          'BackgroundColor', [1 1 1], 'EdgeColor', [0 0 0], 'Margin', 3);
     
     gs = imag(g);
     [counts_gs, centers_gs] = hist(gs, 50);
     bin_width_gs = centers_gs(2) - centers_gs(1);
     pdf_hist_gs = counts_gs / (sum(counts_gs) * bin_width_gs);
     
     % 理论高斯分布
     x_gs = linspace(min(gs), max(gs), 500);
     mu_gs = mean(gs);
     sigma_gs = std(gs);
     pdf_gauss_gs = (1/(sigma_gs*sqrt(2*pi))) * exp(-(x_gs-mu_gs).^2/(2*sigma_gs^2));
     
     bar(centers_gs, pdf_hist_gs, 'FaceColor', [0.9 0.7 0.7], 'EdgeColor', 'none', 'BarWidth', 1);
     plot(x_gs, pdf_gauss_gs, 'r-', 'LineWidth', ln_wdt);
     xlabel('虚部 $g_s$', 'Interpreter', 'latex', 'FontSize', f_size);
     ylabel('概率密度 $f_{g_s}(x)$', 'Interpreter', 'latex', 'FontSize', f_size);
     legend({'仿真直方图', '理论高斯分布'}, 'FontSize', f_size-2);
     title('虚部分布验证', 'FontSize', f_size);
     grid on;
     
     fprintf('已生成10个验证图表！\n');
 end
 
 