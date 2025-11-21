%% ============================================================================
% Jakes模型Rayleigh衰落信道仿真器 - 独立运行版本
% ============================================================================
% 功能：实现经典的Jakes模型来模拟Rayleigh衰落信道
% 特点：完全独立运行，包含所有必要的函数和绘图功能
% 
% 作者：基于原始项目代码整合
% 日期：2024
% ============================================================================
% 
% Jakes模型数学原理：
%   - 通过多个正弦波叠加模拟多径衰落
%   - 总路径数：N = 4*M + 2（M为基础参数）
%   - 复信道：g(t) = g_c(t) + j*g_s(t)
%   - 包络服从瑞利分布，相位均匀分布
%   - 自相关函数：R(τ) = J_0(2πf_d*τ)
%   - 功率谱：U型多普勒功率谱
% ============================================================================
% 
% 仿真参数说明：
% ============================================================================
% 
% 【时间参数】
%   t_sim     仿真持续时间（秒）
%              默认值：15秒
%              说明：足够长的仿真时间以确保统计特性收敛
% 
%   delta_t   时间步长（秒）
%              默认值：1e-3秒（1毫秒）
%              说明：满足奈奎斯特准则，采样频率 fs = 1/delta_t = 1000 Hz
% 
%   N_step    时间步数
%              计算公式：N_step = t_sim / delta_t
%              默认值：15000
% 
% 【Jakes模型参数】
%   M         基础振荡器数量参数
%              默认值：50
%              说明：总振荡器数量 N = 4*M + 2 = 202
%                    - M个正频率路径（f > 0）
%                    - M个负频率路径（f < 0，镜像）
%                    - 2*M个特殊路径（确保I/Q独立性）
%                    - 2个特殊路径（DC和最高频）
%              注意：实际实现中通常只用M个路径即可达到相同统计特性
% 
%   v         移动速度（m/s）
%              默认值：7000 m/s（7 km/s）
%              说明：固定移动速度，用于计算载波频率
%              应用场景：卫星通信、高速移动平台等
% 
%   fd        最大多普勒频移（Hz）
%              默认值：50 Hz
%              计算公式：fd = v * fc / c
%              其中：
%                v  = 移动速度（m/s）
%                fc = 载波频率（Hz）
%                c  = 光速（3×10^8 m/s）
%              当前设置：
%                v = 7 km/s = 7000 m/s
%                fd = 50 Hz
%                fc = fd × c / v = 50 × 3×10^8 / 7000 ≈ 2.143 GHz
% 
% 【绘图参数】
%   ln_wdt    线宽（Line Width）
%              默认值：2
% 
%   f_size    字体大小（Font Size）
%              默认值：20
% 
% 【输出说明】
%   本程序生成10个验证图表：
%     Figure 1:  衰落包络PDF（瑞利分布验证）
%     Figure 2:  相位PDF（均匀分布验证）
%     Figure 3:  复衰落自相关函数（J₀函数验证）
%     Figure 4:  同相分量自相关函数
%     Figure 5:  I/Q互相关函数（正交性验证）
%     Figure 6:  时域波形（幅度和相位）
%     Figure 7:  幅度直方图（瑞利分布）
%     Figure 8:  多普勒功率谱（U型谱验证）
%     Figure 9:  自相关函数J₀验证（详细视图）
%     Figure 10: 实部/虚部分布（高斯分布验证）
% 
% ============================================================================

function jakes_simulator_standalone()
    %% ========================================================================
    % 1. 仿真参数设置
    % ========================================================================
    clear;      % 清除工作空间变量
    close all;  % 关闭所有图形窗口
    clc;        % 清空命令窗口
    
    %--------------------------------------------------------------------------
    % 仿真时间参数
    %--------------------------------------------------------------------------
    t_sim = 15;          % 仿真持续时间（秒）- 足够长以获取统计特性
    delta_t = 1e-3;      % 时间步长（秒）- 1毫秒，满足奈奎斯特准则
    N_step = t_sim/delta_t;  % 时间步数 = 15000
    time = linspace(0, t_sim, N_step)';  % 时间向量（列向量）
    
    %--------------------------------------------------------------------------
    % Jakes模型参数
    %--------------------------------------------------------------------------
    M = 50;              % 基础振荡器数量参数
    % 总振荡器数量：N = 4*M + 2（Jakes标准公式）
    % - 第一个M：正频率路径（f > 0）
    % - 第二个M：负频率路径（f < 0，镜像）
    % - 2*M：特殊路径（确保I/Q独立性）
    % - +2：两个特殊路径（DC和最高频）
    % 实际实现中通常只用M个路径即可达到相同统计特性
    
    %--------------------------------------------------------------------------
    % 物理参数
    %--------------------------------------------------------------------------
    v = 7000;            % 移动速度（m/s）= 7 km/s
    % 应用场景：卫星通信、高速移动平台等
    % 7 km/s ≈ 25200 km/h（低地球轨道卫星典型速度）
    
    c = 3e8;             % 光速（m/s）
    
    fd = 50;             % 最大多普勒频移（Hz）
    % 计算公式：fd = v * fc / c
    % 根据给定的v和fd，计算对应的载波频率fc
    
    fc = fd * c / v;     % 载波频率（Hz）
    % fc = 50 × 3×10^8 / 7000 ≈ 2.143 GHz
    
    fprintf('========================================\n');
    fprintf('Jakes模型Rayleigh衰落信道仿真\n');
    fprintf('========================================\n');
    fprintf('仿真参数：\n');
    fprintf('  M = %d (N = %d个振荡器)\n', M, 4*M+2);
    fprintf('  移动速度 v = %.1f m/s = %.1f km/s\n', v, v/1000);
    fprintf('  载波频率 fc = %.3f GHz\n', fc/1e9);
    fprintf('  多普勒频率 fd = %.1f Hz\n', fd);
    fprintf('  仿真时长 = %.1f 秒\n', t_sim);
    fprintf('  时间步数 = %d\n', N_step);
    fprintf('========================================\n\n');
    
        %% ========================================================================
    % 2. 生成Jakes模型信道
    % ========================================================================
    % 调用gen_jakes_channel函数生成Jakes信道
    % 输出：
    %   g   - 复衰落过程 g(t) = g_c(t) + j*g_s(t)
    %   pdf - 概率密度函数（包络和相位）
    %   R   - 相关函数（6种不同类型的相关函数）
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
% 功能：根据Jakes模型生成Rayleigh衰落信道
% 
% 数学原理：
%   g(t) = sqrt(2/N) * [g_c(t) + j*g_s(t)]
%   其中：
%     g_c(t) = a_0*cos(w_0*t) + sum_{n=1}^M a_n*cos(w_n*t)
%     g_s(t) = b_0*cos(w_0*t) + sum_{n=1}^M b_n*cos(w_n*t)
%     w_n = 2π*fd*cos(2π*n/N)  - 多普勒频移
%     a_n, b_n - 幅度系数，由角度beta_n决定
% ============================================================================
function [g, pdf, R] = gen_jakes_channel(M, fd, time)
    % 输入参数：
    %   M    - 基础振荡器数量参数（总路径数 N = 4*M+2）
    %   fd   - 最大多普勒频移（Hz）
    %   time - 时间向量（列向量，单位：秒）
    % 
    % 输出参数：
    %   g    - 复衰落过程 g(t) = g_c(t) + j*g_s(t)
    %   pdf  - 概率密度函数（cell数组）
    %          pdf{1}: 包络幅度PDF
    %          pdf{2}: 相位PDF
    %   R    - 相关函数（cell数组，包含6种相关函数）
    %          R{1}: 同相分量自相关 R_{g_c,g_c}
    %          R{2}: 正交分量自相关 R_{g_s,g_s}
    %          R{3}: 互相关 R_{g_c,g_s}
    %          R{4}: 互相关 R_{g_s,g_c}
    %          R{5}: 复衰落自相关 R_{g,g}
    %          R{6}: 功率自相关 R_{|g|^2,|g|^2}
    
    N = 4*M + 2;  % 总振荡器数量（Jakes标准公式）
    
    %--------------------------------------------------------------------------
    % 步骤1: 初始化特殊振荡器（n=0）
    %--------------------------------------------------------------------------
    % n=0对应最高多普勒频移路径（f = ±fd）
    beta_0 = pi/4;              % 特殊角度（π/4）
    a_0 = sqrt(2)*cos(beta_0);  % I分量幅度系数
    b_0 = sqrt(2)*sin(beta_0);  % Q分量幅度系数
    w_0 = 2*pi*fd;              % 最大多普勒角频率（rad/s）
    phi = 0;                     % 初始相位（可设为随机值）
    
    % 初始化同相分量（I分量）和正交分量（Q分量）
    gr = a_0*cos(w_0*time + phi);  % g_c(t)的初始值
    gi = b_0*cos(w_0*time + phi);  % g_s(t)的初始值
    
    %--------------------------------------------------------------------------
    % 步骤2: 叠加M个振荡器（n=1,2,...,M）
    %--------------------------------------------------------------------------
    % 每个振荡器对应不同的多普勒频移，模拟不同到达角度的路径
    for n = 1:M
        % 计算角度：beta_n = n*π/M（等间隔分布）
        beta_n = n*pi/M;
        
        % 计算幅度系数（确保功率归一化）
        a_n = 2*cos(beta_n);  % I分量幅度系数
        b_n = 2*sin(beta_n);  % Q分量幅度系数
        
        % 计算多普勒角频率：w_n = w_0*cos(2π*n/N)
        % 这对应到达角度α_n，使得f_n = fd*cos(α_n)
        w_n = w_0*cos((2*pi*n)/N);
        
        % 叠加到I和Q分量
        gr = gr + a_n*cos(w_n*time + phi);
        gi = gi + b_n*cos(w_n*time + phi);
    end
    
    %--------------------------------------------------------------------------
    % 步骤3: 归一化并生成复衰落过程
    %--------------------------------------------------------------------------
    % 归一化系数：2/sqrt(N)，确保平均功率为1
    % 复衰落过程：g(t) = g_c(t) + j*g_s(t)
    g = 2/sqrt(N)*(gr + 1i*gi);
    
    %--------------------------------------------------------------------------
    % 步骤4: 计算概率密度函数（PDF）
    %--------------------------------------------------------------------------
    % 使用核密度估计（KDE）方法估计PDF，比直方图更平滑
    F = abs(g);  % 衰落包络幅度 |g(t)|
    [pdf_abs, ~] = ksdensity(F, linspace(0, 4, 100));  % 幅度PDF（瑞利分布）
    
    TH = angle(g);  % 相位 arg{g(t)}
    [pdf_ph, ~] = ksdensity(TH, linspace(-pi, pi, 100));  % 相位PDF（均匀分布）
    
    % 存储PDF结果
    pdf = cell(1, 2);
    pdf{1} = pdf_abs;  % 包络幅度PDF
    pdf{2} = pdf_ph;   % 相位PDF
    
    %--------------------------------------------------------------------------
    % 步骤5: 计算相关函数
    %--------------------------------------------------------------------------
    % 计算6种不同类型的相关函数，用于验证Jakes模型的统计特性
    R = compute_correlations(g, gr, gi);
end

%% ============================================================================
% 理想参考模型生成函数
% ============================================================================
% 功能：生成理想Rayleigh衰落信道的理论参考值
% 
% 数学原理：
%   - 包络PDF：f_R(r) = (r/σ²)*exp(-r²/(2σ²))，瑞利分布
%   - 相位PDF：f_Θ(θ) = 1/(2π)，均匀分布
%   - 自相关函数：R(τ) = J_0(2πf_d*τ)，零阶贝塞尔函数
% ============================================================================
function [pdf_id, R_id] = gen_ideal_reference(fd, delta_t, N_step)
    % 输入参数：
    %   fd     - 最大多普勒频移（Hz）
    %   delta_t - 时间步长（秒）
    %   N_step  - 时间步数
    % 
    % 输出参数：
    %   pdf_id - 理论PDF（cell数组：pdf_id{1}幅度，pdf_id{2}相位）
    %   R_id   - 理论相关函数（cell数组，6种相关函数）
    
    %--------------------------------------------------------------------------
    % 步骤1: 计算归一化时间延迟
    %--------------------------------------------------------------------------
    % 归一化时间：τ_norm = fd*τ（无量纲）
    % 用于理论相关函数的计算
    lags = (0:N_step-1)*(fd*delta_t);  % 归一化时间延迟向量
    lags = lags';  % 转为列向量
    
    %--------------------------------------------------------------------------
    % 步骤2: 计算理论相关函数（基于零阶贝塞尔函数J₀）
    %--------------------------------------------------------------------------
    % Jakes模型的理论自相关函数：R(τ) = J_0(2πf_d*τ)
    R_id = cell(1, 6);
    
    % R{1}: 同相分量自相关 R_{g_c,g_c}(τ) = J_0(2πf_d*τ)
    R_id{1} = besselj(0, 2*pi*lags);
    
    % R{2}: 正交分量自相关 R_{g_s,g_s}(τ) = J_0(2πf_d*τ)
    R_id{2} = besselj(0, 2*pi*lags);
    
    % R{3}: 互相关 R_{g_c,g_s}(τ) = 0（I和Q分量应该正交）
    R_id{3} = zeros(N_step, 1);
    
    % R{4}: 互相关 R_{g_s,g_c}(τ) = 0（I和Q分量应该正交）
    R_id{4} = zeros(N_step, 1);
    
    % R{5}: 复衰落自相关 R_{g,g}(τ) = 2*J_0(2πf_d*τ)
    % 因为 g = g_c + j*g_s，所以 R_{g,g} = R_{g_c,g_c} + R_{g_s,g_s} = 2*J_0
    R_id{5} = 2*besselj(0, 2*pi*lags);
    
    % R{6}: 功率自相关 R_{|g|^2,|g|^2}(τ)
    % 理论值：R_{|g|^2,|g|^2}(τ) = 4 + 4*J_0²(2πf_d*τ)
    R_id{6} = 4 + 4*besselj(0, 2*pi*lags).^2;
    
    %--------------------------------------------------------------------------
    % 步骤3: 计算理论概率密度函数
    %--------------------------------------------------------------------------
    r = linspace(0, 4, 500);  % 幅度范围
    
    % 理论Rayleigh分布PDF：f_R(r) = (r/σ²)*exp(-r²/(2σ²))
    % 假设σ²=1（E[|g|²] = 2σ² = 2，所以σ²=1
    pdf_abs_id = r.*exp(-r.^2/2);
    
    % 理论均匀分布PDF：f_Θ(θ) = 1/(2π)，θ ∈ [-π, π]
    pdf_ph_id = (1/(2*pi))*ones(1, 500);
    
    % 存储PDF结果
    pdf_id = cell(1, 2);
    pdf_id{1} = pdf_abs_id;  % 包络幅度PDF（瑞利分布）
    pdf_id{2} = pdf_ph_id;   % 相位PDF（均匀分布）
end

%% ============================================================================
% 相关函数计算函数
% ============================================================================
% 功能：计算Jakes信道的6种相关函数，用于验证统计特性
% 
% 相关函数定义：
%   R_{x,y}(τ) = E[x(t)*y*(t+τ)]  （互相关）
%   R_{x,x}(τ) = E[x(t)*x*(t+τ)]  （自相关）
% 
% 归一化相关函数（'coeff'选项）：
%   R_norm(τ) = R(τ) / R(0)
% ============================================================================
function R = compute_correlations(g, gr, gi)
    % 输入参数：
    %   g  - 复衰落过程 g(t) = g_c(t) + j*g_s(t)
    %   gr - 同相分量（实部）g_c(t) = real(g)
    %   gi - 正交分量（虚部）g_s(t) = imag(g)
    % 
    % 输出参数：
    %   R  - 相关函数cell数组（6种）
    %        R{1}: R_{g_c,g_c}(τ) - 同相分量自相关
    %        R{2}: R_{g_s,g_s}(τ) - 正交分量自相关
    %        R{3}: R_{g_c,g_s}(τ) - 互相关（应该≈0）
    %        R{4}: R_{g_s,g_c}(τ) - 互相关（应该≈0）
    %        R{5}: R_{g,g}(τ)     - 复衰落自相关
    %        R{6}: R_{|g|^2,|g|^2}(τ) - 功率自相关
    
    N_step = length(g);
    R = cell(6, 1);
    
    %--------------------------------------------------------------------------
    % R{1}: 同相分量自相关 R_{g_c,g_c}(τ)
    %--------------------------------------------------------------------------
    % 理论值：R_{g_c,g_c}(τ) = J_0(2πf_d*τ)
    Rgcgc = xcorr(gr, gr, 'coeff');  % 'coeff'表示归一化到R(0)=1
    Rgcgc = Rgcgc(N_step:end);       % 只取τ≥0的部分（单边相关函数）
    
    %--------------------------------------------------------------------------
    % R{2}: 正交分量自相关 R_{g_s,g_s}(τ)
    %--------------------------------------------------------------------------
    % 理论值：R_{g_s,g_s}(τ) = J_0(2πf_d*τ)
    Rgsgs = xcorr(gi, gi, 'coeff');
    Rgsgs = Rgsgs(N_step:end);
    
    %--------------------------------------------------------------------------
    % R{3}: 互相关 R_{g_c,g_s}(τ)
    %--------------------------------------------------------------------------
    % 理论值：R_{g_c,g_s}(τ) = 0（I和Q分量应该正交）
    Rgcgs = xcorr(gr, gi, 'coeff');
    Rgcgs = Rgcgs(N_step:end);
    
    %--------------------------------------------------------------------------
    % R{4}: 互相关 R_{g_s,g_c}(τ)
    %--------------------------------------------------------------------------
    % 理论值：R_{g_s,g_c}(τ) = 0（I和Q分量应该正交）
    Rgsgc = xcorr(gi, gr, 'coeff');
    Rgsgc = Rgsgc(N_step:end);
    
    %--------------------------------------------------------------------------
    % R{5}: 复衰落自相关 R_{g,g}(τ)
    %--------------------------------------------------------------------------
    % 理论值：R_{g,g}(τ) = 2*J_0(2πf_d*τ)
    % 因为 g = g_c + j*g_s，所以 R_{g,g} = R_{g_c,g_c} + R_{g_s,g_s}
    Rgg = 2*xcorr(g, g, 'coeff');  % 乘以2是为了匹配理论值
    Rgg = Rgg(N_step:end);
    
    %--------------------------------------------------------------------------
    % R{6}: 功率自相关 R_{|g|^2,|g|^2}(τ)
    %--------------------------------------------------------------------------
    % 理论值：R_{|g|^2,|g|^2}(τ) = 4 + 4*J_0²(2πf_d*τ)
    Rg2g2 = 8*xcorr(abs(g).^2, abs(g).^2, 'coeff');  % 乘以8是为了匹配理论值
    Rg2g2 = Rg2g2(N_step:end);
    
    % 存储所有相关函数
    R{1} = Rgcgc;  % 同相分量自相关
    R{2} = Rgsgs;  % 正交分量自相关
    R{3} = Rgcgs;  % 互相关（I-Q）
    R{4} = Rgsgc;  % 互相关（Q-I）
    R{5} = Rgg;    % 复衰落自相关
    R{6} = Rg2g2;  % 功率自相关
end

%% ============================================================================
% 绘图函数
% ============================================================================
% 功能：绘制Jakes模型的所有验证图表，包括PDF、相关函数、时域波形等
% 
% 生成10个Figure：
%   Figure 1: 衰落包络PDF（瑞利分布验证）
%   Figure 2: 相位PDF（均匀分布验证）
%   Figure 3: 复衰落自相关函数（J₀函数验证）
%   Figure 4: 同相分量自相关函数
%   Figure 5: I/Q互相关函数（正交性验证）
%   Figure 6:
%   Figure 7: 幅度直方图（瑞利分布）
%   Figure 8: 多普勒功率谱（U型谱验证）
%   Figure 9: 自相关函数J₀验证（详细视图）
%   Figure 10: 实部/虚部分布（高斯分布验证）
% ============================================================================
function plot_jakes_results(g, pdf, R, pdf_id, R_id, time, fd, delta_t, t_sim, M)
    % 输入参数：
    %   g      - 复衰落过程
    %   pdf    - 仿真PDF（包络和相位）
    %   R      - 仿真相关函数（6种）
    %   pdf_id - 理论PDF（包络和相位）
    %   R_id   - 理论相关函数（6种）
    %   time   - 时间向量
    %   fd     - 最大多普勒频移（Hz）
    %   delta_t - 时间步长（秒）
    %   t_sim  - 仿真时长（秒）
    %   M      - 基础振荡器数量参数
    
    %--------------------------------------------------------------------------
    % 绘图参数设置
    %--------------------------------------------------------------------------
    ln_wdt = 2;      % 线宽（Line Width）
    f_size = 20;     % 字体大小（Font Size）
    N_step = length(g);
    
    %--------------------------------------------------------------------------
    % 计算归一化时间延迟
    %--------------------------------------------------------------------------
    % 归一化时间：τ_norm = fd*τ（无量纲）
    % 用于相关函数的x轴显示
    lags = (0:N_step-1)*(fd*delta_t);  % 归一化时间延迟向量
    lags = lags';  % 转为列向量
    
    %% Figure 1: 衰落包络概率密度函数（PDF）
    %--------------------------------------------------------------------------
    % 验证：包络幅度是否服从瑞利分布
    % 理论PDF：f_R(r) = (r/σ²)*exp(-r²/(2σ²))，r ≥ 0
    %--------------------------------------------------------------------------
    figure('Name', 'Jakes - Fading Envelope', 'NumberTitle', 'on');
    set(gcf, 'Position', [200, 100, 800, 600]);
    grid on; hold on;
    text(0.02, 0.98, 'Figure 1', 'Units', 'normalized', 'FontSize', 18, ...
         'FontWeight', 'bold', 'VerticalAlignment', 'top', ...
         'BackgroundColor', [1 1 1], 'EdgeColor', [0 0 0], 'Margin', 3);
    
    % 插值仿真PDF到理论PDF的采样点，便于对比
    x = linspace(0, 4, 500);           % 理论PDF的x轴采样点（幅度范围）
    x_pdf = linspace(0, 4, 100);      % 仿真PDF的x轴采样点（幅度范围）
    pdf_interp = interp1(x_pdf, pdf{1}, x, 'linear', 'extrap');  % 线性插值
    
    % 绘制理论曲线（黑色实线）和仿真曲线（红色虚线）
    plot(x, pdf_id{1}, 'k', 'LineWidth', ln_wdt);        % 理论Rayleigh分布
    plot(x, pdf_interp, 'r--', 'LineWidth', ln_wdt);     % 仿真结果
    xlabel('$r$', 'Interpreter', 'latex', 'FontSize', f_size);
    ylabel('$f_{|g|}(r)$', 'Interpreter', 'latex', 'FontSize', f_size);
    xlim([0, 4]);
    legend({'Reference (Ideal)', sprintf('Jakes ($M=%d$)', M)}, ...
           'Interpreter', 'latex', 'FontSize', f_size);
    title('衰落包络概率密度函数', 'FontSize', f_size);
    
    %% Figure 2: 相位概率密度函数（PDF）
    %--------------------------------------------------------------------------
    % 验证：相位是否均匀分布
    % 理论PDF：f_Θ(θ) = 1/(2π)，θ ∈ [-π, π]
    %--------------------------------------------------------------------------
    figure('Name', 'Jakes - Phase PDF', 'NumberTitle', 'on');
    set(gcf, 'Position', [200, 100, 800, 600]);
    grid on; hold on;
    text(0.02, 0.98, 'Figure 2', 'Units', 'normalized', 'FontSize', 18, ...
         'FontWeight', 'bold', 'VerticalAlignment', 'top', ...
         'BackgroundColor', [1 1 1], 'EdgeColor', [0 0 0], 'Margin', 3);
    
    % 相位归一化到[-1, 1]（单位：π）
    x_ph_id = linspace(-1, 1, 500);   % 理论PDF的x轴（归一化相位）
    x_ph = linspace(-1, 1, 100);      % 仿真PDF的x轴（归一化相位）
    
    % 绘制理论曲线（黑色实线）和仿真曲线（红色虚线）
    plot(x_ph_id, pdf_id{2}, 'k', 'LineWidth', ln_wdt);  % 理论均匀分布
    plot(x_ph, pdf{2}, 'r--', 'LineWidth', ln_wdt);       % 仿真结果
    xlabel('$\theta$ ($\times \pi$)', 'Interpreter', 'latex', 'FontSize', f_size);
    ylabel('$f_{\theta}(\theta)$', 'Interpreter', 'latex', 'FontSize', f_size);
    xlim([-1, 1]);
    legend({'Reference (Ideal)', sprintf('Jakes ($M=%d$)', M)}, ...
           'Interpreter', 'latex', 'FontSize', f_size);
    title('相位概率密度函数', 'FontSize', f_size);
    
    %% Figure 3: 复衰落过程自相关函数
    %--------------------------------------------------------------------------
    % 验证：自相关函数是否等于 J₀(2πf_d*τ)
    % 理论值：R_{g,g}(τ) = 2*J_0(2πf_d*τ)
    %--------------------------------------------------------------------------
    figure('Name', 'Jakes - Second Order Statistic', 'NumberTitle', 'on');
    set(gcf, 'Position', [200, 100, 800, 600]);
    grid on; hold on;
    text(0.02, 0.98, 'Figure 3', 'Units', 'normalized', 'FontSize', 18, ...
         'FontWeight', 'bold', 'VerticalAlignment', 'top', ...
         'BackgroundColor', [1 1 1], 'EdgeColor', [0 0 0], 'Margin', 3);
    
    % 绘制理论曲线（黑色实线）和仿真曲线（红色虚线）
    plot(lags, R_id{5}, 'k', 'LineWidth', ln_wdt);   % 理论：2*J_0(2πf_d*τ)
    plot(lags, R{5}, 'r--', 'LineWidth', ln_wdt);    % 仿真结果
    xlabel('Normalized Time: $f_d \tau$', 'Interpreter', 'latex', 'FontSize', f_size);
    ylabel('$R_{g,g}(\tau)$', 'Interpreter', 'latex', 'FontSize', f_size);
    xlim([0, min(15, t_sim*fd)]);
    legend({'Reference (Ideal)', sprintf('Jakes ($M=%d$)', M)}, ...
           'Interpreter', 'latex', 'FontSize', f_size);
    title('复衰落过程自相关函数', 'FontSize', f_size);
    
    %% Figure 4: 同相分量自相关函数
    %--------------------------------------------------------------------------
    % 验证：同相分量自相关函数是否等于 J₀(2πf_d*τ)
    % 理论值：R_{g_c,g_c}(τ) = J_0(2πf_d*τ)
    %--------------------------------------------------------------------------
    figure('Name', 'Jakes - Second Order Statistic', 'NumberTitle', 'on');
    set(gcf, 'Position', [200, 100, 800, 600]);
    grid on; hold on;
    text(0.02, 0.98, 'Figure 4', 'Units', 'normalized', 'FontSize', 18, ...
         'FontWeight', 'bold', 'VerticalAlignment', 'top', ...
         'BackgroundColor', [1 1 1], 'EdgeColor', [0 0 0], 'Margin', 3);
    
    % 绘制理论曲线（黑色实线）和仿真曲线（红色虚线）
    plot(lags, R_id{1}, 'k', 'LineWidth', ln_wdt);   % 理论：J_0(2πf_d*τ)
    plot(lags, R{1}, 'r--', 'LineWidth', ln_wdt);    % 仿真结果
    xlabel('Normalized Time: $f_d \tau$', 'Interpreter', 'latex', 'FontSize', f_size);
    ylabel('$R_{g_c,g_c}(\tau)$', 'Interpreter', 'latex', 'FontSize', f_size);
    xlim([0, min(15, t_sim*fd)]);
    legend({'Reference (Ideal)', sprintf('Jakes ($M=%d$)', M)}, ...
           'Interpreter', 'latex', 'FontSize', f_size);
    title('同相分量自相关函数', 'FontSize', f_size);
    
    %% Figure 5: 同相与正交分量互相关函数
    %--------------------------------------------------------------------------
    % 验证：I和Q分量是否正交（互相关应该为0）
    % 理论值：R_{g_c,g_s}(τ) = 0（I和Q分量应该正交）
    %--------------------------------------------------------------------------
    figure('Name', 'Jakes - Second Order Statistic', 'NumberTitle', 'on');
    set(gcf, 'Position', [200, 100, 800, 600]);
    grid on; hold on;
    text(0.02, 0.98, 'Figure 5', 'Units', 'normalized', 'FontSize', 18, ...
         'FontWeight', 'bold', 'VerticalAlignment', 'top', ...
         'BackgroundColor', [1 1 1], 'EdgeColor', [0 0 0], 'Margin', 3);
    
    % 绘制理论曲线（黑色实线，应该为0）和仿真曲线（红色虚线）
    plot(lags, R_id{3}, 'k', 'LineWidth', ln_wdt);   % 理论：0（理想正交）
    plot(lags, R{3}, 'r--', 'LineWidth', ln_wdt);    % 仿真结果（应该接近0）
    xlabel('Normalized Time: $f_d \tau$', 'Interpreter', 'latex', 'FontSize', f_size);
    ylabel('$R_{g_c,g_s}(\tau)$', 'Interpreter', 'latex', 'FontSize', f_size);
    ylim([-0.5, 0.5]);
    xlim([0, min(15, t_sim*fd)]);
    legend({'Reference (Ideal)', sprintf('Jakes ($M=%d$)', M)}, ...
           'Interpreter', 'latex', 'FontSize', f_size);
    title('同相与正交分量互相关函数', 'FontSize', f_size);
    
    %% Figure 6: 时域波形
    %--------------------------------------------------------------------------
    % 显示：Jakes信道的时域特性
    % 上子图：包络幅度 |g(t)| 随时间的变化（瑞利衰落）
    % 下子图：相位 θ(t) 随时间的变化（均匀分布）
    %--------------------------------------------------------------------------
    figure('Name', 'Jakes - Time Domain Waveform', 'NumberTitle', 'on');
    set(gcf, 'Position', [200, 100, 1000, 700]);
    
    %--------------------------------------------------------------------------
    % 上子图：幅度时域波形
    %--------------------------------------------------------------------------
    subplot(2, 1, 1);
    hold on;
    text(0.02, 0.98, 'Figure 6', 'Units', 'normalized', 'FontSize', 18, ...
         'FontWeight', 'bold', 'VerticalAlignment', 'top', ...
         'BackgroundColor', [1 1 1], 'EdgeColor', [0 0 0], 'Margin', 3);
    % 只显示前5000个采样点（或全部，取较小值），避免图形过于密集
    plot(time(1:min(5000, length(time))), abs(g(1:min(5000, length(g)))), ...
         'b-', 'LineWidth', 1);
    grid on;
    xlabel('时间 t (s)', 'FontSize', f_size);
    ylabel('幅度 |g(t)|', 'FontSize', f_size);
    title('Jakes信道时域波形 - 幅度', 'FontSize', f_size);
    xlim([0, min(5, t_sim)]);  % 显示前5秒（或全部仿真时长）
    
    %--------------------------------------------------------------------------
    % 下子图：相位时域波形
    %--------------------------------------------------------------------------
    subplot(2, 1, 2);
    hold on;
    text(0.02, 0.98, 'Figure 6-cont', 'Units', 'normalized', 'FontSize', 18, ...
         'FontWeight', 'bold', 'VerticalAlignment', 'top', ...
         'BackgroundColor', [1 1 1], 'EdgeColor', [0 0 0], 'Margin', 3);
    % 相位归一化到[-1, 1]（单位：π）
    plot(time(1:min(5000, length(time))), angle(g(1:min(5000, length(g))))/pi, ...
         'r-', 'LineWidth', 1);
    grid on;
    xlabel('时间 t (s)', 'FontSize', f_size);
    ylabel('相位 \theta(t)/\pi', 'FontSize', f_size);
    title('Jakes信道时域波形 - 相位', 'FontSize', f_size);
    ylim([-1, 1]);              % 相位范围：[-π, π]（归一化后为[-1, 1]）
    xlim([0, min(5, t_sim)]);   % 显示前5秒（或全部仿真时长）
    
    %% Figure 7: 幅度直方图
    %--------------------------------------------------------------------------
    % 验证：通过直方图验证包络幅度是否服从瑞利分布
    % 理论PDF：f_R(r) = (r/σ²)*exp(-r²/(2σ²))，r ≥ 0
    %--------------------------------------------------------------------------
    figure('Name', 'Jakes - Amplitude Histogram', 'NumberTitle', 'on');
    set(gcf, 'Position', [200, 100, 800, 600]);
    grid on; hold on;
    text(0.02, 0.98, 'Figure 7', 'Units', 'normalized', 'FontSize', 18, ...
         'FontWeight', 'bold', 'VerticalAlignment', 'top', ...
         'BackgroundColor', [1 1 1], 'EdgeColor', [0 0 0], 'Margin', 3);
    
    %--------------------------------------------------------------------------
    % 计算直方图
    %--------------------------------------------------------------------------
    abs_g = abs(g);  % 包络幅度
    [counts, centers] = hist(abs_g, 50);  % 50个bin的直方图
    bin_width = centers(2) - centers(1);   % bin宽度
    pdf_hist = counts / (sum(counts) * bin_width);  % 归一化为PDF
    
    %--------------------------------------------------------------------------
    % 计算理论Rayleigh分布
    %--------------------------------------------------------------------------
    r = linspace(0, max(abs_g)*1.2, 500);  % 幅度范围
    sigma_r = sqrt(mean(abs_g.^2)/2);      % 估计Rayleigh分布的参数σ
    % 理论Rayleigh PDF：f_R(r) = (r/σ²)*exp(-r²/(2σ²))
    pdf_rayleigh = (r / sigma_r^2) .* exp(-r.^2 / (2*sigma_r^2));
    pdf_rayleigh(r < 0) = 0;  % 确保r < 0时PDF为0
    
    % 绘制直方图（蓝色柱状图）和理论曲线（红色实线）
    bar(centers, pdf_hist, 'FaceColor', [0.7 0.7 0.9], 'EdgeColor', 'none', 'BarWidth', 1);
    plot(r, pdf_rayleigh, 'r-', 'LineWidth', ln_wdt);
    xlabel('幅度 r', 'FontSize', f_size);
    ylabel('概率密度 f_{|g|}(r)', 'FontSize', f_size);
    legend({'仿真直方图', '理论Rayleigh分布'}, 'FontSize', f_size-2);
    title('幅度分布验证 - 瑞利分布', 'FontSize', f_size);
    xlim([0, max(abs_g)*1.1]);
    
    %% Figure 8: 多普勒频谱（核心验证）
    %--------------------------------------------------------------------------
    % 验证：功率谱是否为U型谱
    % 理论PSD：S(f) = 1/(π*fd*sqrt(1-(f/fd)²))，|f| < fd
    % U型特征：在f=0处功率最大，向±fd处逐渐减小
    %--------------------------------------------------------------------------
    figure('Name', 'Jakes - Doppler Spectrum', 'NumberTitle', 'on');
    set(gcf, 'Position', [200, 100, 800, 600]);
    grid on; hold on;
    text(0.02, 0.98, 'Figure 8', 'Units', 'normalized', 'FontSize', 18, ...
         'FontWeight', 'bold', 'VerticalAlignment', 'top', ...
         'BackgroundColor', [1 1 1], 'EdgeColor', [0 0 0], 'Margin', 3);
    
    %--------------------------------------------------------------------------
    % 计算功率谱密度（PSD）
    %--------------------------------------------------------------------------
    N_fft = 2^nextpow2(length(g));  % FFT点数（2的幂次，提高计算效率）
    G_fft = fft(g, N_fft);          % FFT变换得到频域表示
    PSD = abs(G_fft).^2;            % 功率谱密度 = |G(f)|²
    PSD = fftshift(PSD);             % 零频移到中心（双边谱）
    
    %--------------------------------------------------------------------------
    % 计算频率轴
    %--------------------------------------------------------------------------
    % 频率分辨率：Δf = 1/(N_fft * delta_t)
    f = (-N_fft/2:N_fft/2-1) / (N_fft * delta_t);  % 频率轴（Hz）
    f_norm = f / fd;                                % 归一化频率（f/fd）
    
    %--------------------------------------------------------------------------
    % 计算理论U型频谱
    %--------------------------------------------------------------------------
    % 理论公式：S(f) = 1/(π*fd*sqrt(1-(f/fd)²))，|f| < fd
    f_theory = linspace(-0.99, 0.99, 1000);  % 归一化频率范围（避免奇点）
    PSD_theory = 1 ./ (pi * sqrt(1 - f_theory.^2));  % 理论U型谱
    PSD_theory = PSD_theory / max(PSD_theory);       % 归一化到最大值
    
    %--------------------------------------------------------------------------
    % 归一化仿真频谱
    %--------------------------------------------------------------------------
    PSD_norm = PSD / max(PSD);  % 归一化到最大值，便于对比
    
    %--------------------------------------------------------------------------
    % 绘制频谱（只显示多普勒频率范围内）
    %--------------------------------------------------------------------------
    idx_plot = abs(f_norm) <= 1;  % 只显示|f/fd| ≤ 1的范围
    plot(f_norm(idx_plot), 10*log10(PSD_norm(idx_plot) + eps), 'b-', 'LineWidth', ln_wdt);  % 仿真结果（蓝色）
    plot(f_theory, 10*log10(PSD_theory + eps), 'r--', 'LineWidth', ln_wdt);                % 理论U型谱（红色虚线）
    xlabel('归一化频率 $f/f_d$', 'Interpreter', 'latex', 'FontSize', f_size);
    ylabel('归一化功率谱密度 (dB)', 'FontSize', f_size);
    xlim([-1, 1]);
    ylim([-40, 5]);
    legend({'仿真结果', '理论U型频谱'}, 'FontSize', f_size-2);
    title('多普勒功率谱密度 - U型频谱验证', 'FontSize', f_size);
    
    %% Figure 9: 自相关函数J₀验证（详细视图）
    %--------------------------------------------------------------------------
    % 验证：自相关函数的前几个零点是否与J₀函数一致
    % J₀函数的第一个零点在2.4048，第二个在5.5201，第三个在8.6537...
    % 详细视图便于观察前几个零点的位置
    %--------------------------------------------------------------------------
    figure('Name', 'Jakes - Autocorrelation J0 Verification', 'NumberTitle', 'on');
    set(gcf, 'Position', [200, 100, 800, 600]);
    grid on; hold on;
    text(0.02, 0.98, 'Figure 9', 'Units', 'normalized', 'FontSize', 18, ...
         'FontWeight', 'bold', 'VerticalAlignment', 'top', ...
         'BackgroundColor', [1 1 1], 'EdgeColor', [0 0 0], 'Margin', 3);
    
    idx_detail = lags <= 5;  % 只显示归一化时间延迟≤5的部分
    % 绘制理论曲线（黑色实线）和仿真曲线（红色虚线）
    plot(lags(idx_detail), R_id{5}(idx_detail), 'k', 'LineWidth', ln_wdt);  % 理论：2*J_0
    plot(lags(idx_detail), R{5}(idx_detail), 'r--', 'LineWidth', ln_wdt);    % 仿真结果
    xlabel('归一化时间 $f_d \tau$', 'Interpreter', 'latex', 'FontSize', f_size);
    ylabel('自相关函数 $R_{g,g}(\tau)$', 'Interpreter', 'latex', 'FontSize', f_size);
    xlim([0, 5]);
    legend({'理论 $J_0(2\pi f_d \tau)$', '仿真结果'}, ...
           'Interpreter', 'latex', 'FontSize', f_size-2);
    title('自相关函数J₀验证（详细视图）', 'FontSize', f_size);
    
    %% Figure 10: 实部/虚部分布
    %--------------------------------------------------------------------------
    % 验证：I和Q分量是否分别服从高斯分布
    % 理论：g_c(t)和g_s(t)应该都是零均值高斯过程
    % 理论PDF：f(x) = (1/(σ√(2π)))*exp(-(x-μ)²/(2σ²))
    %--------------------------------------------------------------------------
    figure('Name', 'Jakes - Real/Imaginary Distribution', 'NumberTitle', 'on');
    set(gcf, 'Position', [200, 100, 1200, 600]);
    
    %--------------------------------------------------------------------------
    % 左子图：实部分布（I分量）
    %--------------------------------------------------------------------------
    subplot(1, 2, 1);
    hold on;
    text(0.02, 0.98, 'Figure 10', 'Units', 'normalized', 'FontSize', 18, ...
         'FontWeight', 'bold', 'VerticalAlignment', 'top', ...
         'BackgroundColor', [1 1 1], 'EdgeColor', [0 0 0], 'Margin', 3);
    
    gc = real(g);  % 同相分量（实部）
    [counts_gc, centers_gc] = hist(gc, 50);  % 50个bin的直方图
    bin_width_gc = centers_gc(2) - centers_gc(1);
    pdf_hist_gc = counts_gc / (sum(counts_gc) * bin_width_gc);  % 归一化为PDF
    
    % 计算理论高斯分布
    x_gc = linspace(min(gc), max(gc), 500);
    mu_gc = mean(gc);      % 均值（应该接近0）
    sigma_gc = std(gc);    % 标准差
    % 理论高斯PDF：f(x) = (1/(σ√(2π)))*exp(-(x-μ)²/(2σ²))
    pdf_gauss_gc = (1/(sigma_gc*sqrt(2*pi))) * exp(-(x_gc-mu_gc).^2/(2*sigma_gc^2));
    
    % 绘制直方图（绿色柱状图）和理论曲线（红色实线）
    bar(centers_gc, pdf_hist_gc, 'FaceColor', [0.7 0.9 0.7], 'EdgeColor', 'none', 'BarWidth', 1);
    plot(x_gc, pdf_gauss_gc, 'r-', 'LineWidth', ln_wdt);
    xlabel('实部 $g_c$', 'Interpreter', 'latex', 'FontSize', f_size);
    ylabel('概率密度 $f_{g_c}(x)$', 'Interpreter', 'latex', 'FontSize', f_size);
    legend({'仿真直方图', '理论高斯分布'}, 'FontSize', f_size-2);
    title('实部分布验证', 'FontSize', f_size);
    grid on;
    
    %--------------------------------------------------------------------------
    % 右子图：虚部分布（Q分量）
    %--------------------------------------------------------------------------
    subplot(1, 2, 2);
    hold on;
    text(0.02, 0.98, 'Figure 10-cont', 'Units', 'normalized', 'FontSize', 18, ...
         'FontWeight', 'bold', 'VerticalAlignment', 'top', ...
         'BackgroundColor', [1 1 1], 'EdgeColor', [0 0 0], 'Margin', 3);
    
    gs = imag(g);  % 正交分量（虚部）
    [counts_gs, centers_gs] = hist(gs, 50);  % 50个bin的直方图
    bin_width_gs = centers_gs(2) - centers_gs(1);
    pdf_hist_gs = counts_gs / (sum(counts_gs) * bin_width_gs);  % 归一化为PDF
    
    % 计算理论高斯分布
    x_gs = linspace(min(gs), max(gs), 500);
    mu_gs = mean(gs);      % 均值（应该接近0）
    sigma_gs = std(gs);    % 标准差
    % 理论高斯PDF：f(x) = (1/(σ√(2π)))*exp(-(x-μ)²/(2σ²))
    pdf_gauss_gs = (1/(sigma_gs*sqrt(2*pi))) * exp(-(x_gs-mu_gs).^2/(2*sigma_gs^2));
    
    % 绘制直方图（红色柱状图）和理论曲线（红色实线）
    bar(centers_gs, pdf_hist_gs, 'FaceColor', [0.9 0.7 0.7], 'EdgeColor', 'none', 'BarWidth', 1);
    plot(x_gs, pdf_gauss_gs, 'r-', 'LineWidth', ln_wdt);
    xlabel('虚部 $g_s$', 'Interpreter', 'latex', 'FontSize', f_size);
    ylabel('概率密度 $f_{g_s}(x)$', 'Interpreter', 'latex', 'FontSize', f_size);
    legend({'仿真直方图', '理论高斯分布'}, 'FontSize', f_size-2);
    title('虚部分布验证', 'FontSize', f_size);
    grid on;
    
    fprintf('已生成10个验证图表！\n');
end

