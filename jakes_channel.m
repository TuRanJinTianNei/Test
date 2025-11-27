function h_jakes = jakes_channel(fd, N0, channel_length, fs)
%JAKES_CHANNEL 生成Jakes模型Rayleigh衰落信道
%   从jakes_simulator_standalone.m提取的Jakes信道生成函数
%
%   输入参数：
%       fd           - 最大多普勒频移（Hz）
%       N0           - 散射体数量（Jakes模型中的基础振荡器数量）
%       channel_length - 信道样本长度（样本数）
%       fs           - 采样频率（Hz）
%
%   输出参数：
%       h_jakes      - 复信道冲激响应（时域），行向量
%
%   实现方法：
%       自定义函数实现，基于Jakes模型，使用多个正弦振荡器叠加模拟瑞利衰落
%       核心算法（gen_jakes_channel_core函数）：
%       1. 初始化特殊振荡器（n=0）：对应最高多普勒频移（f = ±fd），计算I/Q分量的初始值
%       2. 叠加M个振荡器（n=1,2,...,M）：每个振荡器对应不同的多普勒频移
%          - 角度：beta_n = n*π/M（等间隔分布）
%          - 频率：w_n = w_0*cos(2π*n/N)，其中N = 4*M + 2
%          - 叠加到I和Q分量
%       3. 归一化：归一化系数2/sqrt(N)，生成复衰落过程g(t) = g_c(t) + j*g_s(t)
%
%   特点：
%       - 返回复信道冲激响应向量（时域）
%       - 通过多个振荡器叠加模拟多径传播
%       - 最终归一化确保平均功率为1
%       - 与jakesChannel.m（基于comm.RayleighChannel）不同，这是自定义实现
%
%   说明：
%       本函数基于jakes_simulator_standalone.m中的gen_jakes_channel函数
%       进行了适配，使其能够直接用于test5.m等OFDM仿真程序
%
%   使用示例：
%       fd = 100;
%       N0 = 34;
%       channel_length = 1000;
%       fs = 15.36e6;
%       h = jakes_channel(fd, N0, channel_length, fs);

    % 将N0映射到M（基础振荡器数量参数）
    % 在jakes_simulator_standalone中，总路径数 N = 4*M + 2
    % 为了保持兼容性，我们使用 M = N0（因为实际实现中通常只用M个路径）
    M = N0;
    
    % 生成时间向量
    dt = 1/fs;  % 采样间隔
    time = (0:channel_length-1)' * dt;  % 时间向量（列向量，单位：秒）
    
    % 调用核心Jakes信道生成函数
    g = gen_jakes_channel_core(M, fd, time);
    
    % 转换为行向量（与test5.m中的使用方式一致）
    h_jakes = g(:)';
    
    % 归一化信道功率（确保平均功率为1）
    h_jakes = h_jakes / sqrt(mean(abs(h_jakes).^2));
end

function g = gen_jakes_channel_core(M, fd, time)
%GEN_JAKES_CHANNEL_CORE Jakes模型核心生成函数
%   从jakes_simulator_standalone.m提取的核心实现
%
%   输入参数：
%       M    - 基础振荡器数量参数（总路径数 N = 4*M + 2）
%       fd   - 最大多普勒频移（Hz）
%       time - 时间向量（列向量，单位：秒）
%
%   输出参数：
%       g    - 复衰落过程 g(t) = g_c(t) + j*g_s(t)

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
end

