function [B_welch, B_ar] = method(sig1_chnl, fs, snr)
    %**************************************************************************
    % AR和Welch方法估计OFDM信号带宽
    % 
    % 输入参数:
    %   sig1_chnl  - 输入的OFDM信号（可以是复数或实数信号）
    %   fs         - 采样频率 (Hz)
    %   snr        - 信噪比 (dB)，仅用于AR模型法带宽提取时的阈值选择
    %
    % 输出参数:
    %   B_welch    - Welch算法估计的带宽 (Hz)
    %   B_ar       - AR模型法估计的带宽 (Hz)
    %
    % 功能说明:
    %   1. 将输入信号转换为实信号（如果是复数，取实部）
    %   2. 使用AR模型法（Burg算法）和Welch算法进行功率谱估计
    %   3. 从功率谱中提取带宽（根据SNR选择不同的阈值点）
    %
    % 创建日期: 2025.12.10
    % 提取自: test1.m 和 PSD_OFDM_rayleigh.m
    % 修改日期: 2025.12.10 - 集成所有依赖函数，使其成为独立函数
    %           2025.12.10 - 删除上变频、信道仿真和添加噪声部分
    %**************************************************************************
    
    % 将输入信号转换为实信号（如果是复数，取实部）
    if ~isreal(sig1_chnl)
        sig3_chnl = real(sig1_chnl);
    else
        sig3_chnl = sig1_chnl;
    end
    
    % 步骤1: 功率谱估计
    % 1.1 使用AR模型法（Burg算法）进行功率谱估计
    [Pxx1, f, ~] = Burg_local(sig3_chnl, fs, 'AIC');
    
    % 1.2 使用Welch算法进行功率谱估计
    [Pxx2, f1] = pwelch(sig3_chnl, hanning(100), 55, 4096*2, fs);
    
    % 步骤2: 功率谱归一化和单位转换（Welch算法）
    Pxx22 = Pxx2;
    Pxx22 = Pxx22 / min(Pxx22);        % 归一化处理
    Pxx22 = 10*log10(Pxx22);            % 转换为dB单位
    Pxx22 = Pxx22 - max(Pxx22);         % 归一化到最大值
    
    % 步骤3: 从功率谱中提取带宽
    % 6.1 Welch算法带宽提取（使用-3dB点）
    L1 = ceil(length(Pxx22)/2);
    P1 = Pxx22(1:L1, 1);
    P2 = Pxx22(L1:end, 1);
    [as1, ~] = Proximate_local(-3, P1);    % 取最接近-3dB处信号的f值
    band1 = f1(as1);
    [as2, ~] = Proximate_local(-3, P2);
    band2 = f1(as2 + L1 - 1);
    B_welch = abs(band1 - band2);
    
    % 6.2 AR模型法带宽提取（根据SNR选择不同的阈值点）
    L2 = ceil(length(Pxx1)/2);
    P3 = Pxx1(1:L2, 1);
    P4 = Pxx1(L2:end, 1);
    
    if snr > 4
        % SNR > 4dB: 使用-6dB点
        [as3, ~] = Proximate_local(-6, P3);
        band3 = f(as3);
        [as4, ~] = Proximate_local(-6, P4);
        band4 = f(as4 + L2 - 1);
        B_ar = abs(band4 - band3);
    elseif snr > 0 && snr <= 4
        % 0dB < SNR <= 4dB: 使用-5dB点
        [as3, ~] = Proximate_local(-5, P3);
        band3 = f(as3);
        [as4, ~] = Proximate_local(-5, P4);
        band4 = f(as4 + L2 - 1);
        B_ar = abs(band4 - band3);
    else
        % SNR <= 0dB: 使用-3dB点
        [as3, ~] = Proximate_local(-3, P3);
        band3 = f(as3);
        [as4, ~] = Proximate_local(-3, P4);
        band4 = f(as4 + L2 - 1);
        B_ar = abs(band4 - band3);
    end
    
    end
    
    %===========================================================================
    % 以下为集成的本地函数（Local Functions）
    %===========================================================================
    
    function [psdviaBurg, f, p] = Burg_local(x, Fs, varargin)
    % 使用Burg算法实现的AR模型功率谱估计
    % 输入参数:
    %   x        - 输入信号
    %   Fs       - 采样频率
    %   varargin - 可以为数值型（AR模型阶数）或字符串（准则，如'AIC'）
    % 输出参数:
    %   psdviaBurg - 使用Burg算法计算的功率谱值
    %   f          - 频率采样点
    %   p          - 模型阶数
    
    % 根据输入参数类型判断
    if isnumeric(varargin{1}) && isscalar(varargin{1})
        p = varargin{1};
    elseif ischar(varargin{1}) || isstring(varargin{1})
        criterion = varargin{1};
    else
        error('第2个参数必须为数值型或字符串');
    end
    
    x = x(:);
    N = length(x);
    
    % 模型参数估计
    if exist('p', 'var') % p变量是否存在，如果存在则不需要估计，直接使用p值
        [a, E] = computeARpara_local(x, p);
    else % p不存在，需要估计，根据准则criterion
        p = ceil(N/3); % 阶数一般不超出信号长度的1/3
        
        % 计算1到p阶的误差
        [a, E] = computeARpara_local(x, p);
        
        % 计算目标函数的最小值
        kc = 1:p + 1;
        switch criterion
            case 'FPE'
                goalF = E.*(N + (kc + 1))./(N - (kc + 1));
            case 'AIC'
                goalF = N.*log(E) + 2.*kc;
        end
        [~, p] = min(goalF); % p是目标函数最小值的位值，也即准则确定的阶数
        
        % 使用p值重新计算AR模型参数
        [a, E] = computeARpara_local(x, p);
    end
    
    [h, f] = freqz(1, a, 20e5, Fs);
    psdviaBurg = E(end)*abs(h).^2./Fs;
    psdviaBurg = psdviaBurg/abs(max(psdviaBurg));
    psdviaBurg = (10*log10(abs(psdviaBurg)));
    end
    
    function [a, E] = computeARpara_local(x, p)
    % 根据输入信号x和阶数p计算AR模型参数估计（Burg算法）
    % 输入参数:
    %   x - 输入信号
    %   p - AR模型阶数
    % 输出参数:
    %   a - AR模型系数
    %   E - 误差功率序列
    
    N = length(x);
    % 初始值
    ef = x; % 前向预测误差
    eb = x; % 后向预测误差
    a  = 1; % 初始模型参数
    E  = x'*x/N; % 初始误差
    k  = zeros(1, p); % 为反射系数预分配空间，加快循环速度
    E  = [E k]; % 为误差预分配空间，加快速度
    
    for m = 1:p
        % 按照Burg算法步骤，首先计算m阶的反射系数
        efm = ef(2:end); % 前一阶次的前向预测误差
        ebm = eb(1:end - 1); % 前一阶次的后向预测误差
        num = -2.*ebm'*efm;  % 反射系数的分子项
        den = efm'*efm + ebm'*ebm; % 反射系数的分母项
        k(m) = num./den; % 当前阶次的反射系数
        
        % 更新前后向预测误差
        ef = efm + k(m)*ebm;
        eb = ebm + conj(k(m))*efm;
        
        % 更新模型系数a
        a = [a; 0] + k(m)*[0; conj(flipud(a))];
        
        % 当前阶次的误差功率
        E(m + 1) = (1 - conj(k(m))*k(m))*E(m);
    end
    end
    
    function [as1, as11] = Proximate_local(b, aa)
    % 查找数组中与目标值最接近的元素索引
    % 输入参数:
    %   b  - 要找的参数值
    %   aa - 要搜索的数组
    % 输出参数:
    %   as1  - 最接近b的第一个元素的位置索引
    %   as11 - 所有最接近b的元素位置索引
    
    % 判断aa是几维数组，统一转换为一维
    a = aa(:);                    % 将数组转换为一维数组
    ab = (a(:) - b)';            % 计算数组a和b的差值
    abc = abs(ab);               
    abc = sort(abc);             % 绝对值取最小值，排序
    
    % 找到与b值最接近的第一个数组元素的位置（最接近b的）
    [as1, as11] = find(abs((a(:) - b)) == abc(1, 1));
    end
    