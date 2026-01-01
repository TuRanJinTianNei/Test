%% -----------------------------------------------------
% 函数：基于循环前缀相关性的有效符号长度估计
% 输入： % r是接收信号，应注意截取完整的OFDM信号，且减少纯噪声段
        % L是信号可能的有效符号长度的最大值
        % fs是信号的采样率
% 输出： % N_Ts是有效符号长度估计值对应的采样点个数
%% -----------------------------------------------------
function [N_Ts] = esti_Ts(r, L, mode)
    
    % 参数默认设置
    if nargin < 3
        mode = 'run';
        if nargin < 2
            L = floor(length(r)/6);
        end
    end
    res = zeros(1, L);  % 相关性结果记录
    L_win = 5*L;        % 参数设置
    r = r(:);           % 转为列向量

    % 循环移位计算相关性
    for i = round(L/50) + 1: length(res) 
        x = r(1 : L_win);
        y = r(1 + i : L_win + i);
        corr_tmp = x'* y;
        if corr_tmp == 0
            res(i) = 0; % 防止由于0/0导致相关系数出现NAN
        else
            res(i) = corr_tmp/vecnorm(x, 2)/vecnorm(y, 2);    % 计算Pearson相关系数
        end
    end
    res = abs(real(res));
    if strcmp(mode, 'debug')  % debug模式下画图
        figure;plot(101:length(res), res(101:end));
        xlabel('\fontname{宋体}位移量', 'Fontsize', 12);
        ylabel('\fontname{宋体}相关函数值', 'Fontsize', 12);
        set(gca, 'fontname', 'Times New Roman', 'fontsize', 12);
    end
    
    [value, index] = max(res);  % 找到相关性最大的值和位置
    p_locs = AMPD(res, round(index/4), value/2);  % 以最大相关性及其所处位置为参数搜索峰值
    N_Ts = p_locs(end);   % 选择对应值最大（数轴上最靠右）的峰值作为有效符号长度估计结果（对于OFDM1信号，有效符号长度为8192，但在4096处同样存在峰值，且峰的高度超过8192处）
end