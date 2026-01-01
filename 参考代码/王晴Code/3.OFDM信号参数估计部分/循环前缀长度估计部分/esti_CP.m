%% --------------------------------------------------------------------------------------------- %%
% 函数：估计OFDM信号循环前缀长度
% 输入：ReData：OFDM信号
%       N_Ts：有效符号长度
%       L_win：相关窗长，数值上等于前一次CP估计结果，用于指导下一次估计时的相关计算，若首次估计则赋0
%       N_win：相关窗口个数，仅在精估计时使用
%       mode：函数运行模式，分为run（运行）和debug（调试）两种，debug模式下会输出中间结果图
% 输出：N_CPs：每个OFDM符号中CP的长度

function [N_CPs_ori, N_CPs] = esti_CP(ReData, N_Ts, L_win, N_win, mode)
%     addpath('../utils/');
    if nargin < 5
        mode = 'run';
        if nargin < 4
            N_win = 5;
            if nargin < 3
                L_win = round(N_Ts*0.05);
            end
        end
    end
    if L_win == 0
        L_win = round(N_Ts*0.05); % 如果首次估计则窗长设小一点
    end
    ReData = ReData(:);             % 转为列向量
    Padding = zeros(N_Ts, 1);       % 补零，便于后续寻峰计算
    r = [Padding; ReData; Padding]; % 拼接
    
    %% 直接计算相关性，得到循环前缀长度的粗略估计结果
    res = zeros(1,length(r) - N_Ts - L_win + 1);  % 相关性结果记录
    for i = 1: length(res)
        x = r(i : i + L_win - 1);               % 从r上截取L_win个数据
        y = r(i + N_Ts : i + N_Ts + L_win - 1); % 相隔N_Ts个位置截取L_win个数据
        corr_tmp = x'* y;   % 计算相关性
        if corr_tmp == 0
            res(i) = 0; % 防止由于0/0导致相关系数出现NAN
        else
            res(i) = corr_tmp/vecnorm(x, 2)/vecnorm(y, 2);  % 计算Pearson相关系数
        end
    end
    res = abs(real(res)); % 取实部
    p_locs = AMPD(res, floor(N_Ts*0.9), 0.6*max(res)); % 寻峰
    if strcmp(mode, 'debug')
        figure;plot(res,'-o','MarkerIndices',p_locs,'MarkerFaceColor', 'red', 'MarkerSize',3); % 画图，红点标识峰值位置
        xlabel('\fontname{宋体}位移量', 'Fontsize', 12);
        ylabel('\fontname{宋体}相关系数', 'Fontsize', 12);
        set(gca, 'fontname', 'Times New Roman', 'fontsize', 12, 'looseInset', [0 0 0.03 0.05]);
        axis([0 length(res) 0 1]);
        set(gca,'xtick',0:floor(0.2*length(res)):length(res));
        set(gca,'ytick',0:0.1:1);
    end
    N_CPs_ori = diff(p_locs) - N_Ts; % 做差，减去有效符号长度得到CP估计结果
    %% 多段累加，精细估计循环前缀长度
    N_CP = round(mean(N_CPs_ori));    % 参数设置，令窗长与CP长度相近
    L_win = round(N_CP*0.5);
    wins = zeros(1, L_win*N_win);
    % 各窗间隔等于符号总长度，即有效符号长度加循环前缀长度，确保一个窗对应一个符号
    for i = 1: N_win
        wins((i - 1) * L_win + 1: i * L_win) = (i - 1) * (N_Ts + N_CP) + 1 : (i - 1) * (N_Ts + N_CP) + L_win;
    end
    res = zeros(1,length(r) - (N_win - 1) * (N_Ts + N_CP) - N_Ts - L_win + 1);  % 相关性结果记录
    for i = 1: length(res)
        x = r(i + wins - 1);
        y = r(i + N_Ts + wins - 1);
        corr_tmp = x'* y;
        if corr_tmp == 0
            res(i) = 0; % 防止由于0/0导致相关系数出现NAN
        else
            res(i) = corr_tmp/vecnorm(x, 2)/vecnorm(y, 2); % 计算Pearson相关系数
        end
    end
    res = abs(real(res)); % 取实部
    p_locs = AMPD(res, floor(N_Ts*0.9), 0.6*max(res)); % 寻峰
    if strcmp(mode, 'debug')  % debug模式下画图
        figure;plot(res,'-o','MarkerIndices',p_locs,'MarkerFaceColor', 'red', 'MarkerSize',3); % 画图，红点标识峰值位置
        xlabel('\fontname{宋体}位移量', 'Fontsize', 12);
        ylabel('\fontname{宋体}相关系数', 'Fontsize', 12);
        set(gca, 'fontname', 'Times New Roman', 'fontsize', 12, 'looseInset', [0 0 0.03 0.05]);
        axis([0 length(res) 0 1]);
        set(gca,'xtick',0:floor(0.2*length(res)):length(res));
        set(gca,'ytick',0:0.1:1);
    end
    N_CPs = diff(p_locs) - N_Ts;  % 相邻元素做差并减去有效符号长度
end