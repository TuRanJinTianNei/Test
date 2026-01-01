%% 寻峰函数
% 输入参数： % res：待寻峰信号
            % max_index：寻找极值点的窗口范围
            % flag：给定的阈值条件，所找极值需要超过的门限值
% 输出：合并的突发信号
function [P_locs] = AMPD(res, max_index, flag)
    if nargin < 3
        flag = 0;
        if nargin < 2
            max_index = 0;
        end
    end
    p_data = zeros(length(res), 1);                  % 峰值个数记录向量
    if max_index == 0                                % 如果没给出对应参数，则自适应计算该参数
        arr_sum = zeros(ceil(length(res)/2) - 1, 1); % 极值个数向量
        for k = 1: ceil(length(res)/2) - 1           % 遍历窗长k，寻找极值
            for i = k + 1 : length(p_data) - k
                if abs(res(i)) > abs(res(i - k)) && abs(res(i)) > abs(res(i + k)) %朴实的极值条件
                    arr_sum(k) = arr_sum(k) + 1;
                end
            end
        end
        [~, max_index] = max(arr_sum); %找到极值点个数最多对应的窗长
    end
    if flag == 0
        for k = 1:max_index
            for i = max_index + 1 : length(p_data) - max_index
                if abs(res(i)) > abs(res(i - k)) && abs(res(i)) > abs(res(i + k))
                    p_data(i) = p_data(i) + 1; %找到均满足极值条件的点，即峰值点
                end
            end
        end
    else
        for k = 1:max_index
            for i = max_index + 1 : length(p_data) - max_index
                if abs(res(i)) - abs(res(i - k)) >= 0 && abs(res(i)) - abs(res(i + k)) > 0 && abs(res(i)) > flag
                    p_data(i) = p_data(i) + 1; %找到均满足极值条件的点，即峰值点
                end
            end
        end
    end
    % 找到所有满足条件的点：
    % 对于p_data(i), 若在[i-max_index, i+max_index]内p_data(i)最大，则是峰值。
    P_locs = find(p_data == max_index);
end