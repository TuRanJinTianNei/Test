function [center] = cluster_test(y, ra, rb, ratio)
    center = zeros(1,1);%聚类中心
    seq = zeros(1,1);
    len = length(y);
    D = zeros(len, 1);%为每一点量化一个值，当该点附近存在点越多时，该值越大
    y = abs(y);
    My = max(y);
    y = y./My;
    for i = 1:len
        for j = 1:len
            if i ~= j
%                 D(i) = D(i) + (ra^4/((y(i) - y(j))^4 + 1e-6));%定义函数
                D(i) = D(i) + exp(-(y(i) - y(j))^2/ra^2);
            end
        end
    end
    if ratio == 0
        th = 0;
    else
        th = max(D)/ratio;
    end
    index = 0;
    while index < 2
        [tmp_max, tmp] = max(D);
        if tmp_max < th
            break;
        end
        index = index + 1;
        seq(index) = tmp(1);
        center(index) = y(seq(index));
        for i = 1:len
            D(i) = D(i) - tmp_max * exp(-(y(i) - center(index))^2/rb^2);
%             D(i) = D(i) - tmp_max*(rb^4/(((y(i) - center(index))^2)^2 + 1e-6));
        end
    end
    center = sort(center.*My, "ascend");
end