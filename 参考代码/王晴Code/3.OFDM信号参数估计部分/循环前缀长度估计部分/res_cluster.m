function N_CP = res_cluster(N_CPs)
    centers = cluster_test(N_CPs, 1, 0.5, 16); % 聚类
    [classes, counts] = identify_sig(N_CPs, centers); % 归类
    N_CP = zeros(1, length(centers));
    for i = 1: length(N_CPs)
        N_CP(classes(i)) = N_CP(classes(i)) + N_CPs(i);
    end
    N_CP = round(N_CP ./ counts);  % 每类计算均值
    CPs = zeros(1, length(N_CPs));
    for i = 1: length(CPs)
        CPs(i) = N_CP(classes(i)); % 每个符号CP长度
    end
end