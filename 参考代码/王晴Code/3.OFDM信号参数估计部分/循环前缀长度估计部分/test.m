subCarry = 256;             % 有效子载波数%%四倍过采样
fftLen = 1024;              % FFT长度为1024
symbN = 30;                 % 一帧中OFDM符号个数
cpLen = 256;                % 循环前缀长度
csLen = 20;                 % 循环后缀长度
k = 2;                      % 调制阶数
SNR_start = -3;
SNR_end = 5;
repeat = 100;
res_ori = zeros(SNR_end - SNR_start + 1, 1);
res_add = zeros(SNR_end - SNR_start + 1, 1);
for SNR = SNR_start: SNR_end
    tmp1 = zeros(1, repeat);
    tmp2 = zeros(1, repeat);
    for i = 1: repeat 
        ReData = generate_OFDM(subCarry, fftLen, cpLen, csLen, symbN, SNR, k);
        L = 5000;     % 参数设置，保证有效符号长度小于此值
        N_Ts = esti_Ts(ReData, L, 'run');
        fprintf('----------------------\n');
        fprintf('# 有效符号长度：%g\n', N_Ts);
        L_win = 0;  % 参数设置
        N_win = 10;  % 多段累加的段数
        [N_CP_ori, N_CP] = esti_CP(ReData, N_Ts, L_win, N_win, 'run');
        CPs_mean = mean(N_CP);
        CPs_ori_mean = mean(N_CP_ori);
        fprintf('# 循环前缀长度：%g\n', round(CPs_ori_mean));
        fprintf('# 循环前缀长度：%g（多段累加）\n', round(CPs_mean));
        fprintf('----------------------\n');
        tmp1(i) = round(CPs_ori_mean);
        tmp2(i) = round(CPs_mean);
    end
    res_ori(SNR - SNR_start + 1) = mean((tmp1(~isnan(tmp1))-cpLen).^2);
    res_add(SNR - SNR_start + 1) = mean((tmp2(~isnan(tmp2))-cpLen).^2);
end

figure;
f1 = plot(SNR_start: SNR_end, log10(res_ori), 'ro-');
set(f1,'Color','#da5f4f', 'LineWidth', 2.4, 'MarkerSize', 4, 'MarkerFaceColor','#da5f4f');
hold on;
f2 = plot(SNR_start: SNR_end, log10(res_add), 'rp-');
set(f2,'Color','#fa7f6f', 'LineWidth', 2.4, 'MarkerSize', 4, 'MarkerFaceColor','#fa7f6f');
grid on;
axis([SNR_start SNR_end 0 5]);
xlabel('\fontname{Times New Roman}SNR (dB)','fontsize',12);
ylabel('\fontname{宋体}方差','fontsize',12);
legend("ori", "add", "location", "southeast");
set(gca,'fontname','Times New Roman','fontsize',12,'looseInset',[0 0 0.03 0.05]);
set(gca,'xtick',SNR_start:2:SNR_end);
set(gca,'ytick',0:0.5:5);