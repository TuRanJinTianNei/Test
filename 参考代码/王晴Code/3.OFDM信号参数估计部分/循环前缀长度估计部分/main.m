fprintf('----------------------\n');
%% ----------------------------------------参数设置----------------------------------------------------- %%
subCarry = 256;             % 有效子载波数%%四倍过采样
fftLen = 1024;              % FFT长度为1024

symbN = 100;                 % 一帧中OFDM符号个数
cpLen = 256;                % 循环前缀长度
csLen = 20;                 % 循环后缀长度

% 生成不同循环前缀长度的OFDM符号
% cpLen1 = 72;
% cpLen2 = 88;
% csLen1 = 20;
% csLen2 = 20;
% symbN1 = 13;
% symbN2 = 1;
SNR = -3;                   % 信噪比取值，单位dB
k = 2;                      % 调制阶数
%% 生成信号
ReData = generate_OFDM(subCarry, fftLen, cpLen, csLen, symbN, SNR, k);
% ReData = generate_OFDM2(subCarry, fftLen, cpLen1, cpLen2, csLen1, csLen2, symbN1, symbN2, symbN, SNR, k);
% h = channel_test(cpLen);
% ReData = conv(ReData, h);
% fd = exp(1j*2*pi*20*(1:length(ReData))/8192); 
% Redata = ReData .* fd;
%% 参数估计（有效符号长度、循环前缀长度、FFT点数）
L = 10000;     % 参数设置，保证有效符号长度小于此值
N_Ts = esti_Ts(ReData, L, 'debug');
fprintf('# 有效符号长度：%g\n', N_Ts);

L_win = 0;  % 参数设置
N_win = 10;  % 多段累加的段数
[N_CP_ori, N_CP] = esti_CP(ReData, N_Ts, L_win, N_win, 'debug');
% 取均值
CPs_mean = mean(N_CP);
CPs_ori_mean = mean(N_CP_ori);
% 聚类
% N_CPs_ori = res_cluster(N_CP_ori);
% N_CPs = res_cluster(N_CP);
fprintf('# 循环前缀长度：%g\n', round(CPs_ori_mean));
fprintf('# 循环前缀长度：%g（多段累加）\n', round(CPs_mean));

fprintf('----------------------\n');
