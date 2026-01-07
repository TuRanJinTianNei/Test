function [ff] = cyclic_spectrum(x, N, fs, M,K)
%*******************************************************************
% 功能：平滑周期图法计算循环谱
% x 输入信号
% N 是 循环谱估计的数据长度，必须小于等于输入信号的长度
% fs 是 采样频率, 频率范围为-fs/2到fs/2
% M 平滑长度, 时间分辨率*频率分辨率=M
%*******************************************************************
win = 'hamming'; % 平滑窗函数

d_alpha = fs/(N); % 1/时间分辨率=循环频率分辨率
alpha = -fs:d_alpha:fs; % 循环频率, 分辨率=1/时间分辨率
a_len = length(alpha); % 循环频率取值个数

f_len = floor(N/M-1)+1; % 平滑周期图的个数, 频率分辨率不变
f = -(fs/2-d_alpha*floor(M/2)) + d_alpha*M*(0:f_len-1); % 频率采样点位置

S = zeros(a_len, f_len); % 初始化循环谱矩阵
i = 1; 

%%% 信号fft变换 %%%
X = fftshift(fft(x(1:N))); 
X = X';

%%% 遍历循环频率取值 %%%
for alfa = alpha

    interval_f_N = round(abs(alfa)/d_alpha); % 循环频率对应频率的间隔
    f_N = floor((N-interval_f_N-M)/M)+1; % 平滑周期的个数
    
    %%% 生成平滑窗函数 %%%
    %g = feval(win, M); 
    %window_M = g(:, ones(f_N,1));
    
    %%% 频率域平滑模型 %%%
    t = 1:M*f_N;
    t = reshape(t, M, f_N);
    
    %%% 计算X1,X2 %%%
    %X1 = X(t).*window_M;
    %X2 = X(t+interval_f_N).*window_M; 
    X1 = X(t);
    X2 = X(t+interval_f_N); 
    %%% 计算循环相关 %%%
    St = conj(X1).*X2;
    St = mean(St, 1); % 平滑平均
    S(i, floor((f_len-f_N)/2)+(1:f_N)) = St/N; % 计算平滑周期图并存储以便绘图
    i = i+1;
    
end
%%% 遍历循环频率取值结束 %%%
 
%%% ѭ������ͼ %%%
if K==1
    figure;
    mesh(f, alpha, abs(S)); 
    axis tight;
    title('QAM-OFDM 循环谱');
    xlabel('频率 f (Hz)'); 
    ylabel('循环频率 α');
end

normal_data=S(floor(length(S(:,1))/2)+1,:);
if K==1
    figure;
    plot(f,abs(normal_data))
    xlabel('频率 [Hz]');
    ylabel('S0(f)');
    title('OFDM循环自相关函数α=0的二维切片');
    grid;
end
normal_data =interp(normal_data,100);
ff = interp(f,100);
%***********************************************
%载频估计
%***********************************************
% 针对宽带信号（如OFDM），PSD是平坦的，直接找最大值会导致误差。
% 使用阈值法找到频带范围，取中心作为载频估计。
abs_data = abs(normal_data);
len_half = floor(length(abs_data)/2);

% 负频率部分
neg_part = abs_data(1:len_half);
neg_freq = ff(1:len_half);
max_neg = max(neg_part);
% 使用0.8作为阈值（-2dB左右），取能量集中的区域中心
indices_neg = find(neg_part > max_neg * 0.8);
if ~isempty(indices_neg)
    f1 = abs(mean(neg_freq(indices_neg))); 
else
    [~, idx] = max(neg_part);
    f1 = abs(neg_freq(idx));
end

% 正频率部分
pos_part = abs_data(len_half+1:end);
pos_freq = ff(len_half+1:end);
max_pos = max(pos_part);
indices_pos = find(pos_part > max_pos * 0.8);
if ~isempty(indices_pos)
    f2 = abs(mean(pos_freq(indices_pos)));
else
    [~, idx] = max(pos_part);
    f2 = abs(pos_freq(idx));
end

ff=(f1+f2)/2;

[as2, as22]=Proximate(10e6,f);
data3=S(:,as2);
normal_data3=data3/max(data3);
 if K==1
    figure;
    plot(alpha,abs(normal_data3))
    xlabel('循环频率 [Hz]');
    ylabel('S0(f)');
    title('OFDM循环自相关函数f=fc的二维切片');
    grid;
 end