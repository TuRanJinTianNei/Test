function [ sig_chnl ] = PSD_generate( N)
%生成用于测试的ofdm信号
%N一个帧结构中OFDM信号的个数，每次发送的ofdm信号的个数
%使用窗函数实现符号间平滑过渡，避免频谱泄露

para=128;%设置带内数据子载波数量,有效数据长度
M=16;
Signal=randi([0,M-1],1,para*N);
QAM_out=qammod(Signal,M);  % 使用新的 qammod 函数替代已删除的 modem.qammod
x=reshape(QAM_out,para,N);      %矩阵转换
L=length(x(:,1));
for i=1:N
      x1(:,i)  = [x( 1 : end / 2,i) ;zeros(4*L,1) ;x( end / 2 + 1 : end,i)];  %4倍过采样
end

yy=ifft(x1);                         %对调制后的数据进行ifft变换
yy1=[yy(end-round(length(yy)/4)+1:end,:);yy];%添加循环前缀，比例为1/4

% 移除补零，使用窗函数实现平滑过渡
% 设置过渡窗长度（符号长度的5%）
symbol_len = length(yy1(:,1));
transition_len = max(round(symbol_len * 0.05), 20);  % 过渡区域长度，最小20个点
if transition_len > symbol_len / 4
    transition_len = round(symbol_len / 4);  % 不超过符号长度的1/4
end

% 初始化输出信号
sig_chnl = [];

for i = 1:N
    symbol = yy1(:, i);  % 当前符号
    
    % 创建窗函数实现平滑过渡
    window = ones(length(symbol), 1);
    
    if i == 1
        % 第一个符号：只在末尾加窗（平滑过渡到下一个符号）
        fade_out = hanning(2*transition_len);
        fade_out = fade_out(transition_len+1:end);  % 取后半部分（从0.5到1）
        window(end-transition_len+1:end) = fade_out;
    elseif i == N
        % 最后一个符号：只在开头加窗（从上一个符号平滑过渡）
        fade_in = hanning(2*transition_len);
        fade_in = fade_in(1:transition_len);  % 取前半部分（从0到0.5）
        window(1:transition_len) = fade_in;
    else
        % 中间符号：两端都加窗
        % 开头：从0平滑上升到1（与上一个符号的末尾重叠相加）
        fade_in = hanning(2*transition_len);
        fade_in = fade_in(1:transition_len);
        window(1:transition_len) = fade_in;
        % 末尾：从1平滑下降到0.5（为下一个符号的上升做准备）
        fade_out = hanning(2*transition_len);
        fade_out = fade_out(transition_len+1:end);
        window(end-transition_len+1:end) = fade_out;
    end
    
    % 应用窗函数
    symbol_windowed = symbol .* window;
    
    % 连接符号（平滑过渡，无补零）
    if i == 1
        sig_chnl = symbol_windowed;
    else
        % 重叠相加：将过渡区域重叠相加，实现平滑连接
        overlap_start = length(sig_chnl) - transition_len + 1;
        if overlap_start > 0 && overlap_start <= length(sig_chnl)
            % 重叠区域：上一个符号的末尾 + 当前符号的开头
            sig_chnl(overlap_start:end) = sig_chnl(overlap_start:end) + symbol_windowed(1:transition_len);
            % 添加当前符号的剩余部分
            sig_chnl = [sig_chnl; symbol_windowed(transition_len+1:end)];
        else
            % 如果没有重叠，直接连接
            sig_chnl = [sig_chnl; symbol_windowed];
        end
    end
end

% 转换为行向量
sig_chnl = sig_chnl(:)';
P0=std(sig_chnl);
sig_chnl=sig_chnl/P0;
