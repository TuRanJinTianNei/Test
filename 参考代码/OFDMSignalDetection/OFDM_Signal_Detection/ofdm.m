function y2=ofdm(N,para,ratio)
%**************************************************************************
%功能：ofdm信号生成函数
%N:符号个数
%para:子载波数目
%**************************************************************************

M = 16;   %16QAM调制
Signal = randi([0,M-1],1,para*N);
QAM_out = qammod(Signal,M);  % 使用新的 qammod 函数替代已删除的 modem.qammod
x = reshape(QAM_out,para,N);       %矩阵转换
y = ifft(x);                        %对调制后的数据进行ifft变换
y1 = [y(end-round(length(y(:,1))*ratio)+1:end,:);y];%添加循环前缀，比例为1/4
y2 = reshape(y1,1,(length(y1(:,1)))*N);   %矩阵转换
P = std(y2);
y2 = y2/P;


