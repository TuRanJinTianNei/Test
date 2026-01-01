function [ sig_chnl ] = PSD_generate( N)
%生成用于测试的ofdm信号
%N一个帧结构中OFDM信号的个数，每次发送的ofdm信号的个数

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
yy1=[zeros(round(length(x1)/3),N);yy1; zeros(round(2*length(x1)/3),N)];  %前后补零
sig_chnl=reshape(yy1,1,(length(yy1))*N);
P0=std(sig_chnl);
sig_chnl=sig_chnl/P0;
