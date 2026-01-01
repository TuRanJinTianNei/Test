function y = QAM16()
N=64; % 符号数
df=50; % 频率分辨率
M=16; % 16QAM调制
x=randi([0,M-1],1,N);% 随机生成数据 (信号源)
%============= Use 16-QAM modulation=====================% QAM 调制
% fc=10e6; % 载波频率
% fs=40e6; % 采样频率
% ts=1/fs; % 采样间隔

% t=[0:ts:N*df*ts]; % 时间序列

y = qammod(x,M);  % 使用新的 qammod 函数替代已删除的 modem.qammod
y = repmat(y,df,1);
y = reshape(y,1,df*N);
P=std(y);           %计算y的标准差，使用函数std来计算，mean和abs函数计算的是不同的标准差
y=y/P;              %P1和P虽然相同，但一个是除以标准差，一个是除以模型的均值,一个是除以平方和的开方

% c=exp(j*2*pi*fc.*t); % 载波信号
% y=y.*c(1:length(y));
