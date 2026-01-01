function qpsk= QPSK()
%------------------ 初始化 --------------------
chip=randi([0,1],1,32);
fc = 10e6;          %载波频率
fs = 40e6;          %采样频率
L=200;             % 每个码元的采样点数,即过采样倍数
M=length(chip);    % 码元数
N=M*L;             % 采样点数
Rb=200000;            % 码元速率，即信息速率:2Mbps
Tb=1/Rb;           % 码元周期:0.5us
dt=Tb/L;           % 采样间隔 
T=N*dt;            % 总时间
N_samples=T*fs;    %每个周期的采样点数
dt = 1/fs;
t = 0:dt:T-dt;

for n=1:length(chip)
    if chip(n)==1
        ak(n)=+1;
    elseif chip(n)==0;
        ak(n)=-1;
    end
end
            
% 矩阵变化，使用reshape函数转换
IQ=reshape(ak,2,length(ak)/2);
I=IQ(1,:);   
Q=IQ(2,:);

Isig=[];
Qsig=[];
for ii=1:length(IQ)
    if I(ii)==1;
        if Q(ii)==1;
           Isig(ii)=1;
           Qsig(ii)=-1;
        else Isig(ii)=-1;
             Qsig(ii)=-1;
        end
    else
        if Q(ii)==1;
           Isig(ii)=-1;
           Qsig(ii)=1;
        else Isig(ii)=1;
             Qsig(ii)=1;
        end
    end
end


I=Isig;
Q=Qsig;          
I=[I;I];
I=reshape(I,1,2*length(I));
% I=[1,I];   %延时一个时间单位，图显示的是对Q路延时，对应I的第一位和Q的第二位相加，在本代码中反映的就是这样

Q=[Q;Q];
Q=reshape(Q,1,2*length(Q));

%将码元转换为码元时间序列，使用p;
It=[];Qt=[];                 %Q,I两路信号的时域序列
%码元重复一次，直接加到后面，使得I,Q两路信号码元多一倍
m=min(length(I),length(Q));
for i=1:L
    It=[It;I];
    Qt=[Qt;Q];
end

It1=It(:);
It2=It1';
It=It2(1:m*L);
Qt1=Qt(:);
Qt2=Qt1';
Qt=Qt2(1:m*L);
c1=1;
s1=1;
c2=cos(2*pi*fc*t);
s2=sin(2*pi*fc*t);
%对两路信号进行加权处理
Q_out1=Qt.*s1;      %加权sin(pi*t/2Ts)的载波信号
I_out1=It.*c1;      %加权cos(pi*t/2Ts)的载波信号
Q_out=Q_out1.*s2;%调制I路的载波信号  
I_out=I_out1.*c2;%调制Q路的载波信号
qpsk=I_out-1i*(Q_out);
