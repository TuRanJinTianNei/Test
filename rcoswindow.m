function[rcosw]=rcoswindow(beta,Ts)

%定义升余弦窗，其中beta为滚降系数，Ts为包含循环前缀的OFDM符号的长度
t=0:(1+beta)*Ts;
rcosw=zeros(1,(1+beta)*Ts);
for i=1:beta*Ts
    rcosw(i)=0.5+0.5*cos(pi+t(i)*pi/(beta*Ts));
end

for i=beta*Ts+1:Ts
    rcosw(i)=1;
end

for i=(Ts+1):(1+beta)*Ts+1
    rcosw(i-1)=0.5+0.5*cos((t(i)-Ts)*pi/(beta*Ts));
end
% 变换为列向量
rcosw=rcosw';


 


