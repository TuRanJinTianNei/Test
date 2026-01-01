%random_binary.m
%生成随机二进制信源数据
function [info]=random_binary(N)
 if nargin == 0,     %如果没有输入参数，则默认指定信息序列为10000个码元
  N=10000;
end;
for i=1:N,
  temp=rand;             
  if (temp<0.5),
    info(i)=0;         % 1/2的概率取值为0
  else
    info(i)=1;         % 1/2的概率取值为1
  end
end;