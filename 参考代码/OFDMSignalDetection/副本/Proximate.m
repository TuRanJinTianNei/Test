function [as1, as11]=Proximate(b,aa)
%%% b是要找的参数值
%%% aa是要找的数组

%%判断aa是几维数组，统一转换为一维
a=aa(:);                                            %%将数组转换为化为一维数组
ab=(a(:)-b)';                                     %%计算数组a和b的差值
abc=abs(ab);                  
abc=sort(abc);                                  %%绝对值取最小值，排序

%%%  二维数组as存储b值最接近的值，在abc中的位置 
%    找到与b值最接近的第一个数组元素的位置（最接近b的）
[as1 ,as11] =find(abs((a(:)-b))==abc(1,1));
