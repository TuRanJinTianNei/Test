function c21=rt_C21(rt)%计算累积量C21
crt=conj(rt);%rt的共轭
c21=mean(rt.*crt);%方差