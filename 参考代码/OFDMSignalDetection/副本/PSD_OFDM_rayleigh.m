function [B_welch,B_ar] =PSD_OFDM_rayleigh(sig1_chnl,fc,fs,snr,k,itau,power,fmax,itn)
%k=1时绘制PSD
%B的标准值8MHz
sig3_chnl = (sig1_chnl.*exp(1j*2*pi*fc/fs*(0:length(sig1_chnl)-1)));
sig3_chnl = real(MUL_RAYLEIGH(sig3_chnl,itau,power,itn,length(itau),length(sig3_chnl),1/fs,fmax,0));
sig3_chnl = awgn(sig3_chnl,snr,'measured');
%q= Burg(sig1_chnl,'AIC');
%P=mean(abs(sig1_chnl));           %计算信号的标准差，使用函数std来计算，mean和abs函数计算的是不同的标准差
%sig1_chnl=sig1_chnl/P;            %P1和P虽然相同，但一个是除以标准差，一个是除以模型的均值,一个是除以平方和的开方
% [Pxx1,f]=pburg(sig3_chnl,100,4096*2,fs);  %使用AR模型方法进行功率谱估计
[Pxx1, f, p] = Burg(sig3_chnl,fs, 'AIC');   %使用AR模型方法进行功率谱估计
[Pxx2,f1]=pwelch(sig3_chnl,hanning(100),55,4096*2,fs);  %使用welch算法进行功率谱估计

Frc=0:fs/(length(sig3_chnl)):fs/2-1;
OfdmSymComput = 20 * log10(abs(fft(sig3_chnl)));
OfdmSymPSDy = fftshift(OfdmSymComput) - max(OfdmSymComput);

Pxx22 = Pxx2;
Pxx22=Pxx22/min(Pxx22);%找出功率密度中的最小值，归一化处理 
Pxx22=10*log10(Pxx22);%将功率密度转换为dB单位
Pxx22=Pxx22-max(Pxx22);
if k==1
  figure('Name', 'Rayleigh信道-AR模型功率谱密度估计')
  plot(f,Pxx1);
  grid on;
  xlabel('频率 f (Hz)'); 
  ylabel('PSD (dB)');
  title('Rayleigh信道-AR模型方法的功率谱密度估计'); 
  
  figure('Name', 'Rayleigh信道-Welch算法功率谱密度估计')
  plot(f1,Pxx22);
  grid on;
  xlabel('频率 f (Hz)');
  ylabel('PSD (dB)');
  title('Rayleigh信道-Welch算法估计的功率谱密度估计'); 
  
  figure('Name', 'Rayleigh信道-OFDM信号频谱')
  %plot(Frc,OfdmSymPSDy(1,1:end/2));
  plot(Frc,OfdmSymPSDy(1,1:end/2));
  xlabel('频率 f (Hz)');
  ylabel('PSD (dB)');
  title('Rayleigh信道-OFDM信号频谱'); 
end
%****************************
%计算信号的带宽
%****************************

L1=ceil(length(Pxx22)/2);
P1=Pxx22(1:L1,1);
P2=Pxx22(L1:end,1);
[as1, as11]=Proximate(-3,P1);  %取最接近-3dB处信号的f值
band1=f1(as1);
[as2, as22]=Proximate(-3,P2);
band2=f1(as2+L1-1);
B_welch =abs(band1-band2);

L2=ceil(length(Pxx1)/2);
P3=Pxx1(1:L2,1);
P4=Pxx1(L2:end,1);
if snr>4
    [as3, as33]=Proximate(-6,P3);  %取最接近-6dB处信号的f值
    band3=f(as3);
    [as4, as44]=Proximate(-6,P4);
    band4=f(as4+L2-1);
    B_ar =abs(band4-band3);
elseif (snr>0&&snr<=4)
    [as3, as33]=Proximate(-5,P3);  %取最接近-5dB处信号的f值
    band3=f(as3);
    [as4, as44]=Proximate(-5,P4);
    band4=f(as4+L2-1);
    B_ar =abs(band4-band3);
else
    [as3, as33]=Proximate(-3,P3);  %取最接近-3dB处信号的f值
    band3=f(as3);
    [as4, as44]=Proximate(-3,P4);
    band4=f(as4+L2-1);
    B_ar =abs(band4-band3);
end
end