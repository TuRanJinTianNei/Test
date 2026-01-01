function [BB_AR, Carrynum_B ] = solve_carrynum_rate( signal,fc,fs,snr,N )
%**************************************************************************
%功能：利用带宽的方式估计子载波数量
%signal：输入信号
%fc：载波频率
%fs：采样频率
%snr：信噪比
%N：符号个数
%**************************************************************************
for i = 1:length(snr)
 [BB_W(i),BB_AR(i)] = PSD_OFDM(signal,fc,fs,snr(i),0);  %检测功率谱(25db)
 sig_22 = awgn(signal,snr(i),'measured');
 [Txu(i)] = overfind_num(sig_22,length(sig_22)/N,N,5);
 Vu(i)= Txu(i)*(1/fs);
 Carrynum_B(i) =round(BB_AR(i)*Vu(i));
end
%carrynum_rate =abs(round(Carrynum_B)-128)/128;
figure
plot(snr,Carrynum_B,'b-o');
hold on
plot([snr(1) snr(end)],[128 128],'k');
 legend('估计值','真实值');
xlabel('SNR (dB)');
ylabel('子载波数量');
title('子载波数量估计值');