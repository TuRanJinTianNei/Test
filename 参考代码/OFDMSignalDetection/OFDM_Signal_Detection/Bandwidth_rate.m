function [ B_rate_welch,B_rate_ar ] = Bandwidth_rate(sig1_chnl,fc,fs, snr)
%**************************************************************************
%功能：不同信道下的带宽检测率
%sig1_chnl:输入信号
%fc:载波频率
%fs:采样频率
%snr:信噪比
%B_rate_welch:Welch算法下的带宽检测准确率
%B_rate_ar：AR模型下的带宽检测准确率
%**************************************************************************
 
 B_ideal = 8e6;  %OFDM信号带宽理想值
 numb = 10; %蒙特卡洛仿真的次数
 LL = length(snr);
 B_welch = zeros(1,LL);
 B_ar = zeros(1,LL);
 fprintf('正在计算 AWGN 信道带宽检测率...\n');
 fprintf('总进度: 0/%d SNR 值\n', LL);
 fprintf('参数: 每个SNR值重复 %d 次\n', numb);
 for i=1:LL
     snr_start = tic;
     fprintf('  处理 SNR = %.1f dB (%d/%d)...\n', snr(i), i, LL);
     for j=1:numb
         [b_welch(i,j),b_ar(i,j)]=PSD_OFDM(sig1_chnl,fc,fs,snr(i),0);
         if mod(j,5)==0
             fprintf('    已完成 %d/%d 次重复\n', j, numb);
         end
     end
     B_welch(i)=sum(b_welch(i,:))/numb;
     B_ar(i)=sum(b_ar(i,:))/numb;
     B_rate_welch(i)=1-abs((B_welch(i)-B_ideal))/B_ideal;
     B_rate_ar(i)=1-abs((B_ar(i)-B_ideal))/B_ideal;
     snr_time = toc(snr_start);
     fprintf('  完成 SNR = %.1f dB，耗时: %.2f 秒\n', snr(i), snr_time);
 end
 fprintf('AWGN 信道计算完成！\n');
 
figure(1)
plot(snr,B_rate_welch,'b-o');
hold on
plot(snr,B_rate_ar,'k-*');
legend('Welch','AR');
grid on;
xlabel('SNR (dB)')
ylabel('Percentage (%)');
title('AWGN信道带宽检测准确率');
