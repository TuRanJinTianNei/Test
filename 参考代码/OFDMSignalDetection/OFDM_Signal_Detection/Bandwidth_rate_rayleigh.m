function [ B_rate_welch,B_rate_ar ] = Bandwidth_rate_rayleigh(sig1_chnl,fc,fs, snr,itau,power,fmax,itn)
 %不同信道下的带宽检测率
 %瑞利衰落信道
 B_ideal = 8e6;  %OFDM信号带宽理想值
 numb = 1; %蒙特卡洛仿真的次数
 LL = length(snr);
 fprintf('正在计算 Rayleigh 信道带宽检测率...\n');
 fprintf('总进度: 0/%d SNR 值\n', LL);
 fprintf('参数: 每个SNR值重复 %d 次\n', numb);
 fprintf('Rayleigh信道参数: fmax=%.1f Hz, 多径数=%d\n', fmax, length(itau));
 for i = 1:LL
     snr_start = tic;
     fprintf('  处理 SNR = %.1f dB (%d/%d)...\n', snr(i), i, LL);
     for j = 1:numb
         [b_welch(i,j),b_ar(i,j)]=PSD_OFDM_rayleigh(sig1_chnl,fc,fs,snr(i),0,itau,power,fmax,itn);
     end
     B_welch(i) = sum(b_welch(i,:))/numb;
     B_ar(i) = sum(b_ar(i,:))/numb;
     B_rate_welch(i) = 1-abs((B_welch(i)-B_ideal))/B_ideal;
     B_rate_ar(i) = 1-abs((B_ar(i)-B_ideal))/B_ideal;
     snr_time = toc(snr_start);
     fprintf('  完成 SNR = %.1f dB，耗时: %.2f 秒\n', snr(i), snr_time);
 end
 fprintf('Rayleigh 信道计算完成！\n');

