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
 numb = 1; %蒙特卡洛仿真的次数
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
    end
     B_welch(i)=sum(b_welch(i,:))/numb;
     B_ar(i)=sum(b_ar(i,:))/numb;
     B_rate_welch(i)=1-abs((B_welch(i)-B_ideal))/B_ideal;
     B_rate_ar(i)=1-abs((B_ar(i)-B_ideal))/B_ideal;
     
     % 计算估计误差
     error_welch_abs = abs(B_welch(i) - B_ideal);
     error_ar_abs = abs(B_ar(i) - B_ideal);
     error_welch_rel = error_welch_abs / B_ideal * 100;
     error_ar_rel = error_ar_abs / B_ideal * 100;
     
     snr_time = toc(snr_start);
     fprintf('  完成 SNR = %.1f dB，耗时: %.2f 秒\n', snr(i), snr_time);
     fprintf('     Welch算法: 估计值=%.2f MHz, 绝对误差=%.2f MHz, 相对误差=%.2f%%\n', ...
         B_welch(i)/1e6, error_welch_abs/1e6, error_welch_rel);
     fprintf('     AR模型:   估计值=%.2f MHz, 绝对误差=%.2f MHz, 相对误差=%.2f%%\n', ...
         B_ar(i)/1e6, error_ar_abs/1e6, error_ar_rel);
 end
 fprintf('AWGN 信道计算完成！\n');
 
 % 输出估计误差总结
 fprintf('\n========================================\n');
 fprintf('AWGN信道带宽估计误差总结\n');
 fprintf('========================================\n');
 fprintf('理想带宽: %.2f MHz\n', B_ideal/1e6);
 fprintf('%-8s | %-12s | %-15s | %-15s | %-15s | %-15s\n', ...
     'SNR(dB)', 'Welch估计(MHz)', 'Welch绝对误差(MHz)', 'Welch相对误差(%)', ...
     'AR估计(MHz)', 'AR相对误差(%)');
 fprintf('%-8s-+-%-12s-+-%-15s-+-%-15s-+-%-15s-+-%-15s\n', ...
     repmat('-',1,8), repmat('-',1,12), repmat('-',1,15), repmat('-',1,15), ...
     repmat('-',1,15), repmat('-',1,15));
 for i = 1:LL
     error_welch_abs = abs(B_welch(i) - B_ideal);
     error_ar_abs = abs(B_ar(i) - B_ideal);
     error_welch_rel = error_welch_abs / B_ideal * 100;
     error_ar_rel = error_ar_abs / B_ideal * 100;
     fprintf('%-8.1f | %-12.4f | %-15.4f | %-15.2f | %-15.4f | %-15.2f\n', ...
         snr(i), B_welch(i)/1e6, error_welch_abs/1e6, error_welch_rel, ...
         B_ar(i)/1e6, error_ar_rel);
 end
 fprintf('========================================\n\n');
 
figure(1)
plot(snr,B_rate_welch,'b-o');
hold on
plot(snr,B_rate_ar,'k-*');
legend('Welch','AR');
grid on;
xlabel('SNR (dB)')
ylabel('Percentage (%)');
title('AWGN信道带宽检测准确率');
