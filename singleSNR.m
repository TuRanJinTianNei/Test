% 对OFDM系统进行建模，采用16QAM星座映射
tic
clc;
clear all;
close all;
%============= OFDM系统参数设置================

carrier_count=200;
symbols_per_carrier=12;
bits_per_symbol=4;
IFFT_bin_length=512;
PrefixRatio=1/4;
GI=PrefixRatio*IFFT_bin_length;  % 保护前缀长度
beta=1/16; %滚降系数
GIP=beta*(IFFT_bin_length+GI);  %加窗扩展时域信号长度
% 只在单一SNR下仿真（修改下行数值即可）
targetSNRdB = 15;

%=============发送端=======================
%=============信号产生=====================

baseband_out_length=carrier_count*symbols_per_carrier*bits_per_symbol;
%%进行子载波共轭映射，使得OFDM符号经过IFFT之后是实信号
carriers=(1:carrier_count)+(floor(IFFT_bin_length/4)-floor(carrier_count/2));
conjugate_carriers=IFFT_bin_length-carriers+2; 
rand('twister',0);
baseband_out=round(rand(1,baseband_out_length));

%===============16QAM调制=============

complex_carrier_matrix=qam16(baseband_out);
complex_carrier_matrix=reshape(complex_carrier_matrix',carrier_count,symbols_per_carrier)';

%================IFFT==================

IFFT_modulation=zeros(symbols_per_carrier,IFFT_bin_length);
IFFT_modulation(:,carriers)=complex_carrier_matrix;
IFFT_modulation(:,conjugate_carriers)=conj(complex_carrier_matrix);

%================画出其中一个OFDM符号的幅度和相位==============================
% Figure 1: IFFT各频点的幅度（展示激活的子载波与其共轭镜像的幅度分布）
% - 可直观看到载波分配与未用频点（幅度为0），理解频域输入结构
figure(1);
stem(0:IFFT_bin_length-1, abs(IFFT_modulation(2,1:IFFT_bin_length)),'b*-')
grid on
axis ([0 IFFT_bin_length -0.5 4.5]);
ylabel('Magnitude');
xlabel('IFFT Bin');
title('OFDM Carrier Frequency Magnitude');
 
% Figure 2: IFFT各频点的相位（度）
% - 体现共轭映射带来的相位对称性，验证实值时域信号条件
figure(2);
plot(0:IFFT_bin_length-1, (180/pi)*angle(IFFT_modulation(2,1:IFFT_bin_length)), 'go')
hold on
stem(0:carriers-1, (180/pi)*angle(IFFT_modulation(2,1:carriers)),'b*-');
stem(0:conjugate_carriers-1, (180/pi)*angle(IFFT_modulation(2,1:conjugate_carriers)),'b*-');
axis ([0 IFFT_bin_length -200 +200])
grid on
ylabel('Phase (degrees)')
xlabel('IFFT Bin')
title('OFDM Carrier Phase')

%=================未添加保护间隔的OFDM符号===================================
signal_after_IFFT=ifft(IFFT_modulation,IFFT_bin_length,2);
time_wave_matrix=signal_after_IFFT;
% Figure 3: 单个OFDM符号（未加循环前缀/后缀）的时域波形
% - 展示IFFT输出的一个符号周期包络与幅度范围
figure(3)
plot(0:IFFT_bin_length-1,time_wave_matrix(2,:));
axis([0 IFFT_bin_length -0.2 0.2]);
grid on
ylabel('Amplitude');
xlabel('Time');
title('OFDM Time Signal, One Symbol Period');

%================添加循环前缀与后缀==================
XX=zeros(symbols_per_carrier,IFFT_bin_length+GI+GIP);
for k=1:symbols_per_carrier
    for i=1:IFFT_bin_length
        XX(k,i+GI)=signal_after_IFFT(k,i);
    end
    for i=1:GI
        XX(k,i)=signal_after_IFFT(k,i+IFFT_bin_length-GI);
    end
    for j=1:GIP
         XX(k,IFFT_bin_length+GI+j)=signal_after_IFFT(k,j);
    end
end  
time_wave_matrix_cp=XX;
% Figure 4: 单个OFDM符号添加循环前缀(CP)与后缀后的时域波形
% - 左侧CP、右端后缀（用于窗口过渡），观察边沿延拓
figure(4);
plot(0:length(time_wave_matrix_cp)-1,time_wave_matrix_cp(2,:));
axis([0, length(time_wave_matrix_cp), -0.2, 0.2]);
grid on;
ylabel('Amplitude');
xlabel('Time');
title('OFDM Time Signal with CP, One Symbol Period');

%==============OFDM符号加窗=================================================
windowed_time_wave_matrix_cp=zeros(1,IFFT_bin_length+GI+GIP);
for i = 1:symbols_per_carrier
    windowed_time_wave_matrix_cp(i,:)=real(time_wave_matrix_cp(i,:)).*rcoswindow(beta,IFFT_bin_length+GI)';
end
% Figure 5: 加窗后的单个OFDM符号时域波形（含CP/后缀）
% - 余弦窗平滑边沿，降低频谱旁瓣
figure(5)
plot(0:IFFT_bin_length-1+GI+GIP,windowed_time_wave_matrix_cp(2,:)); 
axis([0, IFFT_bin_length-1+GI+GIP, -0.2, 0.2]);
grid on;
ylabel('Amplitude');
xlabel('Time');
title('OFDM Time Signal Apply a Window , One Symbol Period');

%========================生成发送信号，并串变换====================================
%这样做的目的在于将加窗的OFDM符号经过并串转换之后时，对于循环前缀是每个子信道都需要添加
%但是对于循环后缀不需要每个子信道都添加，因此串并之后只需要对整个串行符号添加就行
windowed_Tx_data=zeros(1,symbols_per_carrier*(IFFT_bin_length+GI)+GIP);
windowed_Tx_data(1:IFFT_bin_length+GI+GIP)=windowed_time_wave_matrix_cp(1,:);
for i = 1:symbols_per_carrier-1 ;
    windowed_Tx_data((IFFT_bin_length+GI)*i+1:(IFFT_bin_length+GI)*(i+1)+GIP)=windowed_time_wave_matrix_cp(i+1,:);
 end
%================================================================================
Tx_data_withoutwindow =reshape(time_wave_matrix_cp',(symbols_per_carrier)*(IFFT_bin_length+GI+GIP),1)';%未添加窗信号
Tx_data =reshape(windowed_time_wave_matrix_cp',(symbols_per_carrier)*(IFFT_bin_length+GI+GIP),1)';%加窗发送信号
 %===============================================================================
temp_time1 = (symbols_per_carrier)*(IFFT_bin_length+GI+GIP); 
% Figure 6: 发送端串行时域信号（两种串接方式对比）
% - 子图(1) Tx_data: 将每个符号(含CP/后缀)直接串接
% - 子图(2) windowed_Tx_data: 仅在整帧末尾添加一次后缀以优化过渡
figure (6)
subplot(2,1,1);
plot(0:temp_time1-1,Tx_data );
grid on
ylabel('Amplitude (volts)')
xlabel('Time (samples)')
title('OFDM Time Signal')
temp_time2 =symbols_per_carrier*(IFFT_bin_length+GI)+GIP;
subplot(2,1,2);
plot(0:temp_time2-1,windowed_Tx_data);
grid on
ylabel('Amplitude (volts)')
xlabel('Time (samples)')
title('OFDM Time Signal')
%=================未加窗发送信号频谱===============================================
symbols_per_average = ceil(symbols_per_carrier/5);
 avg_temp_time = (IFFT_bin_length+GI+GIP)*symbols_per_average;
averages = floor(temp_time1/avg_temp_time);
average_fft(1:avg_temp_time) = 0;
for a = 0:(averages-1)
 subset_ofdm = Tx_data_withoutwindow (((a*avg_temp_time)+1):((a+1)*avg_temp_time));
 subset_ofdm_f = abs(fft(subset_ofdm));
 average_fft = average_fft + (subset_ofdm_f/averages);
end
average_fft_log = 20*log10(average_fft);
% Figure 7: 未加窗串行发送信号的平均功率谱（dB），多段平均降低方差
figure (7)

plot((0:(avg_temp_time-1))/avg_temp_time, average_fft_log)
hold on
plot(0:1/IFFT_bin_length:1, -35, 'rd')
grid on
axis([0 0.5 -40 max(average_fft_log)])
ylabel('Magnitude (dB)')
xlabel('Normalized Frequency (0.5 = fs/2)')
title('OFDM Signal Spectrum without windowing')
%===============加窗的发送信号频谱=================================================
symbols_per_average = ceil(symbols_per_carrier/5);
avg_temp_time = (IFFT_bin_length+GI+GIP)*symbols_per_average;
averages = floor(temp_time1/avg_temp_time);
average_fft(1:avg_temp_time) = 0;
for a = 0:(averages-1)
 subset_ofdm = Tx_data(((a*avg_temp_time)+1):((a+1)*avg_temp_time));
 subset_ofdm_f = abs(fft(subset_ofdm));
 average_fft = average_fft + (subset_ofdm_f/averages);
end
average_fft_log = 20*log10(average_fft);
% Figure 8: 加窗后的串行发送信号平均功率谱（dB），旁瓣应较Figure 7更低
figure(8)
plot((0:(avg_temp_time-1))/avg_temp_time, average_fft_log) ;
hold on
plot(0:1/IFFT_bin_length:1, -35, 'rd')
grid on
axis([0 0.5 -40 max(average_fft_log)])
ylabel('Magnitude (dB)')
xlabel('Normalized Frequency (0.5 = fs/2)')
title('Windowed OFDM Signal Spectrum')

%===============添加噪声(AWGN信道)==================
Tx_signal_power=var(windowed_Tx_data);
linear_SNR=10^(targetSNRdB/10);
noise_sigma=Tx_signal_power/linear_SNR;
noise_scale_factor=sqrt(noise_sigma);
noise=randn(1,((symbols_per_carrier)*(IFFT_bin_length+GI)+GIP))*noise_scale_factor;
Rx_data=windowed_Tx_data+noise;

%==============接收端进行串并转换 去除循环前缀和后缀（对加窗信号进行处理）===========
Rx_data_matrix=zeros(symbols_per_carrier,IFFT_bin_length+GI+GIP);
for i=1:symbols_per_carrier
    Rx_data_matrix(i,:)=Rx_data(1,(i-1)*(IFFT_bin_length+GI)+1:i*(IFFT_bin_length+GI)+GIP);
end
Rx_data_complex_matrix=Rx_data_matrix(:,GI+1:IFFT_bin_length+GI); 
%===============================================================================

%====================OFDM解调， 16QAM解星座映射==================================
Y1=fft(Rx_data_complex_matrix,IFFT_bin_length,2);
Rx_carriers=Y1(:,carriers);
Rx_mag=abs(Rx_carriers);
Rx_phase=angle(Rx_carriers);
%============================================================================ 
[M, N]=pol2cart(Rx_phase, Rx_mag); 
Rx_complex_carrier_matrix = complex(M, N);
% Figure 9: 接收端子载波星座图（IQ平面），SNR越高散点越集中于16QAM理想点
figure(9);
plot(Rx_complex_carrier_matrix,'*b');
axis([-4, 4, -4, 4]);
title('the constellation of the received signal in XY coordinates ')
grid on

%=============16QAM解调========================================================
Rx_serial_complex_symbols = reshape(Rx_complex_carrier_matrix',size(Rx_complex_carrier_matrix, 1)*size(Rx_complex_carrier_matrix,2),1)' ; 
Rx_decoded_binary_symbols=demoduqam16(Rx_serial_complex_symbols);
%================================================================================
baseband_in = Rx_decoded_binary_symbols;
% Figure 10: 比特流对比（前100比特）——上:发送，下:接收判决
figure(10);
subplot(2,1,1);
stem(baseband_out(1:100));
subplot(2,1,2);
stem(baseband_in(1:100));
%================误码率计算=======================================================
bit_errors=find(baseband_in ~=baseband_out);
bit_error_count = size(bit_errors, 2); 
ber=bit_error_count/baseband_out_length
%================扩展：命令行输出关键指标（与 test1.m 一致的精简版）================
null_subcarriers = IFFT_bin_length - 2*carrier_count; % 包含DC/保护带

fprintf('\n==== Simulation Summary ====\n');
fprintf('SNR target        : %.2f dB\n', targetSNRdB);
fprintf('BER               : %.6g\n', ber);
fprintf('Tx power (var)    : %.6g\n', Tx_signal_power);
fprintf('Noise variance    : %.6g\n', noise_sigma);
fprintf('IFFT length (N)   : %d\n', IFFT_bin_length);
fprintf('Active carriers   : %d\n', carrier_count);
fprintf('Null carriers     : %d\n', null_subcarriers);
fprintf('CP length (GI)    : %d (ratio=%.3f)\n', GI, GI/IFFT_bin_length);
fprintf('Suffix length     : %d\n', GIP);
fprintf('One OFDM symbol   : %d samples (N+GI+GIP)\n', IFFT_bin_length+GI+GIP);
fprintf('Window roll-off   : 1/%d\n', round(1/beta));
fprintf('Symbols per frame : %d\n', symbols_per_carrier);
fprintf('==============================\n\n');
% Figure 11: BER-SNR曲线（逐SNR点累积绘制），纵轴对数刻度
figure(11);
semilogy(targetSNRdB,ber,'r*');
xlabel('SNR(dB)');
ylabel('BER');
title('the BER  performance of the OFDM system');
grid on
toc
%================end of file====================











