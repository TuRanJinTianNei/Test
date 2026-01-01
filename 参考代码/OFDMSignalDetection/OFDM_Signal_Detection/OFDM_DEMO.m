%**************************************************************************
%OFDM信号识别和子载波检测的功能
%created by Songzhiyong
%2017.03.01
%**************************************************************************
% clc;
% close all;
% clear all;

% 记录仿真开始时间
sim_start_time = tic;
fprintf('========================================\n');
fprintf('OFDM信号识别和子载波检测仿真开始\n');
fprintf('========================================\n');
fprintf('开始时间: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

% %**************************************************************************
% % 使用信号高阶累积量实现信号识别
% % c40_和c42_表示的高阶累积量表示不同信号
% %**************************************************************************
% K = 1;           %当K值为1时绘制图形
N = 20;          %构成一个帧结构的OFDM信号的个数
% para = 128;%设置带内数据子载波数量
% ratio = 1/4 ; %循环前缀比例
% S = 12.5e3  ;   %符号速率12.5Kbps
% snr1 = -5:5:25;   %高斯白噪声的信噪比，表示不同数值的信噪比
snr = -4:2:10;   %高斯白噪声的信噪比
% 
fc = 10e6;   %信号的载波频率
fs = 40e6;   % 采样频率
% 
itau = floor([0,1e-8,2e-8,5e-8,2e-7,5e-7].*fs);  %多径延时
 power = [0,-1.0,-7.0,-10.0,-12.0,-17.0];  %多径信道的每径功率
 fmax = 20 ;   %最大多普勒频率
 itn = [10000,20000,30000,40000,50000,60000];    %瑞利信道的记录次数

% 输出仿真参数
fprintf('\n========== 仿真参数 ==========\n');
fprintf('OFDM符号数量 N = %d\n', N);
fprintf('载波频率 fc = %.2e Hz (%.2f MHz)\n', fc, fc/1e6);
fprintf('采样频率 fs = %.2e Hz (%.2f MHz)\n', fs, fs/1e6);
fprintf('信噪比范围 SNR = [%d:%d:%d] dB\n', snr(1), snr(2)-snr(1), snr(end));
fprintf('最大多普勒频率 fmax = %d Hz\n', fmax);
fprintf('多径延时 itau = [%s] 样本\n', num2str(itau));
fprintf('多径功率 power = [%s] dB\n', num2str(power));
fprintf('================================\n\n');
% 
% %**************************************************************************
% [c42_ofdm,c42_qpsk,c42_qam16,c42_qam64,c42_fsk8] = cumulant(snr1,N,para,ratio,K);
% [VAR_FSK8,VAR_QAM16,VAR_QAM64,VAR_QPSK,VAR_OFDM] = sc_ofdm_wavelet(N,snr1,para,ratio,K);
% [VAR_FSK8_Ray,VAR_QAM16_Ray,VAR_QAM64_Ray,VAR_QPSK_Ray,VAR_OFDM_Ray] = sc_ofdm_wavelet_Ray(N,snr1,para,ratio,K,itau,power,itn,fmax,fs);
% %**************************************************************************
% %�߽������ķ�ʽ��ͬ��������źŵ�ʶ����
% %rate_ofdm��ʾ�źŵ���ȷʶ����
% %**************************************************************************
% [rate_ofdm] = OFDM_rate(snr,N,para,ratio);
% [rate_ofdm1] = OFDM_rate_xiaobo(snr,N,para,ratio);
% [rate_ofdm_Ray] = OFDM_rate_xiaobo_Ray(snr,N,para,ratio,itau,power,itn,fmax,fs);
% %**************************************************************************
% %OFDM不同子载波数目对高阶累积量估计值的不确定性影响
% %c42_ofdm64,c42_ofdm128,c42_ofdm256分别表示
% %子载波数目为128,256,512时的累积量
% %**************************************************************************
% [V_ofdm128,V_ofdm256,V_ofdm512] = ofdm_dif_para(snr,N,ratio,K);
%  
% %**************************************************************************
% %ofdm�ز�Ƶ�ʲ�������
% %ofdmѭ���׹���ofdm�ź��ز�Ƶ��
% %**************************************************************************
% trst_rate = 5e5;  % �źŷ�������,����Ƭ����
% M = 16; % ƽ������ 
% sig = ofdm(N,para,ratio);
% Ns = length(sig);       % ѭ���׼��������ȣ�����С�ڵ����ź����г���
%  
% %**************************************************************************
% %��������OFDM�źŲ���ѭ��ƽ����
% %**************************************************************************
% s_n = ceil(fs/trst_rate); % �������ʽ���ΪOFDM�ź����ʵ������� 
% sign = sig(ones(s_n,1),:); % ��Ԫ����
% sign = reshape(sign, 1, s_n*length(sign));
% data = sign.*exp(1i*2*pi*fc/fs*(0:length(sign)-1));
% %�߰��ŵ���ѭ���׹�����Ƶ
% data_awgn = awgn(data,25,'measured');
% [f_awgn] = cyclic_spectrum(real(data_awgn), Ns, fs, M,K);  % ѭ������Ƶ��f
% 
% %�����ŵ���ѭ���׹���
% %data_rayleigh = RayFade(itau,power,fmax,fs,data,25);
% data_rayleigh = (MUL_RAYLEIGH(data,itau,power,itn,length(itau),length(data),1/fs,fmax,0));
% data_rayleigh = data_rayleigh/std(data_rayleigh);
% data_rayleigh = awgn(data_rayleigh,25,'measured');
% [f_rayleigh] = cyclic_spectrum(real(data_rayleigh), Ns, fs, M,K);  % ѭ������Ƶ��f
%**************************************************************************
%主程序
%**************************************************************************
fprintf('开始生成 OFDM 信号...\n');
signal_gen_start = tic;
[sig_ofdm] = PSD_generate(N);
signal_gen_time = toc(signal_gen_start);
fprintf('OFDM 信号生成完成！耗时: %.2f 秒\n\n', signal_gen_time);

%[band_welch,band_ar] = PSD_OFDM(sig_ofdm,fc,fs,20,K);  %检测功率谱(25db)
% %**************************************************************************
% % AWGN信道下的带宽检测率计算（已注释，不再单独仿真AWGN情况）
% %**************************************************************************
% fprintf('开始计算 AWGN 信道带宽检测率...\n');
% awgn_start = tic;
% [B_rate_welch,B_rate_ar ] = Bandwidth_rate( sig_ofdm,fc,fs,snr);   %计算信号的带宽检测率
% awgn_time = toc(awgn_start);
% fprintf('AWGN 信道计算完成！总耗时: %.2f 秒\n\n', awgn_time);

%**************************************************************************
% Rayleigh信道下的带宽检测率计算
%**************************************************************************
fprintf('开始计算 Rayleigh 信道带宽检测率...\n');
rayleigh_start = tic;
[B_rate_welch_rayleigh,B_rate_ar_rayleigh ] = Bandwidth_rate_rayleigh(sig_ofdm,fc,fs,snr,itau,power,fmax,itn);   %计算信号的带宽检测率
rayleigh_time = toc(rayleigh_start);
fprintf('Rayleigh 信道计算完成！总耗时: %.2f 秒\n\n', rayleigh_time);

fprintf('正在生成 Rayleigh 信道 PSD 图形...\n');
psd_start = tic;
snr_plot = 20;  % 用于生成 PSD 图形的 SNR 值（标量）
[B_welch,B_ar] =PSD_OFDM_rayleigh(sig_ofdm,fc,fs,snr_plot,1,itau,power,fmax,itn);
psd_time = toc(psd_start);
fprintf('PSD 图形生成完成！耗时: %.2f 秒\n\n', psd_time);

fprintf('所有计算完成！正在绘制结果...\n');
plot_start = tic;
figure('Name', 'Rayleigh信道下的带宽检测率对比')
plot(snr,B_rate_ar_rayleigh,'r-o');
hold on
plot(snr,B_rate_welch_rayleigh,'b--s');
xlabel('SNR (dB)');
ylabel('Percentage (%)');
legend('AR模型法','Welch算法', 'Location', 'best');
title('Rayleigh信道下的带宽检测率');
grid on;
plot_time = toc(plot_start);
fprintf('图形绘制完成！耗时: %.2f 秒\n\n', plot_time);

% 输出总仿真时间
total_time = toc(sim_start_time);
fprintf('========================================\n');
fprintf('仿真完成！\n');
fprintf('结束时间: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf('总仿真时间: %.2f 秒 (%.2f 分钟)\n', total_time, total_time/60);
fprintf('========================================\n');

%**************************************************************************
%AWGN信道下循环自相关估计有效数据长度和总长度以及循环前缀长度
%**************************************************************************
% sig_a = sig.*exp(1i*2*pi*fc/fs*(0:length(sig)-1));  %上变频到载频 
% [Tu_128,Ts_128,Tg_128 ] = effectivelength(sig_a,fs,20,N,K);
% [ Tu_rate_128,Ts_rate_128,Tg_rate_128] = Length_rate(sig_a,fs,snr,N,1 );
% 
% %**************************************************************************
% %Rayleigh信道下估计有效数据长度和总长度以及循环前缀长度
% %**************************************************************************
% [Tu_128_rayleigh,Ts_128_rayleigh,Tg_128_rayleigh ] = effectivelength_rayleigh(sig_a,fs,10,N,1,itau,power,fmax,itn);
% [Tu_rate_128_rayleigh,Ts_rate_128_rayleigh,Tg_rate_128_rayleigh] = Length_rate_rayleigh(sig_a,fs,snr,N,0,itau,power,fmax,itn );
% figure
% plot(snr,Tu_rate_128_rayleigh,'r-o');
% hold on

% plot(snr,Tu_rate_128,'k-x');
% xlabel('snr/db');
% ylabel('percentage/%');
% legend('Rayleigh','Awgn');
% title('不同信道对有效数据长度估计的影响')
% figure
% plot(snr,Ts_rate_128_rayleigh,'r-o');
% hold on
% plot(snr,Ts_rate_128,'k-x');
% legend('Rayleigh','Awgn');
% xlabel('snr/db');
% ylabel('percentage/%');
% title('不同信道对符号总长度估计的影响')
% figure
% plot(snr,Tg_rate_128_rayleigh,'r-o');
% hold on
% plot(snr,Tg_rate_128,'k-x');
% legend('Rayleigh','Awgn');
% xlabel('snr/db');
% ylabel('percentage/%');
% title('不同信道对循环前缀长度估计的影响')
% 
% %**************************************************************************
% %AWGN信道下子载波数目对信号参数估计性能的影响
% %**************************************************************************
% % sig_256 = ofdm(N,256,ratio);
% % sig256 = sig_256.*exp(1i*2*pi*fc/fs*(0:length(sig_256)-1));  %上变频到载频 
% % [ Tu_rate_256,Ts_rate_256,Tg_rate_256] = Length_rate_difpara(sig256,fs,snr,N,0 );
% % figure
% % plot(snr,Tu_rate_256,'r-o');
% % hold on
% % plot(snr,Tu_rate_128,'k-x');
% % xlabel('snr/db');
% % ylabel('percentage/%');
% % legend('256个子载波','128个子载波');
% % title('子载波数目对有效数据长度估计的影响')
% % figure
% % plot(snr,Ts_rate_256,'r-o');
% % hold on
% % plot(snr,Ts_rate_128,'k-x');
% % legend('256个子载波','128个子载波');
% % xlabel('snr/db');
% % ylabel('percentage/%');
% % title('子载波数目对符号总长度估计的影响')
% % figure
% % plot(snr,Tg_rate_256,'r-o');
% % hold on
% % plot(snr,Tg_rate_128,'k-x');
% % legend('256个子载波','128个子载波');
% % xlabel('snr/db');
% % ylabel('percentage/%');
% % title('子载波数目对循环前缀长度估计的影响')
% 
% %**************************************************************************
% %AWGN信道下符号个数对信号参数估计性能的影响
% %**************************************************************************
% sig_10 = ofdm(10,para,ratio);
% sig10 = sig_10.*exp(1i*2*pi*fc/fs*(0:length(sig_10)-1));  %上变频到载频 
% [ Tu_rate_10,Ts_rate_10,Tg_rate_10] = Length_rate(sig10,fs,snr,10,0 );
% sig_30 = ofdm(30,para,ratio);
% sig30 = sig_30.*exp(1i*2*pi*fc/fs*(0:length(sig_30)-1));  %上变频到载频 
% [Tu_rate_30,Ts_rate_30,Tg_rate_30] = Length_rate(sig30,fs,snr,30,0 );
% figure
% plot(snr,Tu_rate_10,'r-o');
% hold on
% plot(snr,Tu_rate_128,'k-x');
% hold on
% plot(snr,Tu_rate_30,'b-*');
% legend('10个符号','20个符号','30个符号');
% xlabel('snr/db');
% ylabel('percentage/%');
% title('符号个数对有效数据长度估计的影响')
% figure
% plot(snr,Ts_rate_10,'r-o');
% hold on
% plot(snr,Ts_rate_128,'k-x');
% hold on
% plot(snr,Ts_rate_30,'b-*');
% legend('10个符号','20个符号','30个符号');
% xlabel('snr/db');
% ylabel('percentage/%');
% title('符号个数对符号总长度估计的影响')
% figure
% plot(snr,Tg_rate_10,'r-o');
% hold on
% plot(snr,Tg_rate_128,'k-x');
% hold on
% plot(snr,Tg_rate_30,'b-*');
% legend('10个符号','20个符号','30个符号');
% xlabel('snr/db');
% ylabel('percentage/%');
% title('符号个数对循环前缀长度估计的影响')
% 
% %**************************************************************************
% %AWGN信道下循环前缀比例对信号参数估计性能的影响
% %**************************************************************************
% sig_ratio = ofdm(N,para,1/8);
% sig_r = sig_ratio.*exp(1i*2*pi*fc/fs*(0:length(sig_ratio)-1));  %上变频到载频 
% [ Tu_rate_ratio,Ts_rate_ratio,Tg_rate_ratio] = Length_rate_difratio(sig_r,fs,snr,N,0 );
% 
% sig_ratio1 = ofdm(N,para,3/16);
% sig_r1= sig_ratio1.*exp(1i*2*pi*fc/fs*(0:length(sig_ratio1)-1));  %上变频到载频 
% [Tu_rate_ratio1,Ts_rate_ratio1,Tg_rate_ratio1] = Length_rate_difratio(sig_r1,fs,snr,N,0 );
% figure
% plot(snr,Tu_rate_ratio,'r-o');
% hold on
% plot(snr,Tu_rate_ratio1,'b-*');
% hold on
% plot(snr,Tu_rate_128,'k-x');
% legend('ratio=1/8','ratio=3/16','ratio=1/4');
% xlabel('snr/db');
% ylabel('percentage/%');
% title('循环前缀比例对有效数据长度估计的影响')
% figure
% plot(snr,Ts_rate_ratio,'r-o');
% hold on
% plot(snr,Ts_rate_ratio1,'b-*');
% hold on
% plot(snr,Ts_rate_128,'k-x');
% legend('ratio=1/8','ratio=3/16','ratio=1/4');
% xlabel('snr/db');
% ylabel('percentage/%');
% title('循环前缀比例对符号总长度估计的影响')
% figure
% plot(snr,Tg_rate_ratio,'r-o');
% hold on
% plot(snr,Tg_rate_ratio1,'b-*');
% hold on
% plot(snr,Tg_rate_128,'k-x');
% legend('ratio=1/8','ratio=3/16','ratio=1/4');
% xlabel('snr/db');
% ylabel('percentage/%');
% title('循环前缀比例对循环前缀长度估计的影响')
% 
% %**************************************************************************
% %AWGN信道下OFDM子载波数目估计
% %**************************************************************************
% [carry_num] = carrier_number(sig_a ,N,25,K);
% [Carry_num_rate_128 ] = Rate_Carrynum(sig_a,snr,para,N,K);
% 
% %**************************************************************************
% %AWGN信道下每帧的符号个数对子载波数目估计的影响
% %**************************************************************************
% [Carry_num_rate_10 ] = Rate_Carrynum(sig_10,snr,para,10,0);
% %[Carry_num_rate_30 ] = Rate_Carrynum(sig_30,snr,para,30,0);
% Carry_num_rate_30 = Tu_rate_30 ;
% figure
% plot(snr,Carry_num_rate_10,'r-o');
% hold on
% plot(snr,Carry_num_rate_128,'k-x');
% hold on
% plot(snr,Carry_num_rate_30,'b-*');
% xlabel('snr/db');
% ylabel('percentage/%');
% legend('10个符号','20个符号','30个符号');
% title('符号个数对子载波数目估计的影响')
% 
% %**************************************************************************
% %AWGN信道下循环前缀的比例对子载波数目估计的影响
% %**************************************************************************
% [Carry_num_rate_ratio ] = Rate_Carrynum(sig_ratio,snr,para,N,0);
% [Carry_num_rate_ratio1 ] = Rate_Carrynum(sig_ratio1,snr,para,N,0);
% figure
% plot(snr,Carry_num_rate_ratio,'r-o');
% hold on
% plot(snr,Carry_num_rate_ratio1,'b-*');
% hold on
% plot(snr,Carry_num_rate_128,'k-x');
% xlabel('snr/db');
% ylabel('percentage/%');
% legend('ratio=1/8','ratio=3/16','ratio=1/4');
% title('循环前缀比例对子载波数目估计的影响')
% 
% %**************************************************************************
% %AWGN信道下子载波数目对子载波数目估计的影响
% %**************************************************************************
% %[Carry_num_rate_256 ] = Rate_Carrynum(sig_256,snr,256,N,0);
% % Carry_num_rate_256 = Tu_rate_256;
% % figure
% % plot(snr,Carry_num_rate_256,'r-o');
% % hold on
% % plot(snr,Carry_num_rate_128,'k-x');
% % legend('256个子载波','128个子载波');
% % xlabel('snr/db');
% % ylabel('percentage/%');
% % title('子载波数目对子载波数目估计的影响')
% 
% %**************************************************************************
% %Rayleigh信道下OFDM子载波数目估计准确率
% %**************************************************************************
% [Carry_num_rate_rayleigh ]= Rate_Carrynum_rayleigh(sig_a, snr,N ,para,itau,power,fmax,fs,itn);
% figure
% plot(snr,Carry_num_rate_rayleigh,'r-o');
% hold on
% plot(snr,Carry_num_rate_128,'k-x');
% legend('Rayleigh','Awgn');
% xlabel('snr/db');
% ylabel('percentage/%');
% title('不同信道对子载波数目估计准确率');
% 
% %**************************************************************************
% %利用带宽估计子载波数量
% %**************************************************************************
% [BB, Carrynum_B] =solve_carrynum_rate(sig_ofdm,fc,fs,snr,N );    %利用带宽的方式估计子载波数量
