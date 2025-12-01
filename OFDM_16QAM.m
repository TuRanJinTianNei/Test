% ============================================================================
% �ļ���: test1.m
% ����: OFDMϵͳ�����������������
% ����: 
%   1. ʵ��������OFDMϵͳ��·�����Ͷˣ�16QAM���ơ�IFFT��CP���Ӵ����� 
%      �ŵ���AWGN���� ���նˣ�ȥCP��FFT��16QAM�����
%   2. ����13��Figure���п��ӻ�������������
%      - �����ŷ�������ӻ���Figure 1-5��
%      - ֡���������ź���Ƶ�׷�����Figure 6-8��
%      - �ŵ�Ӱ������ն�������Figure 12-13��
%      - ����о������ܣ�Figure 9-11������BER-SNR���ߣ�
%   3. ����BER-SNR���ߣ�10-30dB������2dB��
%   4. ��15dB SNR�½�����ϸ����������SNR������BER����
% ============================================================================

tic
clc;
clear all;
close all;
%===================== Figure ����˵�� =====================
%
% ��һ�������ŷ�������ӻ�
%   Figure 1 : IFFT ��Ƶ�����
%   Figure 2 : IFFT ��Ƶ����λ
%   Figure 3 : ������ʱ��δ�� CP/��׺��
%   Figure 4 : ������ʱ�򣨼� CP ���׺��
%   Figure 5 : ������ʱ�򣨼Ӵ���
%
% ������֡���������ź���Ƶ�׷���
%   Figure 6 : ���з���ʱ��Աȣ��ϣ�����ţ��£�֡��+ĩβ��׺��
%   Figure 7 : �Ӵ�ǰ��ֲ�����ʱ��Աȣ��� Figure 12 ��ͬ��Χ��
%   Figure 8 : �Ӵ�ǰ��ֲ�����Ƶ��Աȣ��� Figure 12 ��ͬ��Χ��
%   ע��Figure 7 չʾ�Ӵ�ǰ���ź��ھֲ������ڵ�ʱ���Σ����ڹ۲촰�������źű�Ե��ƽ�����á�
%       Figure 8 չʾ�Ӵ�ǰ���ź��ھֲ������ڵķ����ף�dBֵ�������ڶԱȼӴ����԰������Ч����
%       ���߾�ʹ���� Figure 12 ��ͬ�ľֲ����ڽ��з������������������ο���
%
% �������ŵ�Ӱ������ն�����
%   Figure 12: �ֲ�����ʱ��Աȣ���SNR/MSE�����д�ӡ��
%   Figure 13: �ֲ�����˫���׶Ա�
%
% ���ġ�����о�������
%   Figure 9 : ���ն�16QAM����ͼ
%   Figure 10: ��������/��ǰ100λ�Ա�
%   Figure 11: BER-SNR���ߣ�10-30dB������2dB��
%
%==========================================================
%===============================================================================
% ��ϵͳ�������á�
%===============================================================================
% ˵�������²���������OFDMϵͳ�ĺ������ã��ɸ���ʵ���������

carrier_count=200;                 % ��Ч�������ز�������������������ز���
symbols_per_frame=50;              % ÿ֡ OFDM ������
total_symbols=1000;                % �ܹ�Ҫ����� OFDM ������
symbols_per_carrier=total_symbols; % Ϊ���ݺ�������ߴ磬����ԭ������ʾ"�ܷ�����"
bits_per_symbol=4;                 % ÿ�����ز����صı�������4=16QAM��

% IFFT/FFT����
IFFT_bin_length=512;               % IFFT ������������Ч���ա�������DC��

% �����������
PrefixRatio=1/4;                   % ѭ��ǰ׺������GI/N��������1/4��ʾGI=N/4
GI=PrefixRatio*IFFT_bin_length;    % ����ǰ׺���ȣ���������

% �Ӵ�������LTE���
% ============================================================================
% LTE������·OFDMʱ����żӴ�ԭ����
% 1. Ŀ�ģ����Ʒ��ű�Ե�Ķ������䣬���ٴ�����䣨Ƶ��й¶��
% 2. �ص����ƣ�����ƴ��ʱ���ڿɿص��ص����ص�����������ѭ��ǰ׺��CP����
% 3. ��Ӱ����Ч���ݣ��ص�ֻ��CP�ڣ����ն�ȥ��CPʱ����Ӱ��
% ============================================================================
beta=1/16;                         % �Ӵ�����ϵ�������ɴ�ռ�ȣ�ԽС����Խ�̣�
GIP=beta*(IFFT_bin_length+GI);     % �Ҷ˺�׺���ȣ���ϴ���β������
GIP=min(GIP, GI);                  % LTEҪ�󣺺�׺���Ȳ�����CP���ȣ�ȷ���ص�������CP��
GIP=floor(GIP);                    % ȷ��Ϊ����

% ���Ž��紦����������������LTE��׼������
% ============================================================================
% LTE������ʹ�������Ҵ�ƽ�����ű�Ե������ƴ��ʱ��CP��Χ���ص����
% - ÿ�����Žṹ��[CP(GI) | ����(N) | ��׺(GIP)]
% - �ص����򣺵�ǰ���ŵĺ�׺(GIP)����һ�����ŵ�CPǰGIP�������ص�
% - �ص���ӣ����ص������������ŵķ�����ӣ�ȷ��ƽ������
% - ���նˣ�ȥ��CPʱֻȡ���岿�֣������ص�Ӱ��
% ============================================================================
solution_method = 3;               % ʹ��LTE��׼�������ص���ӣ�����3��

% �ŵ�����
targetSNRdB = 15;                  % Ŀ������ȣ�dB��������15dB�µ���ϸ����

%===============================================================================
% �����Ͷ˴������̡�
%===============================================================================

%------------------------------------------------------------------------------
% ����1: �������������
%------------------------------------------------------------------------------
baseband_out_length=carrier_count*symbols_per_carrier*bits_per_symbol;
% ʹ�����������������ȷ��ÿ�����в�ͬ�����踴��ʵ��ɸ�Ϊ rng(fixed_seed)
rng('default'); rng('shuffle');
baseband_out = randi([0 1], 1, baseband_out_length);

%------------------------------------------------------------------------------
% ����2: 16QAM����
%------------------------------------------------------------------------------
complex_carrier_matrix=qam16(baseband_out);
complex_carrier_matrix=reshape(complex_carrier_matrix',carrier_count,symbols_per_carrier)';

%------------------------------------------------------------------------------
% ����3: Ƶ�����ز�ӳ�䣨���찣�����ع���Գƣ�ʹIFFT���Ϊʵ�źţ�
%------------------------------------------------------------------------------
% ���ز��������㣺�������ز�����ӳ�䣬ʹ��OFDM���ž���IFFT֮����ʵ�ź�
carriers=(1:carrier_count)+(floor(IFFT_bin_length/4)-floor(carrier_count/2));
conjugate_carriers=IFFT_bin_length-carriers+2; 

IFFT_modulation=zeros(symbols_per_carrier,IFFT_bin_length);
IFFT_modulation(:,carriers)=complex_carrier_matrix;
IFFT_modulation(:,conjugate_carriers)=conj(complex_carrier_matrix);

%------------------------------------------------------------------------------
% ����4: IFFT�任��Ƶ�� �� ʱ��
%------------------------------------------------------------------------------
% Figure 1: IFFT��Ƶ��ķ��ȣ�չʾ��������ز����乲���ķ��ȷֲ���
% - ��ֱ�ۿ����ز�������δ��Ƶ�㣨����Ϊ0��������Ƶ������ṹ
figure(1);
stem(0:IFFT_bin_length-1, abs(IFFT_modulation(2,1:IFFT_bin_length)),'b*-')
grid on
axis ([0 IFFT_bin_length -0.5 4.5]);
ylabel('Magnitude');
xlabel('IFFT Bin');
title('OFDM���ز�Ƶ�ʷ���');
 
% Figure 2: IFFT��Ƶ�����λ���ȣ�
% - ���ֹ���ӳ���������λ�Գ��ԣ���֤ʵֵʱ���ź�����
figure(2);
plot(0:IFFT_bin_length-1, (180/pi)*angle(IFFT_modulation(2,1:IFFT_bin_length)), 'go')
hold on
stem(carriers-1, (180/pi)*angle(IFFT_modulation(2,carriers)),'b*-');
stem(conjugate_carriers-1, (180/pi)*angle(IFFT_modulation(2,conjugate_carriers)),'b*-');
axis ([0 IFFT_bin_length -200 +200])
grid on
ylabel('Phase (degrees)')
xlabel('IFFT Bin')
title('OFDM���ز���λ')

%------------------------------------------------------------------------------
% ����5: IFFT�任���õ�ʱ��OFDM���ţ�δ�ӱ��������
%------------------------------------------------------------------------------
signal_after_IFFT=ifft(IFFT_modulation,IFFT_bin_length,2);
time_wave_matrix=signal_after_IFFT;

% Figure 3: ����OFDM���ţ�δ��ѭ��ǰ׺/��׺����ʱ����
% - չʾIFFT�����һ���������ڰ�������ȷ�Χ
figure(3)
plot(0:IFFT_bin_length-1,time_wave_matrix(2,:));
axis([0 IFFT_bin_length -0.2 0.2]);
grid on
ylabel('Amplitude');
xlabel('Time');
title('OFDMʱ���źţ�����������');

%------------------------------------------------------------------------------
% ����6: ����ѭ��ǰ׺(CP)���׺
%------------------------------------------------------------------------------
% CP���ã��Կ��ྶ���ţ����ַ��ż�������
% ��׺���ã����ڼӴ���ƽ�����ɣ�����Ƶ���԰�
XX=zeros(symbols_per_carrier,IFFT_bin_length+GI+GIP);
for k=1:symbols_per_carrier
    % �������岿�֣��м䣩
    for i=1:IFFT_bin_length
        XX(k,i+GI)=signal_after_IFFT(k,i);
    end
    % ѭ��ǰ׺��������β�����Ƶ���ͷ
    for i=1:GI
        XX(k,i)=signal_after_IFFT(k,i+IFFT_bin_length-GI);
    end
    % ��׺��������ͷ�����Ƶ�ĩβ�����ڴ����Ҳ���ɣ�
    for j=1:GIP
         XX(k,IFFT_bin_length+GI+j)=signal_after_IFFT(k,j);
    end
end  
time_wave_matrix_cp=XX;

% Figure 4: ����OFDM��������ѭ��ǰ׺(CP)���׺���ʱ����
% - ���CP���Ҷ˺�׺�����ڴ��ڹ��ɣ����۲��������
figure(4);
plot(0:size(time_wave_matrix_cp,2)-1, time_wave_matrix_cp(2,:));
axis([0, size(time_wave_matrix_cp,2), -0.2, 0.2]);
grid on;
ylabel('Amplitude');
xlabel('Time');
title('OFDMʱ���źţ���ѭ��ǰ׺��������������');

%------------------------------------------------------------------------------
% ����7: OFDM���żӴ�������LTE���
%------------------------------------------------------------------------------
% LTEԭ����ʹ�������Ҵ�ƽ�����ű�Ե�����Ʒ��ű�Ե�Ķ������䣬���ٴ������
% ���������Ƿ�Χ��[CP(GI) | ����(N)]���ܳ���Ϊ N+GI
% �����������Ե��CP��ʼ�����ұ�Ե������������ṩƽ������
% ��׺���֣�GIP����������һ�����ŵ�CP�ص���ʵ��ƽ��ƴ��
windowed_time_wave_matrix_cp=zeros(symbols_per_carrier,IFFT_bin_length+GI+GIP);

% ���������Ҵ�����������CP+���岿�֣�����ΪN+GI��
% ע�⣺rcoswindow���صĳ�����(1+beta)*(N+GI)��������ҪֻȡǰN+GI��Ԫ��
rcos_win_full = rcoswindow(beta, IFFT_bin_length+GI);  % ������������ = (1+beta)*(N+GI)
rcos_win = rcos_win_full(1:IFFT_bin_length+GI)';  % ֻȡǰN+GI��Ԫ�أ�ת��Ϊ������

% ��ÿ�����żӴ�
for i = 1:symbols_per_carrier
    % ���Žṹ��[CP(GI) | ����(N) | ��׺(GIP)]
    % ������Ӧ����ǰ N+GI ��������CP+���壩
    windowed_time_wave_matrix_cp(i, 1:IFFT_bin_length+GI) = ...
        real(time_wave_matrix_cp(i, 1:IFFT_bin_length+GI)) .* rcos_win;
    
    % ��׺���֣�GIP��������������ԭֵ����������һ�����ŵ�CP�ص�
    % ע�⣺��׺���ֲ��Ӵ�����Ϊ��������һ�����ŵ�CP�ص����
    if GIP > 0
        windowed_time_wave_matrix_cp(i, IFFT_bin_length+GI+1:IFFT_bin_length+GI+GIP) = ...
            real(time_wave_matrix_cp(i, IFFT_bin_length+GI+1:IFFT_bin_length+GI+GIP));
    end
end

% Figure 5: �Ӵ���ĵ���OFDM����ʱ���Σ���CP/��׺��
% - LTE��������Ҵ�ƽ�����أ�����Ƶ���԰�
% - ������Ӧ����CP+���岿�֣���׺������������һ�����ŵ�CP�ص�
figure(5)
plot(0:IFFT_bin_length-1+GI+GIP,windowed_time_wave_matrix_cp(2,:)); 
axis([0, IFFT_bin_length-1+GI+GIP, -0.2, 0.2]);
grid on;
ylabel('Amplitude');
xlabel('Time');
title(sprintf('OFDMʱ���źţ�LTE���Ӵ����ص�����������CP�ڣ�GIP=%d��������������', GIP));

%------------------------------------------------------------------------------
% ����8: ���ɷ����źţ������任����֡��֯��LTE����ص���ӣ�
%------------------------------------------------------------------------------
% LTEԭ��������ƴ��ʱ��CP��Χ���ص���ӣ�ʵ��ƽ������
% ֡�ṹ��ÿ֡�������OFDM���ţ�����֡ĩβ����һ�κ�׺
% 
% ���Žṹ��[CP(GI) | ����(N) | ��׺(GIP)]
% �ص����ƣ�
%   - ��ǰ���ŵĺ�׺(GIP������)����һ�����ŵ�CPǰGIP�������ص�
%   - �ص������ڴ��������У��������ŵķ������
%   - �ص�����������CP�ڣ�GIP <= GI������Ӱ����Ч���ݣ����岿�֣�
% 
% �������нṹ��
%   [����1: CP+����+��׺] [����2: CP+����+��׺] ... [����N: CP+����+��׺]
%   ���У�����i�ĺ�׺�����i+1��CPǰGIP�������ص����
num_frames = total_symbols / symbols_per_frame;
frame_len_CP_suffix = symbols_per_frame*(IFFT_bin_length+GI)+GIP; % ÿ֡���г��ȣ���ĩβһ�κ�׺��

% ��֡���죺LTE����ص����
Tx_data = zeros(1, num_frames*frame_len_CP_suffix);

write_offset = 0;
for f = 1:num_frames
    sym_start = (f-1)*symbols_per_frame + 1;
    sym_end   = f*symbols_per_frame;

    % ��ǰ֡�ļӴ����ž���
    frame_windowed = windowed_time_wave_matrix_cp(sym_start:sym_end, :);

    % LTE����ص���ӣ�����֡�������ź�
    frame_serial_windowed = zeros(1, frame_len_CP_suffix);
    
    % ��һ�����ţ�����д�루����CP+����+��׺��
    frame_serial_windowed(1:IFFT_bin_length+GI+GIP) = frame_windowed(1, :);
    
    % �������ţ��ص���Ӵ���
    for i = 1:(symbols_per_frame-1)
        % ��һ�������ڴ��������е���ʼλ��
        % ��i������ռ�ݣ�[1+(i-1)*(N+GI), i*(N+GI)+GIP]
        % ��i+1������Ӧ�ôӣ�i*(N+GI)+1 ��ʼ�����i�����ŵĺ�׺�ص�GIP��������
        next_symbol_start = i*(IFFT_bin_length+GI) + 1;
        next_symbol_end = (i+1)*(IFFT_bin_length+GI) + GIP;
        
        % LTE�ص����ƣ�
        % �ص����򣺵�ǰ���ŵĺ�׺(GIP)����һ�����ŵ�CPǰGIP�������ص�
        % �ص�λ�ã����������е� [next_symbol_start, next_symbol_start+GIP-1]
        % ��ǰ���ŵĺ�׺��frame_windowed(i, IFFT_bin_length+GI+1:IFFT_bin_length+GI+GIP)
        % ��һ�����ŵ�CPǰGIP��������frame_windowed(i+1, 1:GIP)
        
        if GIP > 0 && next_symbol_end <= length(frame_serial_windowed)
            % �ص����򣺵�ǰ���ŵĺ�׺����д�룩+ ��һ�����ŵ�CPǰGIP����������д�룩
            overlap_start = next_symbol_start;
            overlap_end = overlap_start + GIP - 1;
            
            % ��ǰ���ŵĺ�׺���֣���frame_serial_windowed����д�룩
            current_suffix = frame_windowed(i, IFFT_bin_length+GI+1:IFFT_bin_length+GI+GIP);
            
            % ��һ�����ŵ�CPǰGIP������
            next_cp_prefix = frame_windowed(i+1, 1:GIP);
            
            % LTE�ص���ӣ����ص������������ŵķ������
            frame_serial_windowed(overlap_start:overlap_end) = ...
                current_suffix + next_cp_prefix;
            
            % ���ص����֣�д����һ�����ŵ�ʣ�ಿ�֣�CP��ʣ�ಿ��+����+��׺��
            if overlap_end < next_symbol_end
                non_overlap_start = overlap_end + 1;
                non_overlap_end = next_symbol_end;
                % ��һ�����Ŵ�GIP+1��ʼ����β
                frame_serial_windowed(non_overlap_start:non_overlap_end) = ...
                    frame_windowed(i+1, GIP+1:IFFT_bin_length+GI+GIP);
            end
        else
            % ���û�к�׺��GIP=0����ֱ��д����һ������
            if next_symbol_end <= length(frame_serial_windowed)
                frame_serial_windowed(next_symbol_start:next_symbol_end) = frame_windowed(i+1, :);
            end
        end
    end

    % д�뵽��֡��������
    Tx_data(write_offset+1:write_offset+frame_len_CP_suffix) = frame_serial_windowed;
    write_offset = write_offset + frame_len_CP_suffix;
end

% �����ֱ�Ӵ��ӣ���ÿ���ź�׺�������ڶԱ���"δ�Ӵ�"Ƶ�׷���
Tx_data_withoutwindow = reshape(time_wave_matrix_cp', (total_symbols)*(IFFT_bin_length+GI+GIP), 1)';

% ���ִ��ӷ�ʽ�Աȣ�����Figure 6��
temp_time_symbol = length(Tx_data_withoutwindow); % ÿ���Ŵ���׺
temp_time_frame  = length(Tx_data);               % ÿ֡һ�κ�׺

% Figure 6: ���Ͷ˴���ʱ���źţ����ִ��ӷ�ʽ�Աȣ�
% - ��ͼ(1): ��ÿ������(��CP/��׺)ֱ�Ӵ���
% - ��ͼ(2): ����ÿ֡ĩβ����һ�κ�׺�ļӴ�����
figure (6)
subplot(2,1,1);
plot(0:temp_time_symbol-1,Tx_data_withoutwindow);
grid on
ylabel('Amplitude (volts)')
xlabel('Time (samples)')
title('OFDMʱ���ź�')
temp_time2 = temp_time_frame;
subplot(2,1,2);
plot(0:temp_time2-1,Tx_data);
grid on
ylabel('Amplitude (volts)')
xlabel('Time (samples)')
title('OFDMʱ���ź�')
%===============================================================================
% ���źŷ�������ӻ����ֲ�����Ƶ����ʱ��Աȡ�
%===============================================================================

%------------------------------------------------------------------------------
% ����9: �ֲ����ڶ��壨������ϸ�������� Figure 12 ��ͬ��Χ��
%------------------------------------------------------------------------------
zoom_len = 5000; % �ֲ����ڳ��ȣ��� Figure 12 ��ͬ��
L = length(Tx_data);
zs = max(1, floor(L/2) - floor(zoom_len/2));
ze = min(L, zs + zoom_len - 1);
% Ϊδ�Ӵ��źż����Ӧ�Ĵ��ڷ�Χ
L2 = length(Tx_data_withoutwindow);
zs2 = max(1, floor(L2/2) - floor(zoom_len/2));
ze2 = min(L2, zs2 + zoom_len - 1);

% ��ȡ�ֲ�����ʱ���ź�
subset_ofdm_no_window = Tx_data_withoutwindow(zs2:ze2);  % δ�Ӵ��źžֲ�����
subset_ofdm_window = Tx_data(zs:ze);                      % �Ӵ��źžֲ�����

% ����Ƶ������ף�˫���ף�������Ƶ�ʣ�
Nfft_no_window = length(subset_ofdm_no_window);
Nfft_window = length(subset_ofdm_window);
subset_ofdm_f_no_window = fftshift(fft(subset_ofdm_no_window));
subset_ofdm_f_log_no_window = 20*log10(abs(subset_ofdm_f_no_window)/max(abs(subset_ofdm_f_no_window)) + eps);
subset_ofdm_f_window = fftshift(fft(subset_ofdm_window));
subset_ofdm_f_log_window = 20*log10(abs(subset_ofdm_f_window)/max(abs(subset_ofdm_f_window)) + eps);

%------------------------------------------------------------------------------
% Figure 7: �Ӵ�ǰ��ֲ�����ʱ��Աȣ��� Figure 12 ��ͬ���ڷ�Χ��
%------------------------------------------------------------------------------
figure (7)
subplot(2,1,1)
plot(real(subset_ofdm_no_window), 'b-', 'LineWidth', 1.5)
grid on
xlabel('������')
ylabel('����')
title(sprintf('δ�Ӵ��ź�ʱ�� (window [%d:%d])', zs2, ze2))

subplot(2,1,2)
plot(real(subset_ofdm_window), 'r-', 'LineWidth', 1.5)
grid on
xlabel('������')
ylabel('����')
title(sprintf('�Ӵ��ź�ʱ�� (window [%d:%d])', zs, ze))

%------------------------------------------------------------------------------
% Figure 8: �Ӵ�ǰ��ֲ�����Ƶ��Աȣ��� Figure 12 ��ͬ���ڷ�Χ��˫���ף�
%------------------------------------------------------------------------------
% ��һ��Ƶ���᣺-0.5..0.5����Ӧ -fs/2..fs/2��
f_norm_no_window = (-Nfft_no_window/2:(Nfft_no_window/2-1))/Nfft_no_window;
f_norm_window = (-Nfft_window/2:(Nfft_window/2-1))/Nfft_window;

figure(8)
subplot(2,1,1)
plot(f_norm_no_window, subset_ofdm_f_log_no_window, 'b-', 'LineWidth', 1.5)
grid on
axis([-0.5 0.5 -40 max(subset_ofdm_f_log_no_window)])
ylabel('���� (dB)')
xlabel('��һ��Ƶ�� (-0.5..0.5 = ��fs/2)')
title(sprintf('δ�Ӵ��źŷ����ף�˫���ף�window [%d:%d]��', zs2, ze2))

subplot(2,1,2)
plot(f_norm_window, subset_ofdm_f_log_window, 'r-', 'LineWidth', 1.5)
grid on
axis([-0.5 0.5 -40 max(subset_ofdm_f_log_window)])
ylabel('���� (dB)')
xlabel('��һ��Ƶ�� (-0.5..0.5 = ��fs/2)')
title(sprintf('�Ӵ��źŷ����ף�˫���ף�window [%d:%d]��', zs, ze))

%===============================================================================
% ���ŵ����䣺AWGN���Ը�˹�������ŵ���
%===============================================================================
% �������ʼ��㣺����Ŀ��SNR�ͷ����źŹ��ʼ�����������
Tx_signal_power=var(Tx_data); % �ԼӴ���ķ��ʹ���Ϊ��׼���㹦��
linear_SNR=10^(targetSNRdB/10);
noise_sigma=Tx_signal_power/linear_SNR;
noise_scale_factor=sqrt(noise_sigma);
noise=randn(1,length(Tx_data))*noise_scale_factor;
Rx_data=Tx_data+noise;

%------------------------------------------------------------------------------
% ����10: �ֲ�����ʱ��Աȣ�����ǰ��
%------------------------------------------------------------------------------
% ע�����ڷ�Χ (zs, ze) ���� Figure 7/8 �����壬�˴�ֱ��ʹ��
figure(12);
subplot(2,1,1);
plot(zs:ze, Tx_data(zs:ze));
grid on;
ylabel('Amplitude (volts)');
xlabel('Sample index');
title(sprintf('�����źţ��ֲ��Ŵ�[%d:%d]', zs, ze));

subplot(2,1,2);
plot(zs:ze, Rx_data(zs:ze));
grid on;
ylabel('Amplitude (volts)');
xlabel('Sample index');
title(sprintf('�����źţ��ֲ��Ŵ�[%d:%d]', zs, ze));

%------------------------------------------------------------------------------
% ����11: ���㲢��ӡ�ֲ����ڵ� SNR �� MSE
%------------------------------------------------------------------------------
tx_seg = Tx_data(zs:ze);
rx_seg = Rx_data(zs:ze);
noise_seg = rx_seg - tx_seg;

sig_power  = mean(tx_seg.^2);
noise_power = mean(noise_seg.^2);
snr_zoom_db = 10*log10(sig_power / max(noise_power, eps));
mse_zoom = noise_power; % �������=������������

fprintf('\n==== Zoomed Window Metrics (Figure 12) ====\n');
fprintf('Window range       : [%d : %d] (length=%d)\n', zs, ze, numel(tx_seg));
fprintf('Signal power       : %.6g\n', sig_power);
fprintf('Noise power        : %.6g\n', noise_power);
fprintf('SNR (Tx vs Rx)     : %.4f dB\n', snr_zoom_db);
fprintf('MSE (Tx vs Rx)     : %.6g\n', mse_zoom);
fprintf('===========================================\n');

%------------------------------------------------------------------------------
% ����12: Figure 13 - �� Figure 12 ͬһ���ڵ�Ƶ�׶Աȣ�˫���ף�
%------------------------------------------------------------------------------
Nfft_zoom = 2^nextpow2(length(tx_seg));
Tx_Fz = fftshift(fft(tx_seg, Nfft_zoom));
Rx_Fz = fftshift(fft(rx_seg, Nfft_zoom));

% ��һ�����������ֵ��ת dB
Tx_mag_db_z = 20*log10( abs(Tx_Fz)/max(abs(Tx_Fz)) + eps );
Rx_mag_db_z = 20*log10( abs(Rx_Fz)/max(abs(Rx_Fz)) + eps );

% ��һ��Ƶ���᣺-0.5..0.5����Ӧ -fs/2..fs/2��
f_norm_z = (-Nfft_zoom/2:(Nfft_zoom/2-1))/Nfft_zoom;

figure(13);
subplot(2,1,1);
plot(f_norm_z, Tx_mag_db_z, 'b');
grid on;
ylabel('Magnitude (dB)');
xlabel('Normalized Frequency (-0.5..0.5 = \pm fs/2)');
title(sprintf('�����ź�Ƶ�ף�˫���ף��ֲ��Ŵ� [%d:%d]��', zs, ze));

subplot(2,1,2);
plot(f_norm_z, Rx_mag_db_z, 'r');
grid on;
ylabel('Magnitude (dB)');
xlabel('Normalized Frequency (-0.5..0.5 = \pm fs/2)');
title(sprintf('�����ź�Ƶ�ף�˫���ף��ֲ��Ŵ� [%d:%d]��', zs, ze));

%===============================================================================
% �����ն˴������̡�
%===============================================================================

%------------------------------------------------------------------------------
% ����13: ����ת����ȥ��ѭ��ǰ׺�ͺ�׺��LTE���
%------------------------------------------------------------------------------
% LTEԭ�������ն�ȥ��CPʱֻȡ���岿�֣������ص�Ӱ��
% - ���Ͷˣ����ŽṹΪ [CP(GI) | ����(N) | ��׺(GIP)]���ص�������CP��
% - ���նˣ�ȥ��CP�ͺ�׺��ֻ��ȡ���岿�֣�N������������FFT���
% - �ص�����Ӱ����Ч���ݣ���Ϊ�ص�ֻ��CP��Χ�ڣ����岿�ֲ���Ӱ��
symbol_len = IFFT_bin_length+GI+GIP;  % ÿ�����ŵ��ܳ��ȣ���CP+����+��׺��

Rx_data_matrix=zeros(total_symbols,symbol_len);
read_offset = 0;
for f = 1:num_frames
    frame_rx = Rx_data(read_offset+1:read_offset+frame_len_CP_suffix);
    for i = 1:symbols_per_frame
        global_sym_idx = (f-1)*symbols_per_frame + i;
        % ��ȡ���ţ��Ӵ�����������ȡÿ�����ţ���CP+����+��׺��
        start_idx = (i-1)*(IFFT_bin_length+GI) + 1;
        end_idx = i*(IFFT_bin_length+GI) + GIP;
        Rx_data_matrix(global_sym_idx,:) = frame_rx(start_idx:end_idx);
    end
    read_offset = read_offset + frame_len_CP_suffix;
end
% LTEԭ����ȥ��CP�ͺ�׺��ֻ�������岿�֣�N������������FFT���
% �ص�������CP�ڣ�ȥ��CP��Ӱ����Ч����
Rx_data_complex_matrix=Rx_data_matrix(:,GI+1:IFFT_bin_length+GI); 

%------------------------------------------------------------------------------
% ����14: FFT�����ʱ�� �� Ƶ�򣩣���ȡ���ز�����
%------------------------------------------------------------------------------
Y1=fft(Rx_data_complex_matrix,IFFT_bin_length,2);
Rx_carriers=Y1(:,carriers);
Rx_mag=abs(Rx_carriers);
Rx_phase=angle(Rx_carriers);

% ������תֱ������
[M, N]=pol2cart(Rx_phase, Rx_mag); 
Rx_complex_carrier_matrix = complex(M, N);
% Figure 9: ���ն����ز�����ͼ��IQƽ�棩��SNRԽ��ɢ��Խ������16QAM�����
figure(9);
plot(Rx_complex_carrier_matrix,'*b');
axis([-4, 4, -4, 4]);
title('�����ź�16QAM����ͼ')
grid on

%------------------------------------------------------------------------------
% ����15: 16QAM�������Сŷ�Ͼ����о���
%------------------------------------------------------------------------------
Rx_serial_complex_symbols = reshape(Rx_complex_carrier_matrix',size(Rx_complex_carrier_matrix, 1)*size(Rx_complex_carrier_matrix,2),1)' ; 
Rx_decoded_binary_symbols=demoduqam16(Rx_serial_complex_symbols);
baseband_in = Rx_decoded_binary_symbols;
% Figure 10: �������Աȣ�ǰ100���أ�������:���ͣ���:�����о�
figure(10);
subplot(2,1,1);
stem(baseband_out(1:100));
title('���ͱ�������ǰ100λ��')
subplot(2,1,2);
stem(baseband_in(1:100));
title('���ձ�������ǰ100λ��')
%------------------------------------------------------------------------------
% ����16: �����ʼ��㣨15dB��������ʾժҪ��
%------------------------------------------------------------------------------
bit_errors=find(baseband_in ~=baseband_out);
bit_error_count = size(bit_errors, 2); 
ber_15dB=bit_error_count/baseband_out_length;

%------------------------------------------------------------------------------
% ����17: ����������ؼ�ָ��
%------------------------------------------------------------------------------
null_subcarriers = IFFT_bin_length - 2*carrier_count; % ����DC/������

fprintf('\n==== Simulation Summary (SNR = %.2f dB) ====\n', targetSNRdB);
fprintf('Total symbols     : %d\n', total_symbols);
fprintf('Total bits         : %d\n', baseband_out_length);
fprintf('BER               : %.6g\n', ber_15dB);
fprintf('Tx power (var)    : %.6g\n', Tx_signal_power);
fprintf('Noise variance    : %.6g\n', noise_sigma);
fprintf('IFFT length (N)   : %d\n', IFFT_bin_length);
fprintf('Active carriers   : %d\n', carrier_count);
fprintf('Null carriers     : %d\n', null_subcarriers);
fprintf('CP length (GI)    : %d (ratio=%.3f)\n', GI, GI/IFFT_bin_length);
fprintf('Suffix length     : %d (GIP <= GI��������CP��)\n', GIP);
fprintf('One OFDM symbol   : %d samples (N+GI+GIP��LTE���)\n', IFFT_bin_length+GI+GIP);
fprintf('Window roll-off   : 1/%d\n', round(1/beta));
fprintf('Solution method   : %d (LTE����ص���ӣ��ص�����������CP��)\n', solution_method);
fprintf('Overlap length    : %d samples (GIP <= GI��������CP��)\n', GIP);
fprintf('Symbols per frame : %d\n', symbols_per_carrier);
fprintf('==============================\n\n');

%===============================================================================
% ��BER-SNR�������߼��㡿
%===============================================================================
% ˵������10-30dB��Χ�ڣ�����2dB�������SNR���������
fprintf('\n==== Computing BER-SNR Curve (10-30 dB, step 2 dB) ====\n');
SNR_range = 10:2:30;  % ����ȷ�Χ
ber_results = zeros(size(SNR_range));  % �洢BER���

for idx = 1:length(SNR_range)
    snr_dB = SNR_range(idx);
    
    %--------------------------------------------------------------------------
    % ��ÿ��SNR�㣺���¼������������ӵ������ź�
    %--------------------------------------------------------------------------
    linear_SNR = 10^(snr_dB/10);
    noise_sigma_loop = Tx_signal_power/linear_SNR;
    noise_scale_factor_loop = sqrt(noise_sigma_loop);
    noise_loop = randn(1, length(Tx_data)) * noise_scale_factor_loop;
    Rx_data_loop = Tx_data + noise_loop;
    
    % ���ն˴���ת������֡��ȥ��ѭ��ǰ׺�ͺ�׺��LTE���
    symbol_len_loop = IFFT_bin_length+GI+GIP;  % ÿ�����ŵ��ܳ��ȣ���CP+����+��׺��
    
    Rx_data_matrix_loop = zeros(total_symbols, symbol_len_loop);
    read_offset = 0;
    for f = 1:num_frames
        frame_rx = Rx_data_loop(read_offset+1:read_offset+frame_len_CP_suffix);
        for i = 1:symbols_per_frame
            global_sym_idx = (f-1)*symbols_per_frame + i;
            % ��ȡ���ţ��Ӵ�����������ȡÿ�����ţ���CP+����+��׺��
            start_idx = (i-1)*(IFFT_bin_length+GI) + 1;
            end_idx = i*(IFFT_bin_length+GI) + GIP;
            Rx_data_matrix_loop(global_sym_idx,:) = frame_rx(start_idx:end_idx);
        end
        read_offset = read_offset + frame_len_CP_suffix;
    end
    % LTEԭ����ȥ��CP�ͺ�׺��ֻ�������岿�֣�N������������FFT���
    Rx_data_complex_matrix_loop = Rx_data_matrix_loop(:,GI+1:IFFT_bin_length+GI);
    
    % FFT�������ȡ���ز�����
    Y1_loop = fft(Rx_data_complex_matrix_loop, IFFT_bin_length, 2);
    Rx_carriers_loop = Y1_loop(:,carriers);
    Rx_mag_loop = abs(Rx_carriers_loop);
    Rx_phase_loop = angle(Rx_carriers_loop);
    [M_loop, N_loop] = pol2cart(Rx_phase_loop, Rx_mag_loop);
    Rx_complex_carrier_matrix_loop = complex(M_loop, N_loop);
    
    % 16QAM���
    Rx_serial_complex_symbols_loop = reshape(Rx_complex_carrier_matrix_loop', size(Rx_complex_carrier_matrix_loop, 1)*size(Rx_complex_carrier_matrix_loop, 2), 1)';
    Rx_decoded_binary_symbols_loop = demoduqam16(Rx_serial_complex_symbols_loop);
    baseband_in_loop = Rx_decoded_binary_symbols_loop;
    
    % ����BER
    bit_errors_loop = find(baseband_in_loop ~= baseband_out);
    bit_error_count_loop = size(bit_errors_loop, 2);
    ber_results(idx) = bit_error_count_loop / baseband_out_length;
    
    fprintf('SNR = %2d dB: BER = %.6g\n', snr_dB, ber_results(idx));
end
fprintf('==========================================================\n');

%------------------------------------------------------------------------------
% ����18: Figure 11 - BER-SNR�������߻���
%------------------------------------------------------------------------------
% ���� BER Ϊ 0 ��������ú�С��ֵ�����Ա��ڶ�����������ʾ
ber_plot = ber_results;
ber_plot(ber_plot == 0) = 1e-10;  % �� 0 �滻Ϊ��С��ֵ

figure(11);
semilogy(SNR_range, ber_plot, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('SNR (dB)');
ylabel('BER');
title('OFDMϵͳ����������');
grid on

% ���ú������᷶Χ
non_zero_ber = ber_results(ber_results > 0);
if ~isempty(non_zero_ber)
    y_min = min(non_zero_ber) * 0.5;
    y_max = max([ber_results, 1]) * 1.1;
else
    y_min = 1e-6;
    y_max = 1;
end
axis([9 31 y_min y_max])

toc
%===============================================================================
% �ļ�����
%===============================================================================











