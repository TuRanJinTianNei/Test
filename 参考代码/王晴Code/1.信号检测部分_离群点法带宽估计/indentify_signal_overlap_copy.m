%信号重叠性检测 - 参数对齐 double_signal_loc2_test1.m
%-----------------------------------修改说明-------------------------------------------------------
% 1. 参数对齐：与 double_signal_loc2_test1.m 保持一致
%    - SymbN = 60 (原20)
%    - SymbN2 = 60 (原20)
%    - freq_shift_ratio = 0.125 (原0.15)
%    - Amp_Scale_Sig2 = 0.6 (原PowerScale=4，现在是幅度缩小)
% 2. 目的：使用相同的信号参数进行重叠性检测，便于结果对比。
%--------------------------------------------------------------------------------------------------
clc
clear all
close all
SNR=15;                  % 提高一点SNR以保证弱信号在噪声中可见，主要测试大信号干扰下的逻辑
kkk=1;                   
indx=zeros(1,kkk);

% 控制是否显示过程图形
SHOW_FIGURES = 1;        

fprintf('========================================\n');
fprintf('信号重叠性检测 - 参数对齐 double_signal_loc2_test1.m\n');
fprintf('SNR: %.0f dB\n', SNR);
fprintf('SymbN: 60, freq_shift: 0.125, Amp_Scale: 0.6\n');
fprintf('========================================\n\n');

%% 信号2 (强干扰信号)
% --- 参数对齐 double_signal_loc2_test1.m ---
SubCarryNN2=256;      
SubCarryN2=1024;      
ratio=1/4;           
fftLen=1024;         
SymbN2=60;            % 从20改为60，与double_signal_loc2_test1.m一致
GuardLen=SubCarryN2*ratio;    

% 生成信号2
SignalLen2=SubCarryNN2*SymbN2*2;           
Signal2=round(rand(1,SignalLen2));       
ParaBitSig2=reshape(Signal2,SubCarryNN2,SymbN2*2);
for jj=1:SymbN2
    ich22(:,jj)=ParaBitSig2(:,2*jj-1); 
    qch22(:,jj)=ParaBitSig2(:,2*jj);   
end
kmod=1./sqrt(2);   
ich02=ich22.*2-1;     
qch02=qch22.*2-1;
ich12=ich02.*kmod;   
qch12=qch02.*kmod;
x2=ich12+qch12.*sqrt(-1);      
% 映射 (适配256)
x2=[x2(1:SubCarryNN2/2,1:SymbN2);zeros(SubCarryN2-SubCarryNN2,SymbN2);x2(SubCarryNN2/2+1:SubCarryNN2,1:SymbN2)];
y2=ifft(x2);   
ich22=real(y2);   
qch22=imag(y2);   
ich32=[ich22(fftLen-GuardLen+1:fftLen,:);ich22];
qch32=[qch22(fftLen-GuardLen+1:fftLen,:);qch22];
ich42=reshape(ich32,1,(fftLen+GuardLen)*SymbN2);
qch42=reshape(qch32,1,(fftLen+GuardLen)*SymbN2);
TrData42=ich42+qch42.*sqrt(-1);  

for index=1:length(SNR)
    fprintf('  正在仿真... \n');
    tic;
    for p=1:kkk
        %% 信号1 (弱目标信号)
        SubCarryNN=256;      
        SubCarryN=1024;      
        ratio=1/4;           
        fftLen=1024;         
        SymbN=60;            % 从20改为60，与double_signal_loc2_test1.m一致
        GuardLen=SubCarryN*ratio;    
        
        SignalLen=SubCarryNN*SymbN*2;           
        Signal=round(rand(1,SignalLen));       
        ParaBitSig=reshape(Signal,SubCarryNN,SymbN*2);
        for j=1:SymbN
            ich(:,j)=ParaBitSig(:,2*j-1); 
            qch(:,j)=ParaBitSig(:,2*j);   
        end
        kmod=1./sqrt(2);   
        ich0=ich.*2-1;     
        qch0=qch.*2-1;
        ich1=ich0.*kmod;   
        qch1=qch0.*kmod;
        x=ich1+qch1.*sqrt(-1);      
        x=[x(1:SubCarryNN/2,1:SymbN);zeros(SubCarryN-SubCarryNN,SymbN);x(SubCarryNN/2+1:SubCarryNN,1:SymbN)];
        y=ifft(x);   
        ich2=real(y);   
        qch2=imag(y);   
        ich3=[ich2(fftLen-GuardLen+1:fftLen,:);ich2];
        qch3=[qch2(fftLen-GuardLen+1:fftLen,:);qch2];
        ich4=reshape(ich3,1,(fftLen+GuardLen)*SymbN);
        qch4=reshape(qch3,1,(fftLen+GuardLen)*SymbN);
        TrData=ich4+qch4.*sqrt(-1);  
        
        % 加噪
        ReData=awgn(TrData,SNR(index),'measured');
        ReData42_noisy=awgn(TrData42,SNR(index),'measured'); 
        
        % --- 参数对齐 double_signal_loc2_test1.m ---
        freq_shift_ratio = 0.125; % 从0.15改为0.125，与double_signal_loc2_test1.m一致
        Amp_Scale_Sig2 = 0.6;     % 从PowerScale=4改为0.6，与double_signal_loc2_test1.m一致（幅度缩放）
        
        ReData_same = ReData42_noisy .* exp(1j * freq_shift_ratio * pi * 2 .* (1:length(ReData42_noisy)));
        ReData_same = ReData_same * Amp_Scale_Sig2;  % 幅度缩放0.6倍 
        
        ReData = ReData + ReData_same;
        
        %% 频谱计算
        fs=1024;
        f=linspace(-fs/2,fs/2,25600);
        f=f/(2*max(f));
        Fs=1;
        nfft=1024;
        window=hanning(150);
        noverlap=50;
        [Pxx1,f]=pwelch(ReData,window,noverlap,length(ReData),Fs);
        Pxx=10*log10(fftshift(Pxx1));
        
        % 图形1：功率谱
        if SHOW_FIGURES && p == 1
            figure(1);
            set(gcf,'Position',[50,50,1400,900],'Color','white');
            sgtitle('信号重叠性检测 (参数对齐 double_signal_loc2_test1.m)','FontSize',16,'FontWeight','bold');
            subplot(2,3,1);
            plot(f,Pxx,'-','LineWidth',2,'Color',[0.2 0.4 0.8]);
            xlabel('归一化频率'); ylabel('PSD (dB)');
            title('步骤1: 功率谱 (可见明显的台阶)');
            grid on;
        end
        
        %% 功率分布与拐点
        Pxx=Pxx+abs(min(Pxx));
        Pxx_int=zeros(size(Pxx,1),1);
        for x=1:1:size(Pxx,1)
            if x==1
                Pxx_int(x)=(1/(size(Pxx,1)-1))*Pxx(x+1);
            else
                Pxx_int(x)=Pxx_int(x-1)+(1/(size(Pxx,1)-1))*Pxx(x);
            end
        end
        
        y_=Pxx_int;
        x_=(linspace(0,1,length(Pxx_int)))';
        loc=[x_,y_];
        
        % 差分计算
        d = zeros(1, length(y_));
        for i=4:length(y_)-4
            n1(i)=3*(loc(i,2))-(loc(i-1,2))-(loc(i-2,2))-(loc(i-3,2));
            n2(i)=loc(i+3,2)+loc(i+2,2)+loc(i+1,2)-3*(loc(i,2));
            d(i)=n2(i)-n1(i);
        end
        
        % 图形3：差分结果
        if SHOW_FIGURES && p == 1
            figure(1);
            subplot(2,3,3);
            plot(d,'-','LineWidth',1.5,'Color',[0 0.7 0]);
            title('步骤3: 差分结果 (强弱峰值差异明显)');
            grid on; axis tight;
        end
        
        %% VPD与峰值提取
        x=d;
        row_acc = abs(x);
        m = length(row_acc);
        
        % 均值滤波
        row_acc1 = row_acc;
        for i=2:m-1
            row_acc1(i)=(row_acc(i-1) + row_acc(i)+row_acc(i+1))/3;
        end
        for i=m-1:-1:2
            row_acc(i) = (row_acc1(i-1) + row_acc1(i)+row_acc1(i+1))/3;
        end
        
        % 极值查找
        peaks = zeros(1,m); peakindexs = zeros(1,m);
        valleys = zeros(1,m); valleyindexs = zeros(1,m);
        peakindex = 1; valleyindex = 1;
        
        for i = 2:m-1
            if row_acc(i) >row_acc(i-1) && row_acc(i)>=row_acc(i+1)
                peaks(peakindex)=row_acc(i);
                peakindexs(peakindex)=i;
                peakindex = peakindex+1;
            end
            if row_acc(i) < row_acc(i-1) && row_acc(i)<row_acc(i+1)
                valleys(valleyindex)=row_acc(i);
                valleyindexs(valleyindex)=i;
                valleyindex=valleyindex+1;
            end
        end
        
        % 简单VPD筛选 (此处简化逻辑，保留原算法核心)
        pcount = peakindex-1;
        indices = peakindexs(1:pcount);
        indices = indices(indices>0);
        if ~isempty(indices)
            peaks=row_acc(indices);
        else
            peaks=[];
        end

        % 图形4：VPD峰值
        if SHOW_FIGURES && p == 1
            figure(1);
            subplot(2,3,4);
            plot(row_acc,'-','Color',[0.6 0.6 0.6]); hold on;
            plot(indices, peaks, 'r.', 'MarkerSize', 10);
            title('步骤4: 峰值检测');
            grid on; axis tight;
        end

        %% 峰值排序与逻辑判断
        if ~isempty(peaks)
            peaks_dou=abs(peaks);
            [~, ia] = unique(peaks_dou); % 去重
            peaks_2=peaks(ia); 
            peaks_2=sort(peaks_2); % 从小到大排序
        else
            peaks_2 = [];
        end
        
        % 图形5：排序后的峰值幅度
        if SHOW_FIGURES && p == 1
            figure(1);
            subplot(2,3,5);
            bar(peaks_2);
            title('步骤5: 排序后的峰值 (最右两根为强信号)');
            xlabel('排序索引'); ylabel('幅度');
            grid on;
        end
        
        % 强制 N_sum = 4 (寻找4个边界)
        N_sum = 4;
        
        if length(peaks_2) >= N_sum
            % 提取最大的N_sum个峰值的位置
            % 注意：peaks_2是排好序的幅度值。我们需要找回它们在原序列中的位置(loc)
            
            % 这种查找方式在有重复值时可能有瑕疵，但在浮点数差分中通常唯一
            % 提取前4大的幅度
            top4_amps = peaks_2(end-N_sum+1:end); 
            
            loc = zeros(1, N_sum);
            for k=1:N_sum
                % 在原始peaks数组中找对应的索引
                % 注意：peaks对应indices位置
                temp_idx = find(abs(peaks - top4_amps(k)) < 1e-10, 1);
                loc(k) = indices(temp_idx);
            end
            
            % loc现在包含了4个最强峰值的【频率轴位置索引】
            % loc(1), loc(2) 对应较小的两个峰 (弱信号边界)
            % loc(3), loc(4) 对应较大的两个峰 (强信号边界)
            
            % 但上面的loc是按幅度排序提取出来的，顺序可能是乱的
            % 这里的逻辑是：
            % loc_max2: 幅度最大的两个峰的位置 (强信号)
            % loc_pass: 幅度次大的两个峰的位置 (弱信号)
            
            loc_max2 = [loc(3), loc(4)]; % 对应最大的两个幅度
            loc_pass = [loc(1), loc(2)]; % 对应较小的两个幅度
            
            % 逻辑判断核心：
            % 如果无重叠，弱信号区间 [min(pass), max(pass)] 应该完全在强信号区间 [min(max2), max(max2)] 的外面
            % 即 max(pass) < min(max2)  OR  min(pass) > max(max2)
            
            is_overlap = ~(max(loc_pass) < min(loc_max2) || min(loc_pass) > max(loc_max2));
            
            % 图形7：结果展示
            if SHOW_FIGURES && p == 1
                figure(2); clf;
                subplot(2,1,1);
                plot(f, Pxx, 'LineWidth', 1.5); hold on;
                % 标出强边界
                xline(f(loc_max2(1)), '--r', 'LineWidth', 2);
                xline(f(loc_max2(2)), '--r', 'LineWidth', 2);
                % 标出弱边界
                xline(f(loc_pass(1)), '--g', 'LineWidth', 2);
                xline(f(loc_pass(2)), '--g', 'LineWidth', 2);
                legend('功率谱', '强信号边界', '', '弱信号边界');
                title('最终检测边界可视化');
                
                subplot(2,1,2);
                axis off;
                if is_overlap
                    text(0.5, 0.5, '检测结果: 有重叠 (Overlapping)', 'Color', 'red', 'FontSize', 20, 'HorizontalAlignment', 'center');
                else
                    text(0.5, 0.5, '检测结果: 无重叠 (No Overlap)', 'Color', 'green', 'FontSize', 20, 'HorizontalAlignment', 'center');
                end
                text(0.5, 0.2, sprintf('弱信号边界: %d, %d | 强信号边界: %d, %d', ...
                    min(loc_pass), max(loc_pass), min(loc_max2), max(loc_max2)), 'HorizontalAlignment', 'center');
            end
            
            indx(p) = is_overlap;
        else
            indx(p) = 0; % 峰值不足4个，默认无法判断（或无重叠）
        end
    end
    
    fprintf('检测完成。重叠判定率: %.0f%%\n', sum(indx)/kkk*100);
end