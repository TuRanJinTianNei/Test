%信号个数检测
%本程序测试两个信号
%-----------------------------------方法说明-------------------------------------------------------
%对接收信号（不知是否存在信号）画出频谱图
%将频谱图转换成拐点图
%将拐点图进行二次差值
%求出二次差值图的峰值点
%求峰值点中的奇异值
%将峰值点从小到大排列，取前50%的峰值点做差分，将最大的差值作为门限
%当进行计算信号数目时，从后面进行计算
%--------------------------------------------------------------------------------------------------
clc
clear all
close all
SNR=5;                   %信噪比取值，仅仿真5dB
kkk=1;                   %仅仿真1次
indx=zeros(1,kkk);

% 控制是否显示过程图形
SHOW_FIGURES = 1;        % 显示过程图形标志（1=显示，0=不显示）

fprintf('========================================\n');
fprintf('信号重叠性检测 - 单次仿真\n');
fprintf('SNR: %.0f dB\n', SNR);
fprintf('仿真次数: %d\n', kkk);
fprintf('========================================\n\n');

%% 信号2
SubCarryNN2=768;      %有效子载波数%%四倍过采样
SubCarryN2=1024;      %子载波数
ratio=1/4;           %循环前缀比例
fftLen=1024;         %FFT长度为1024
SymbN2=20;            %一帧中OFDM符号个数、每个符号传输的子载波数、6×128=768
GuardLen=SubCarryN2*ratio;    %保护时隙的长度
%%----------------------------------------发射端-----------------------------------------------------%%
%输入比特序列长度=子载波数×每载波符号数×每符号比特数10240
SignalLen2=SubCarryNN2*SymbN2*2;           %比特序列长度，256X20X2
Signal2=round(rand(1,SignalLen2));       %输出带调制的二进制比特流,0,1、round是为了产生0，1
%         ParaBitSig2=[];
ParaBitSig2=reshape(Signal2,SubCarryNN2,SymbN2*2);
%进行QPSK数据调制，将数据分为两个通道
%这是一个QPSK的调制方法
for jj=1:SymbN2
    ich22(:,jj)=ParaBitSig2(:,2*jj-1); %同相分量，ParaBitSig奇数列
    qch22(:,jj)=ParaBitSig2(:,2*jj);   %正交分量，ParaBitSig偶数列
end
kmod=1./sqrt(2);   %做一个根号二
ich02=ich22.*2-1;     %这俩变成正负一的形式，1→1，0→-1
qch02=qch22.*2-1;
ich12=ich02.*kmod;   %根号二
qch12=qch02.*kmod;
x2=ich12+qch12.*sqrt(-1);      %产生复信号,768行，6列
x2=[x2(1:SubCarryNN2/2,1:20);zeros(SubCarryN2-SubCarryNN2,20);x2(SubCarryNN2/2+1:SubCarryNN2,1:20)];
y2=ifft(x2);   %通过傅里叶反变换，将频域数据转化为时域数据，1024行，6列
ich22=real(y2);   %I信道取变换后的实部
qch22=imag(y2);   %Q信道取变换后的虚部
%插入保护间隔
%在数组中加循环前缀的方法
ich32=[ich22(fftLen-GuardLen+1:fftLen,:);ich22];%行数从1024-256+1到1024,在前面加了256个，总共1024+256=1280个,1280×6个
qch32=[qch22(fftLen-GuardLen+1:fftLen,:);qch22];
%并串变换
ich42=reshape(ich32,1,(fftLen+GuardLen)*SymbN2);%并串变换，转换成一串，变成一行，7680列
qch42=reshape(qch32,1,(fftLen+GuardLen)*SymbN2);
TrData42=ich42+qch42.*sqrt(-1);  %形成复数发射数据（加完保护时隙，保护时隙也进行相应操作）
for index=1:length(SNR)
    fprintf('  处理 SNR = %3.0f dB [%d/%d]... ', SNR(index), index, length(SNR));
    tic;
    for p=1:kkk
        %% 信号1
        SubCarryNN=256;      %有效子载波数%%四倍过采样
        SubCarryN=1024;      %子载波数
        ratio=1/4;           %循环前缀比例
        fftLen=1024;         %FFT长度为1024
        SymbN=20;            %一帧中OFDM符号个数、每个符号传输的子载波数、6×128=768
        GuardLen=SubCarryN*ratio;    %保护时隙的长度
        
        %输入比特序列长度=子载波数×每载波符号数×每符号比特数10240
        SignalLen=SubCarryNN*SymbN*2;           %比特序列长度，256X20X2
        Signal=round(rand(1,SignalLen));       %输出带调制的二进制比特流,0,1、round是为了产生0，1
        %         ParaBitSig=[];
        ParaBitSig=reshape(Signal,SubCarryNN,SymbN*2);
        %进行QPSK数据调制，将数据分为两个通道
        %这是一个QPSK的调制方法
        for j=1:SymbN
            ich(:,j)=ParaBitSig(:,2*j-1); %同相分量，ParaBitSig奇数列
            qch(:,j)=ParaBitSig(:,2*j);   %正交分量，ParaBitSig偶数列
        end
        kmod=1./sqrt(2);   %做一个根号二
        ich0=ich.*2-1;     %这俩变成正负一的形式，1→1，0→-1
        qch0=qch.*2-1;
        ich1=ich0.*kmod;   %根号二
        qch1=qch0.*kmod;
        x=ich1+qch1.*sqrt(-1);      %产生复信号,768行，6列
        x=[x(1:SubCarryNN/2,1:20);zeros(SubCarryN-SubCarryNN,20);x(SubCarryNN/2+1:SubCarryNN,1:20)];
        y=ifft(x);   %通过傅里叶反变换，将频域数据转化为时域数据，1024行，6列
        ich2=real(y);   %I信道取变换后的实部
        qch2=imag(y);   %Q信道取变换后的虚部
        %插入保护间隔
        %在数组中加循环前缀的方法
        ich3=[ich2(fftLen-GuardLen+1:fftLen,:);ich2];%行数从1024-256+1到1024,在前面加了256个，总共1024+256=1280个,1280×6个
        qch3=[qch2(fftLen-GuardLen+1:fftLen,:);qch2];
        %并串变换
        ich4=reshape(ich3,1,(fftLen+GuardLen)*SymbN);%并串变换，转换成一串，变成一行，7680列
        qch4=reshape(qch3,1,(fftLen+GuardLen)*SymbN);
        TrData=ich4+qch4.*sqrt(-1);  %形成复数发射数据（加完保护时隙，保护时隙也进行相应操作）
        ReData=awgn(TrData,SNR(index),'measured');%信道，加入高斯白噪声  %每一列是一个ofdm符号
        
        ReData42=awgn(TrData42,SNR(index),'measured');%信道，加入高斯白噪声  %每一列是一个ofdm符号
        ReData_same =ReData42.*exp(1j * (0.1) * pi * 2 .* (1:length(ReData42)));
        ReData=ReData+ReData_same;
        
        % ReData=randn(1,length(ReData));
        %ReData=awgn(ReData,SNR,'measured');%信道，加入高斯白噪声  %每一列是一个ofdm符号
        B_carrier=200e3;  % 子载波间隔
        fs=B_carrier*SubCarryN*4;         % 采样率
        fft_ofdm=fftshift(fft(real(ReData).^2));
        log_ofdm=20*log10(abs(fft_ofdm));
        NN=length(ReData42);           % 采样点个数
        n=0:NN-1;
        freq=(n/NN-0.5)*fs/10e5;  %横坐标
        % figure(1)
        % plot(freq,log_ofdm);
        % xlabel('f/MHz');
        % ylabel('幅度');
        
        %% 两个信号叠加
        %画出功率谱
        fs=1024;%采样率
        f=linspace(-fs/2,fs/2,25600);%频域横坐标，注意奈奎斯特采样定理，最大原信号最大频率不超过采样频率的一半
        f=f/(2*max(f));
        Fs=1;
        nfft=1024;
        window=hanning(150);
        noverlap=50;%每个窗口之间重叠的长度，通常取33%~50%。窗口之间重叠得越多，图像越平滑（blurred）；反之则更粗糙（blocky）
        [Pxx1,f]=pwelch(ReData,window,noverlap,length(ReData),Fs);%功率谱密度
        Pxx=10*log10(fftshift(Pxx1));
        
        % 图形1：功率谱估计
        if SHOW_FIGURES && p == 1
            figure(1);
            set(gcf,'Position',[50,50,1400,900],'Color','white');
            sgtitle('信号重叠性检测过程可视化','FontSize',16,'FontWeight','bold','Color',[0.2 0.2 0.2]);
            subplot(2,3,1);
            plot(f,Pxx,'-','LineWidth',2,'Color',[0.2 0.4 0.8]);
            xlabel('归一化频率','FontSize',11,'FontWeight','bold');
            ylabel('功率谱密度 (dB)','FontSize',11,'FontWeight','bold');
            title(sprintf('Welch法功率谱估计 (SNR=%.0f dB)', SNR(index)),'FontSize',12,'FontWeight','bold');
            grid on;
            grid minor;
            set(gca,'FontSize',10);
            box on;
        end
        
        %% 用一次拐点法求出1号边界和4号边界
        Pxx=Pxx+abs(min(Pxx));%进行平移
        Pxx_int=zeros(size(Pxx,1),1);
        for x=1:1:size(Pxx,1)
            if x==1
                Pxx_int(x)=(1/(size(Pxx,1)-1))*Pxx(x+1);
            else
                Pxx_int(x)=Pxx_int(x-1)+(1/(size(Pxx,1)-1))*Pxx(x);
            end
        end
        % 图形2：功率分布函数
        if SHOW_FIGURES && p == 1
            figure(1);
            subplot(2,3,2);
            plot(linspace(0,1,length(Pxx_int)),Pxx_int,'-','LineWidth',2,'Color',[0.85 0.2 0.2]);
            axis([0 1 0 max(Pxx_int)+1]);
            xlabel('归一化频率','FontSize',11,'FontWeight','bold');
            ylabel('累积功率','FontSize',11,'FontWeight','bold');
            title('Welch法功率分布函数','FontSize',12,'FontWeight','bold');
            grid on;
            grid minor;
            set(gca,'FontSize',10);
            box on;
        end
        y_=Pxx_int;
        x_=(linspace(0,1,length(Pxx_int)))';
        loc=[x_,y_];%坐标（第一列是横坐标，第二列是纵坐标）
        %曲线上选一点K，作垂线KA垂直于x轴，K,A,O构成直角三角形
        m1=0;%循环变量
        
        %% 求拐点的算法（差分）
        for i=4:length(y_)-4
            n1(i)=3*(loc(i,2))-(loc(i-1,2))-(loc(i-2,2))-(loc(i-3,2));
            n2(i)=loc(i+3,2)+loc(i+2,2)+loc(i+1,2)-3*(loc(i,2));
            d(i)=n2(i)-n1(i);
        end
        % 图形3：拐点检测结果（二次差分）
        if SHOW_FIGURES && p == 1
            figure(1);
            subplot(2,3,3);
            p_plot=1:length(d);
            plot(p_plot,d,'-','LineWidth',2,'Color',[0 0.7 0]);
            xlabel('采样点索引','FontSize',11,'FontWeight','bold');
            ylabel('差分值','FontSize',11,'FontWeight','bold');
            title('拐点检测结果（二次差分）','FontSize',12,'FontWeight','bold');
            grid on;
            grid minor;
            set(gca,'FontSize',10);
            box on;
        end
        x=d;
        %% VPD方法求峰值
        %%前三点均值滤波
        row_acc = abs(x);
        % row_acc = x;
        m = length(row_acc);
        row_acc1 = linspace(0,0,m);
        row_acc1(1) = row_acc(1);
        row_acc1(m) = row_acc(m);
        for i=2:m-1
            row_acc1(i)=(row_acc(i-1) + row_acc(i)+row_acc(i+1))/3;
        end
        
        for i=m-1:-1:2
            row_acc(i) = (row_acc1(i-1) + row_acc1(i)+row_acc1(i+1))/3;
        end
        
        %%找到局部最小值和局部最大值及其对应的位置，波峰点、波谷点满足：
        peaks = linspace(0,0,m);
        valleys = linspace(0,0,m);
        peakindexs = linspace(0,0,m);
        valleyindexs = linspace(0,0,m);
        peakindex = 1;
        valleyindex = 1;
        
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
        %
        %%计算VPD，VPD(n)表示第n个波谷点的值与第n个波峰点的值的差，VPD用来去掉那些假的波峰点
        pcount = peakindex-1;
        vcount = valleyindex-1;
        peakindices = linspace(0,0,pcount);
        for i = 1:pcount
            peakindices(i) = peakindexs(i);
        end
        valleyindices = linspace(0,0,vcount);
        for i = 1:vcount
            valleyindices(i) = valleyindexs(i);
        end
        if pcount>2 && vcount>2
            if peakindexs(1) < valleyindexs(1)
                peakindex=2;
            else
                peakindex=1;
            end
            vindex=1;
        end
        
        if peakindex == 2
            for i = 1:m-1
                peaks(i)=peaks(i+1);
            end
            pcount = pcount-1;
            pindex=1;
        end
        
        vpd = linspace(0,0,m);
        vpd1 = linspace(0,0,m);
        for i=1:pcount
            vpd(i) = peaks(i) - valleys(i);
        end
        
        dels = linspace(0,0, pcount);
        peakindexs1 = linspace(0,0,pcount);
        if pcount > 2
            lastcount=pcount;
            curcount = 1;
            while lastcount ~= curcount
                lastcount = curcount;
                del_count = 0;
                for i = 2:pcount-1
                    if vpd(i) <= 0.0001 * (vpd(i-1) + vpd(i)+vpd(i+1)) / 3
                        dels(i)=1;
                    end
                end
                
                count = 1;
                for i = 1:pcount
                    if dels(i) ~= 1
                        vpd1(count) = vpd(i);
                        peakindexs1(count) = peakindexs(i);
                        count = count+1;
                    else
                        del_count = del_count + 1;
                        dels(i) = 0;
                    end
                end
                pcount = pcount - del_count;
                for i = 1:pcount
                    vpd(i) = vpd1(i);
                    peakindexs(i) = peakindexs1(i);
                end
                peakindexs(pcount+1) = 0;
                vpd(pcount+1) = 0;
                
                indices = linspace(0,0,pcount);
                for i = 1:pcount
                    indices(i) = peakindexs1(i);
                end
                
                % 图形4：VPD方法检测的峰值点
                if SHOW_FIGURES && p == 1 && curcount == pcount
                    figure(1);
                    subplot(2,3,4);
                    plot(row_acc,'-','LineWidth',1.5,'Color',[0.2 0.4 0.8]);
                    hold on;
                    plot(indices, row_acc(indices),'o','MarkerFaceColor','red','MarkerEdgeColor',[0.6 0 0],...
                        'MarkerSize',10,'LineWidth',1.5);
                    if exist('valleyindices','var') && ~isempty(valleyindices)
                        valid_valleys = valleyindices(valleyindices > 0 & valleyindices <= length(row_acc));
                        if ~isempty(valid_valleys)
                            plot(valid_valleys, row_acc(valid_valleys),'s','MarkerFaceColor',[0 0.7 0],...
                                'MarkerEdgeColor',[0 0.5 0],'MarkerSize',8,'LineWidth',1.5);
                        end
                    end
                    hold off;
                    xlabel('采样点索引','FontSize',11,'FontWeight','bold');
                    ylabel('幅度','FontSize',11,'FontWeight','bold');
                    title(sprintf('VPD方法峰值检测（检测到%d个峰值）', pcount),'FontSize',12,'FontWeight','bold');
                    legend('滤波后数据','峰值点','谷值点','Location','best','FontSize',9);
                    grid on;
                    grid minor;
                    set(gca,'FontSize',10);
                    box on;
                end
                
                curcount = pcount;
            end
        end
        
        % 确保indices已定义
        if ~exist('indices','var') || isempty(indices) || length(indices) == 0
            if pcount > 0
                indices = peakindexs(1:pcount);
                indices = indices(indices > 0);
            else
                indices = [];
            end
        end
        
        %输出的峰值点
        if ~isempty(indices)
        peaks=row_acc(indices);
        else
            peaks = [];
        end
        % figure(6)
        % h=histogram(peaks,100);
        M=mean(x);
        s=std(x);%标准差
        H=M+3*s;
        if ~isempty(peaks)
        peaks_dou=abs(peaks);
        b=[];
        [D, ia] = setdiff(peaks_dou, b);
        peaks_2=peaks(ia);%从小到大排列的峰值点
        else
            peaks_2 = [];
        end
        
        % 图形5：峰值点排序结果
        if SHOW_FIGURES && p == 1
            figure(1);
            subplot(2,3,5);
            if ~isempty(peaks_2)
                scatter(1:length(peaks_2),peaks_2,'filled','SizeData',80,'MarkerFaceColor',[0.85 0.2 0.2],...
                    'MarkerEdgeColor',[0.6 0 0],'LineWidth',1.5);
                xlabel('峰值点索引（排序后）','FontSize',11,'FontWeight','bold');
                ylabel('峰值幅度','FontSize',11,'FontWeight','bold');
                title(sprintf('峰值点排序结果（共%d个峰值点）', length(peaks_2)),'FontSize',12,'FontWeight','bold');
            else
                axis([0 1 0 1]);
                text(0.5, 0.5, '未检测到峰值点', ...
                    'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold', ...
                    'Color', [0.85 0.2 0.2]);
                title('峰值点排序结果','FontSize',12,'FontWeight','bold');
            end
            grid on;
            grid minor;
            set(gca,'FontSize',10);
            box on;
        end
        
        %再做差分
        if length(peaks_2) >= 2
        for ii=1:length(peaks_2)-1
            dd(ii)=peaks_2(ii+1)-peaks_2(ii);%dd是差分结果，dd横坐标如果是118，则118:length(peaks_2)的数据都是离群点
        end
            
            % 图形6：峰值点差分结果
            if SHOW_FIGURES && p == 1
                figure(1);
                subplot(2,3,6);
                pp=1:length(dd);
                plot(pp,dd,'-','LineWidth',2,'Color',[0.8 0.2 0.8]);
                hold on;
                MM=max(dd);
                [~,loc1]=find(dd==MM);
                plot(loc1, MM,'o','MarkerFaceColor','red','MarkerEdgeColor',[0.6 0 0],...
                    'MarkerSize',12,'LineWidth',2);
                Hr=MM/10;
                yline(Hr,'--','LineWidth',2,'Color',[0 0.7 0],'Label','门限 Hr',...
                    'LabelHorizontalAlignment','left','FontSize',10,'FontWeight','bold');
                hold off;
                xlabel('差分点索引','FontSize',11,'FontWeight','bold');
                ylabel('差分值','FontSize',11,'FontWeight','bold');
                title(sprintf('峰值点差分结果（最大差分点：%d）', loc1),'FontSize',12,'FontWeight','bold');
                legend('差分值','最大差分点','门限','Location','best','FontSize',9);
                grid on;
                grid minor;
                set(gca,'FontSize',10);
                box on;
            end
        
        L=length(peaks_2);%峰值点的个数
        MM=max(dd);
        Hr=MM/10;
            [~,loc1]=find(dd==MM);%找到最大差分点处
        dd=dd(1:loc1-1);
        Num=L-loc1;%在最大差分点外的数据个数
        N_sum=4;
        else
            L = 0;
            MM = 0;
            Hr = 0;
            loc1 = 0;
            Num = 0;
            N_sum = 4;
            if SHOW_FIGURES && p == 1
                figure(1);
                subplot(2,3,6);
                axis([0 1 0 1]);
                text(0.5, 0.5, '峰值点不足，无法进行差分分析', ...
                    'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold', ...
                    'Color', [0.85 0.2 0.2]);
                title('峰值点差分结果','FontSize',12,'FontWeight','bold');
                set(gca,'FontSize',10);
            end
        end
%         b=mod(Num,2);
%         M_indx=MM;
%         while b==0&&M_indx>Hr
%             N_sum=N_sum+Num;
%             M_indx=max(dd);
%             [num,loc]=find(dd==M_indx);%找到最大差分点处
%             dd=dd(1:loc-1);
%             Num=L-N_sum-loc;
%             b=mod(Num,2);
%             
%         end
%         %计算信号数目
%         N_signal=N_sum/2;
        if length(peaks_2) >= N_sum
            k=N_sum;
            loc = zeros(1, N_sum);
            peaks_2_temp = peaks_2;  % 保存副本用于图形显示
            while k
                c=max(peaks_2);
                [~,loc(k)]=find(peaks==c);
                [~, a] = setdiff(peaks_2, c);
                peaks_2=peaks_2(a);
                k=k-1;
            end
            loc_max2=[loc(N_sum),loc(N_sum-1)];
            [~, aa]=setdiff(loc,loc_max2);
            loc_pass=loc(aa);
            
            % 图形7：重叠性检测结果总结
            if SHOW_FIGURES && p == 1
                figure(2);
                set(gcf,'Position',[150,150,1200,800],'Color','white');
                
                % 子图1：功率谱上标注检测结果
                subplot(2,2,1);
                plot(f,Pxx,'-','LineWidth',2,'Color',[0.2 0.4 0.8]);
                xlabel('归一化频率','FontSize',12,'FontWeight','bold');
                ylabel('功率谱密度 (dB)','FontSize',12,'FontWeight','bold');
                title('功率谱','FontSize',13,'FontWeight','bold');
                grid on;
                grid minor;
                set(gca,'FontSize',11);
                box on;
                
                % 子图2：峰值点分布和重叠判断
                subplot(2,2,2);
                if ~isempty(peaks_2_temp)
                    scatter(1:length(peaks_2_temp),peaks_2_temp,'filled','SizeData',100,...
                        'MarkerFaceColor',[0.85 0.2 0.2],'MarkerEdgeColor',[0.6 0 0],'LineWidth',1.5);
                    hold on;
                    % 标注前4个最大峰值的位置
                    if length(peaks_2_temp) >= N_sum
                        for i = 1:N_sum
                            idx_in_peaks2 = find(peaks_2_temp == peaks(loc(i)), 1);
                            if ~isempty(idx_in_peaks2)
                                plot(idx_in_peaks2, peaks_2_temp(idx_in_peaks2), 'o', ...
                                    'MarkerSize', 15, 'MarkerFaceColor', 'yellow', ...
                                    'MarkerEdgeColor', 'black', 'LineWidth', 2);
                            end
                        end
                    end
                    hold off;
                end
                xlabel('峰值点索引','FontSize',12,'FontWeight','bold');
                ylabel('峰值幅度','FontSize',12,'FontWeight','bold');
                title('峰值点分布（标注前4个最大峰值）','FontSize',13,'FontWeight','bold');
                grid on;
                grid minor;
                set(gca,'FontSize',11);
                box on;
                
                % 子图3：重叠性判断示意图
                subplot(2,2,3);
                axis([0 1 0 1]);
                axis off;
                
                % 判断重叠性
                if max(loc_pass)<min(loc_max2) || min(loc_pass)>max(loc_max2)
                    result_str = '无重叠';
                    result_color = [0 0.7 0];  % 绿色
                    result_icon = '✓';
                else
                    result_str = '有重叠';
                    result_color = [0.85 0.2 0.2];  % 红色
                    result_icon = '✗';
                end
                
                % 绘制背景框
                rectangle('Position',[0.05 0.05 0.9 0.9],'FaceColor',[0.98 0.98 0.98],...
                    'EdgeColor',[0.5 0.5 0.5],'LineWidth',2);
                
                % 绘制标题
                text(0.5, 0.95, '重叠性判断', 'FontSize', 14, 'FontWeight', 'bold', ...
                    'HorizontalAlignment', 'center', 'Color', [0.2 0.2 0.2]);
                
                % 绘制检测结果（高亮显示）
                text(0.5, 0.7, sprintf('%s %s', result_icon, result_str), 'FontSize', 16, ...
                    'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Color', result_color);
                
                % 绘制详细信息
                % 先判断条件，然后构建字符串数组
                if max(loc_pass)<min(loc_max2) || min(loc_pass)>max(loc_max2)
                    condition_str1 = '  max(loc_pass) < min(loc_max2) 或';
                    condition_str2 = '  min(loc_pass) > max(loc_max2)';
                    text_str = {
                        sprintf('SNR: %.0f dB', SNR(index));
                        sprintf('峰值点总数: %d', L);
                        sprintf('前4个最大峰值位置:');
                        sprintf('  loc(1) = %d', loc(1));
                        sprintf('  loc(2) = %d', loc(2));
                        sprintf('  loc(3) = %d', loc(3));
                        sprintf('  loc(4) = %d', loc(4));
                        '';
                        sprintf('两个主峰值位置:');
                        sprintf('  loc_max2 = [%d, %d]', loc_max2(1), loc_max2(2));
                        sprintf('其他峰值位置:');
                        sprintf('  loc_pass = [%s]', num2str(loc_pass));
                        '';
                        sprintf('判断条件:');
                        condition_str1;
                        condition_str2;
                        '';
                        sprintf('检测结果: %s', result_str);
                    };
                else
                    condition_str1 = '  loc_pass 与 loc_max2 有交集';
                    text_str = {
                        sprintf('SNR: %.0f dB', SNR(index));
                        sprintf('峰值点总数: %d', L);
                        sprintf('前4个最大峰值位置:');
                        sprintf('  loc(1) = %d', loc(1));
                        sprintf('  loc(2) = %d', loc(2));
                        sprintf('  loc(3) = %d', loc(3));
                        sprintf('  loc(4) = %d', loc(4));
                        '';
                        sprintf('两个主峰值位置:');
                        sprintf('  loc_max2 = [%d, %d]', loc_max2(1), loc_max2(2));
                        sprintf('其他峰值位置:');
                        sprintf('  loc_pass = [%s]', num2str(loc_pass));
                        '';
                        sprintf('判断条件:');
                        condition_str1;
                        '';
                        sprintf('检测结果: %s', result_str);
                    };
                end
                text(0.1, 0.55, text_str, 'FontSize', 10, 'VerticalAlignment', 'top', ...
                     'FontName', 'FixedWidth', 'Color', [0.2 0.2 0.2]);
                
                % 子图4：检测结果总结
                subplot(2,2,4);
                axis([0 1 0 1]);
                axis off;
                
                % 绘制背景框
                rectangle('Position',[0.05 0.05 0.9 0.9],'FaceColor',[0.98 0.98 0.98],...
                    'EdgeColor',[0.5 0.5 0.5],'LineWidth',2);
                
                % 绘制标题
                text(0.5, 0.95, '检测结果信息', 'FontSize', 14, 'FontWeight', 'bold', ...
                    'HorizontalAlignment', 'center', 'Color', [0.2 0.2 0.2]);
                
                % 绘制检测结果（高亮显示）
                text(0.5, 0.65, sprintf('%s %s', result_icon, result_str), 'FontSize', 16, ...
                    'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Color', result_color);
                
                % 绘制详细信息
                text_str2 = {
                    sprintf('SNR: %.0f dB', SNR(index));
                    sprintf('峰值点总数: %d', L);
                    sprintf('最大差分位置: %d', loc1);
                    sprintf('最大差分值: %.4f', MM);
                    sprintf('门限 Hr: %.4f', Hr);
                    '';
                    sprintf('算法参数:');
                    sprintf('  N_sum = %d', N_sum);
                    sprintf('  提取前%d个最大峰值', N_sum);
                    '';
                    sprintf('检测结果: %s', result_str);
                };
                text(0.1, 0.45, text_str2, 'FontSize', 11, 'VerticalAlignment', 'top', ...
                     'FontName', 'FixedWidth', 'Color', [0.2 0.2 0.2]);
                
                sgtitle(sprintf('信号重叠性检测结果 (SNR=%.0f dB)', SNR(index)), ...
                        'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.2 0.2 0.2]);
            end
            
            if max(loc_pass)<min(loc_max2) || min(loc_pass)>max(loc_max2)
                indx(p)=0;
                %     disp("无重叠")
            else
                indx(p)=1;
                %     disp("有重叠")
            end
            else
                indx(p)=0;  % 如果峰值点不足，默认无重叠
                if SHOW_FIGURES && p == 1
                    fprintf('警告：峰值点不足，无法进行重叠性判断\n');
                end
        end
    end
    sum_indx=sum(indx);
    corr(index)=sum_indx./(kkk);
    elapsed_time = toc;
    fprintf('完成 | 正确率: %5.2f%% (%d/%d) | 耗时: %.2f秒\n', ...
        corr(index)*100, sum_indx, kkk, elapsed_time);
end
fprintf('----------------------------------------\n');
fprintf('仿真完成！正在绘制结果曲线...\n');
fprintf('========================================\n');

%% 绘制正确率曲线
figure;
plot(SNR,corr,'r--d');
grid off;
% axis([-5 12 10^-8 10^0])
xlabel('SNR');
ylabel('正确率');
legend('信号重叠性检测正确率');


