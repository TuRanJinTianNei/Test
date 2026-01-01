%信号存在性检测
%本程序用蒙特卡洛仿真，在信噪比为-10dB到10dB的范围内计算信号检测概率
%该程序直接运行，只存在一个信号，检测准确率高
%-----------------------------------方法说明-------------------------------------------------------
%对接收信号（不知是否存在信号）画出频谱图
%将频谱图转换成拐点图
%将拐点图进行二次差值
%求出二次差值图的峰值点
%求峰值点中的奇异值
%若存在至少两个奇异值大于门限，则证明存在信号
%--------------------------------------------------------------------------------------------------
clc
clear all
close all
SNR=0;                   %信噪比取值，仅仿真0dB
kkk=1;                    %仅仿真1次

% 控制是否显示过程图形
SHOW_FIGURES = 1;        % 显示过程图形标志（1=显示，0=不显示）

fprintf('========================================\n');
fprintf('信号存在性检测 - 单次仿真\n');
fprintf('SNR: %.0f dB\n', SNR);
fprintf('仿真次数: %d\n', kkk);
fprintf('========================================\n\n');

for index=1:length(SNR)
    fprintf('正在处理 SNR = %.0f dB...\n', SNR(index));
    tic;  % 开始计时
    for p=1:kkk
        fprintf('  开始仿真...\n');
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
    ParaBitSig=[];
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
    x=[x(1:128,1:20);zeros(768,20);x(129:256,1:20)];
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
    %         TrData=zeros(1,25600);
    ReData=awgn(TrData,SNR(index),'measured');%信道，加入高斯白噪声  %每一列是一个ofdm符号

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
        sgtitle('信号存在性检测过程可视化','FontSize',16,'FontWeight','bold','Color',[0.2 0.2 0.2]);
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
            %     Pxx_int(x)=Pxx(x+1);
        else
            Pxx_int(x)=Pxx_int(x-1)+(1/(size(Pxx,1)-1))*Pxx(x);
            %     Pxx_int(x)=Pxx_int(x-1)+Pxx(x);
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
    %% VPD方法
    %%前三点均值滤波
    row_acc = abs(x);
    % row_acc = x;
    m = length(row_acc);%数据的长度
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
    % figure;
    % plot(row_acc,'-o', 'MarkerIndices',peakindices,'MarkerFaceColor','red','MarkerSize',5);
    %
    % figure;
    % plot(row_acc,'-s', 'MarkerIndices',valleyindices,'MarkerFaceColor','green','MarkerSize',5);
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
    %输出的峰值点
        peaks=row_acc(indices);
        M=mean(x);
        s=std(x);%标准差
        H=M+3*s;
        peaks_dou=abs(peaks);
        a1=find(peaks_dou>H);
        h=length(a1);
        
        % 图形5：峰值点统计和门限判断
        if SHOW_FIGURES && p == 1
            figure(1);
            subplot(2,3,5);
            histogram(peaks_dou, 30, 'FaceColor',[0.3 0.6 0.9],'EdgeColor','none');
            hold on;
            xline(H,'--','LineWidth',2.5,'Color',[0.85 0.2 0.2],'Label','门限 H=M+3σ',...
                'LabelHorizontalAlignment','left','FontSize',10,'FontWeight','bold');
            hold off;
            xlabel('峰值幅度','FontSize',11,'FontWeight','bold');
            ylabel('频数','FontSize',11,'FontWeight','bold');
            title(sprintf('峰值点分布（门限=%.4f）', H),'FontSize',12,'FontWeight','bold');
            legend('峰值分布','门限','Location','best','FontSize',9);
            grid on;
            grid minor;
            set(gca,'FontSize',10);
            box on;
            
            % 图形6：检测结果总结
            subplot(2,3,6);
            axis([0 1 0 1]);
            axis off;
            
            % 判断检测结果
            if h > 1
                result_str = '检测到信号';
                result_color = [0 0.7 0];  % 绿色
                result_icon = '✓';
            else
                result_str = '未检测到信号';
                result_color = [0.85 0.2 0.2];  % 红色
                result_icon = '✗';
            end
            
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
            text_str = {
                sprintf('SNR: %.0f dB', SNR(index));
                sprintf('峰值点总数: %d', length(peaks_dou));
                sprintf('大于门限的峰值数: %d', h);
                sprintf('门限值 H: %.4f', H);
                sprintf('均值 M: %.4f', M);
                sprintf('标准差 σ: %.4f', s);
                '';
                sprintf('判断条件: h > 1');
                sprintf('检测结果: %s', result_str);
            };
            text(0.1, 0.45, text_str, 'FontSize', 11, 'VerticalAlignment', 'top', ...
                 'FontName', 'FixedWidth', 'Color', [0.2 0.2 0.2]);
        end
        
        if h>1
            indx(p)=1;
        else
            indx(p)=0;
        end
        
    end
    elapsed_time = toc;  % 计算耗时
    corr(index)=sum(indx)./kkk;
    % 判断检测结果
    if indx(1) == 1
        result_str = '检测到信号';
    else
        result_str = '未检测到信号';
    end
    fprintf('  完成! 检测结果: %s 耗时: %.2f秒\n', result_str, elapsed_time);
    fprintf('  准确率: %.2f%% (%d/%d)\n\n', corr(index)*100, sum(indx), kkk);
end

fprintf('========================================\n');
fprintf('仿真完成！\n');
fprintf('========================================\n');
%% 绘制准确率曲线
if SHOW_FIGURES
    figure(2);
    set(gcf,'Position',[200,200,800,600],'Color','white');
else
    figure;
end
plot(SNR,corr,'--d','LineWidth',3,'MarkerSize',12,'Color',[0.85 0.2 0.2],...
    'MarkerFaceColor',[0.85 0.2 0.2],'MarkerEdgeColor',[0.6 0 0]);
grid on;
grid minor;
xlabel('信噪比 SNR (dB)','FontSize',14,'FontWeight','bold');
ylabel('检测准确率','FontSize',14,'FontWeight','bold');
title('信号存在性检测准确率','FontSize',16,'FontWeight','bold');
legend('信号检测准确率','Location','best','FontSize',12,'FontWeight','bold');
axis([min(SNR)-0.5 max(SNR)+0.5 0 1.1]);
set(gca,'FontSize',12,'LineWidth',1.5);
box on;
yline(1,'--','LineWidth',1.5,'Color',[0.5 0.5 0.5],'Alpha',0.5);
yline(0.5,'--','LineWidth',1,'Color',[0.7 0.7 0.7],'Alpha',0.3);