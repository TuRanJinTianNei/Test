% %% ----------------------------------------参数设置----------------------------------------------------- %%
% subCarry = 256;             %有效子载波数%%四倍过采样
% fftLen = 1024;              %FFT长度为1024
% symbN = 10;                 %一帧中OFDM符号个数
% cpLen = 256;                %保护时隙的长度
% SNR = 10;                   %信噪比取值，单位dB
% csLen = 20;                 %循环后缀
% k = 2;                      %调制阶数
function ReData = generate_OFDM(subCarry, fftLen, cpLen, csLen, symbN, SNR, k)
    %% ----------------------------------------发射端----------------------------------------------------- %%
    %输入比特序列长度=子载波数×每载波符号数×每符号比特数10240
    sigLen = subCarry*symbN*k;           %比特序列长度，256*20*2
    sig = randi([0 1], 1, sigLen);       %输出待调制的二进制比特流,0,1
    sig = sig.*sqrt(2) - 1/sqrt(2);
    ParaBitSig = reshape(sig, subCarry, symbN*2);
    x = ParaBitSig(:,1:2:end) + ParaBitSig(:,2:2:end).*1i;                %产生复信号
    x = [x(1:subCarry/2, :); zeros(fftLen - subCarry, symbN); x(subCarry/2 + 1:subCarry, :)]; %中间补零
%     x = [x(1:subCarry/2, :); x(subCarry/2 + 1:subCarry, :); zeros(fftLen - subCarry, symbN)]; %后面补零
    y = ifft(x);      %通过傅里叶反变换，将频域数据转化为时域数据
    %插入保护间隔
%     csLen = (fftLen + cpLen)*beta;
    y = [y(fftLen - cpLen + 1:fftLen,:); y; y(1:csLen,:)];  %引入循环前缀和后缀

    % 加窗
    window = repmat(rcoswindow(csLen, fftLen + cpLen)', 1, symbN);
    Tx_ps = y.*window;
    Tx_ss = zeros(1,symbN*(fftLen + cpLen) + csLen);
    for i = 1:symbN
        Tx_ss((fftLen + cpLen)*(i - 1) + 1:(fftLen + cpLen)*i + csLen) = Tx_ss((fftLen + cpLen)*(i - 1) + 1:(fftLen + cpLen)*i + csLen) + Tx_ps(:,i).';
    end
    % ReData = Tx_ss;
    ReData = awgn(Tx_ss, SNR, 'measured');%信道，加入高斯白噪声  %每一列是一个ofdm符号
end
function window = rcoswindow(M, N)
    window = zeros(1, N + M);
    for i = 1: M
        window(i) = 0.5 + 0.5*cos(pi + (i - 1)*pi/(M));
    end
    window(M + 1: N)=1;
    for i = N + 1: N + M
        window(i) = 0.5 + 0.5*cos((i - N)*pi/(M));
    end
end