% subCarry = 256;             %有效子载波数%%四倍过采样
% fftLen = 1024;              %FFT长度为1024
% symbN1 = 13;                 %一帧中OFDM符号个数
% symbN2 = 1;
% cpLen1 = 72;               %保护时隙的长度
% cpLen2 = 88;
% SNR = 10;                   %信噪比取值，单位dB
% csLen = 20;                 %循环后缀
% k = 2;                      %调制阶数
function ReData = generate_OFDM2(subCarry, fftLen, cpLen1, cpLen2, csLen1, csLen2, symbN1, symbN2, symbN, SNR, k)
    L = symbN1 * (fftLen + cpLen1) + symbN2 * (fftLen + cpLen2) + csLen2;
    ReData = zeros(symbN, L);
    for i = 1: symbN
        ReData(i, :) = generate_OFDM_frame(subCarry, fftLen, cpLen1, cpLen2, csLen1, csLen2, symbN1, symbN2, SNR, k);
    end
    ReData2 = zeros(1, L*symbN - (symbN-1)*csLen2);
    for i = 1:symbN
        ReData2((L - csLen2)*(i - 1) + 1:(L - csLen2)*i + csLen2) = ReData2((L - csLen2)*(i - 1) + 1:(L - csLen2)*i + csLen2) + ReData(i,:);
    end
    ReData = awgn(ReData2, SNR, 'measured');
%     ReData = ReData - ReData2;
end

function data_frame = generate_OFDM_frame(subCarry, fftLen, cpLen1, cpLen2, csLen1, csLen2, symbN1, symbN2, SNR, k)
    data1 = generate_OFDM3(subCarry, fftLen, cpLen1, csLen1, symbN1, SNR, k);
    data2 = generate_OFDM3(subCarry, fftLen, cpLen2, csLen2, symbN2, SNR, k);
    L1 = symbN1 * (fftLen + cpLen1) + csLen1;
    L2 = symbN2 * (fftLen + cpLen2) + csLen2;
    sig1 = [data1, zeros(1, L2 - csLen1)];
    sig2 = [zeros(1, L1 - csLen1), data2];
    data_frame = sig1 + sig2;
end