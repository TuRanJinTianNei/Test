%% 函数
function [cumulants] = func_get_cumulants(signalData)
shape = size(signalData);
H = shape(1);%行
N = shape(2);%列s
cumulants = zeros(H, N);

for row = 1:H
    M20 = sum(signalData(row,:).^2)/N;
%     M21 = sum(abs(signalData(row,:)).^2)/N;%实际上是信号的能量
    M21 = sum(abs(signalData(row,:)).^2)/N;%实际上是信号的能量
    M22 = sum(conj((signalData(row,:))).^2)/N;
    M40 = sum(signalData(row,:).^4)/N;
    M41 = sum(abs(signalData(row,:)).^2.*signalData(row,:).^2)/N;
    M42 = sum(abs(signalData(row,:)).^4)/N;
    M43 = sum(abs(signalData(row,:)).^2.*conj((signalData(row,:))).^2)/N;
    M60 = sum(signalData(row,:).^6)/N;
    M61 = sum(abs(signalData(row,:)).^2.*signalData(row,:).^4)/N;
    M62 = sum(abs(signalData(row,:)).^4.*signalData(row,:).^2)/N;
    M63 = sum(abs(signalData(row,:)).^6)/N;
    M84 = sum(abs(signalData(row,:)).^8)/N;
    m20 = 2*M42/(M21^2);
    m30 = 6*M63/(M21^3);
    m40 = 24*M84/(M21^4);
    C20 = M20;
    C21 = M21;
    C40 = M40 - 3*M20^2;
    C41 = M41 - 3*M20*M21;
    C42 = M42 - abs(M20)^2 - 2*M21^2;
    C60 = M60 - 15*M20*M40 + 3*M20^3;
    C61 = M61 - 5*M21*M40 - 10*M20*M41 + 30*M20^2*M21;
    C62 = M62 - 6*M20*M42 - 8*M21*M41 - M22*M40 + 6*M20^2*M22 + 24*M21^2*M20;
    C63 = M63 -9*M21*M42 + 12*M21^3 - 3*M20*M43 - 3*M22*M41 + 18*M20*M21*M22;
    C4two=C42/((C21)^2);
    F=(m40-16*m30+18*m20)/(m20-2)^2;
    cumulants(row, 1) = M84;
    cumulants(row, 2) = M63;
    cumulants(row, 3) = M42;
    cumulants(row, 4) = M21;
    cumulants(row, 5) = M20;
    cumulants(row, 6) = M41;
    cumulants(row, 7) = C21;
    cumulants(row, 8) = abs(M84);
    cumulants(row, 9) = abs(C63);
    cumulants(row, 10) = abs(C4two);
end
end