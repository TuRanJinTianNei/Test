function [psdviaBurg, f, p] = Burg(x, Fs, varargin)
%MYBURG      使用burg算法实现的AR模型功率谱估计
% psdviaBurg 使用burg算法计算的功率谱值
% f          频率采样点
% p          模型阶数
% x          输入信号
% Fs         采样频率
% varargin   可以为数值型，即为AR模型阶数
%            可以为字符串，即为准则准则AR模型阶数由准则确定
%
% 根据输入参数类型判断
if strcmp(class(varargin{1}), 'double')
    p = varargin{1};
elseif ischar(varargin{1})
    criterion = varargin{1};
else
    error('第2个参数必须为数值型或字符串');
end
x = x(:);
N = length(x);
% 模型参数估计
if exist('p', 'var') % p变量是否存在，如果存在则不需要估计，直接使用p值
    [a, E] = computeARpara(x, p);
else % p不存在，需要估计，根据准则criterion
    p = ceil(N/3); % 阶数一般不超出信号长度的1/3
    
    % 计算1到p阶的误差
    [a, E] = computeARpara(x, p);
    
    % 计算目标函数的最小值
    kc = 1:p + 1;
    switch criterion
        case 'FPE'
            goalF = E.*(N + (kc + 1))./(N - (kc + 1));
        case 'AIC'
            goalF = N.*log(E) + 2.*kc;
    end
    [minF, p] = min(goalF); % p是目标函数最小值的位值，也即准则准则确定的阶数
    
    % 使用p值重新计算AR模型参数
    [a, E] = computeARpara(x, p);
end
[h, f] = freqz(1, a, 20e5, Fs);
psdviaBurg = E(end)*abs(h).^2./Fs;
psdviaBurg=psdviaBurg/abs(max(psdviaBurg));
psdviaBurg=(10*log10(abs(psdviaBurg)));