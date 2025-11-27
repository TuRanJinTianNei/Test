function rchan = jakesChannel()
%% Jakes模型瑞利衰落信道创建函数
% 功能：创建并返回一个配置好的瑞利衰落信道对象
% 
% 输出参数：
%   rchan - comm.RayleighChannel对象，配置好的瑞利衰落信道
%
% 信道参数（三径多径配置）：
%   - 采样率：15.36 MHz（与test6.m的OFDM系统采样率一致）
%   - 路径延迟：[0, 2e-6, 4e-6] 秒（三径，时延在2-4微秒范围内）
%   - 路径功率：[0, -3, -6] dB（主径功率最高，多径功率递减）
%   - 最大多普勒频移：100 Hz（与test6.m的fd参数一致）
%   - 散射体数量：34（N0参数，用于Jakes模型内部实现）
%   - 相干时间：约1.59 ms（1/(2π*fd)）
%
% 说明：
%   本函数使用MATLAB的comm.RayleighChannel对象实现Jakes模型
%   三径配置模拟频率选择性衰落（多径传播）
%   时延设置为2-4微秒，小于CP时长（4.6875微秒），避免符号间干扰
%
% 使用示例：
%   rchan = jakesChannel();  % 创建信道对象
%   rx_signal = rchan(tx_signal);  % 通过信道传输信号
%   rchan.reset();  % 重置信道状态

%% 信道参数设置（三径多径配置）
fs = 15.36e6;            % Sample rate (Hz) - 与test6.m的OFDM系统采样率一致
pathDelays = [0, 0.5e-6, 1e-6];  % Path delays (s) - 三径：0秒（主径），2微秒，4微秒（多径）
pathPower = [0, -5, -10];       % Path power (dB) - 三径功率：0 dB（主径），-3 dB，-6 dB（多径）
fD = 100;                % Maximum Doppler shift (Hz) - 与test6.m的fd参数一致

%% 创建瑞利信道对象
rchan = comm.RayleighChannel('SampleRate', fs, ...
    'PathDelays', pathDelays, ...
    'AveragePathGains', pathPower, ...
    'MaximumDopplerShift', fD);

end

