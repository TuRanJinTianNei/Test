% ============================================================================
% 文件名: rcoswindow.m
% 功能: 升余弦窗函数生成器
% 描述:
%   生成升余弦（Raised Cosine）窗函数，用于OFDM符号的平滑加窗
%   窗函数在符号起始和结束部分提供平滑过渡，降低频谱旁瓣
% 输入参数:
%   beta - 滚降系数（过渡带占比，通常为1/16或1/32等）
%   Ts   - 窗函数覆盖长度（通常为包含循环前缀的OFDM符号长度 N+GI）
% 输出参数:
%   rcosw - 升余弦窗向量（列向量），长度为 round((1+beta)*Ts)+1
% 说明:
%   窗函数在[0, beta*Ts]和[Ts, (1+beta)*Ts]区间为过渡段（余弦上升/下降），
%   在[beta*Ts, Ts]区间为平坦段（值为1），用于平滑符号边界
% ============================================================================
function[rcosw]=rcoswindow(beta,Ts)

%------------------------------------------------------------------------------
% 【算法流程：生成升余弦窗函数】
%------------------------------------------------------------------------------
% 窗函数总长度：round((1+beta)*Ts)+1，分为三段

% 步骤1: 初始化参数（确保所有长度都是整数）
Ts_int = round(Ts);  % 确保 Ts 为整数
beta_Ts_int = round(beta*Ts_int);  % 确保 beta*Ts 为整数
win_length = round((1+beta)*Ts_int);  % 确保总长度为整数

% 步骤2: 初始化时间轴和窗函数数组
t = 0:win_length;  % 时间轴，长度为 win_length+1（索引 1 到 win_length+1）
rcosw = zeros(1, win_length+1);  % 窗函数数组，长度为 win_length+1（索引 1 到 win_length+1）

% 步骤3: 左过渡段 [0, beta*Ts]：从0平滑上升到1
% 公式：w(t) = 0.5 + 0.5*cos(π + t*π/(beta*Ts))
for i = 1:beta_Ts_int
    if i <= length(t) && beta_Ts_int > 0
        rcosw(i) = 0.5 + 0.5*cos(pi + t(i)*pi/(beta*Ts_int));
    end
end

% 步骤4: 平坦段 [beta*Ts, Ts]：保持为1（不衰减）
for i = (beta_Ts_int+1):Ts_int
    if i <= length(rcosw)
        rcosw(i) = 1;
    end
end

% 步骤5: 右过渡段 [Ts, (1+beta)*Ts]：从1平滑下降到0
% 公式：w(t) = 0.5 + 0.5*cos((t-Ts)*π/(beta*Ts))
for i = (Ts_int+1):(win_length+1)
    if i <= length(t) && beta_Ts_int > 0
        rcosw(i) = 0.5 + 0.5*cos((t(i)-Ts_int)*pi/(beta*Ts_int));
    end
end

% 步骤6: 转换为列向量（便于矩阵运算）
rcosw = rcosw';
