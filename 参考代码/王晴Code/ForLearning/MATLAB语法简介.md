# MATLAB 语法简介
## 以 `indentify_signal_exist.m` 为例

本文档通过分析 `indentify_signal_exist.m` 文件，简要介绍 MATLAB 的基本语法。

---

## 1. 注释

MATLAB 使用 `%` 符号表示单行注释，`%%` 表示代码节（Section）。

```matlab
%信号存在性检测                    % 单行注释
%本程序用蒙特卡洛仿真...

%% 信号1                          % 代码节，用于组织代码结构
SubCarryNN=256;                  % 行内注释
```

---

## 2. 变量和赋值

MATLAB 是动态类型语言，变量无需声明即可使用。

```matlab
SNR=-10:1:0;                     % 创建数组：从-10到0，步长为1
kkk=100;                         % 标量赋值
SubCarryNN=256;                  % 整数赋值
ratio=1/4;                       % 小数赋值
```

**注意**：
- 变量名区分大小写
- 变量名可以包含字母、数字和下划线，但不能以数字开头
- `=` 用于赋值

---

## 3. 数组和矩阵操作

### 3.1 创建数组

```matlab
SNR=-10:1:0;                     % 创建向量：[ -10, -9, -8, ..., 0 ]
f=linspace(-fs/2,fs/2,25600);   % 线性间隔数组：从-fs/2到fs/2，共25600个点
zeros(768,20);                   % 创建768×20的零矩阵
ones(1,100);                     % 创建1×100的全1向量
```

### 3.2 数组索引

```matlab
x(1:128,1:20)                    % 提取第1-128行，第1-20列
x(129:256,1:20)                  % 提取第129-256行，第1-20列
Pxx(x+1)                         % 访问数组元素（注意：MATLAB索引从1开始）
loc(i,2)                          % 访问矩阵第i行第2列的元素
```

**重要**：MATLAB 数组索引从 **1** 开始（不是0）！

### 3.3 数组拼接

```matlab
x=[x(1:128,1:20);zeros(768,20);x(129:256,1:20)];  % 垂直拼接（用分号）
loc=[x_,y_];                     % 水平拼接（用空格或逗号）
```

### 3.4 数组形状变换

```matlab
reshape(Signal,SubCarryNN,SymbN*2);  % 将数组重塑为指定维度
```

---

## 4. 数学运算

### 4.1 基本运算

```matlab
GuardLen=SubCarryN*ratio;        % 乘法
SignalLen=SubCarryNN*SymbN*2;    % 乘法运算
kmod=1./sqrt(2);                 % 除法（注意：`./` 是元素级除法）
```

### 4.2 元素级运算（点运算）

```matlab
ich0=ich.*2-1;                   % `.*` 表示元素级乘法
qch1=qch0.*kmod;                 % 每个元素分别相乘
x=ich1+qch1.*sqrt(-1);           % 复数运算
```

**区别**：
- `*`：矩阵乘法
- `.*`：元素级乘法（对应元素相乘）

### 4.3 数学函数

```matlab
sqrt(2)                          % 平方根
abs(x)                            % 绝对值
real(y)                           % 取实部
imag(y)                           % 取虚部
mean(x)                           % 平均值
std(x)                            % 标准差
max(f)                            % 最大值
min(Pxx)                          % 最小值
length(ReData)                    % 数组长度
size(Pxx,1)                       % 矩阵维度（第1维的大小）
```

---

## 5. 控制结构

### 5.1 for 循环

```matlab
for index=1:length(SNR)          % 从1到SNR的长度
    % 循环体
end

for p=1:kkk                      % 从1到kkk
    % 循环体
end

for i=4:length(y_)-4            % 从4到length(y_)-4
    % 循环体
end

for i=m-1:-1:2                   % 从m-1递减到2，步长为-1
    % 循环体
end
```

### 5.2 while 循环

```matlab
while lastcount ~= curcount       % `~=` 表示不等于
    % 循环体
    curcount = pcount;
end
```

### 5.3 if-else 条件语句

```matlab
if x==1                           % `==` 表示等于判断
    Pxx_int(x)=(1/(size(Pxx,1)-1))*Pxx(x+1);
else
    Pxx_int(x)=Pxx_int(x-1)+(1/(size(Pxx,1)-1))*Pxx(x);
end

if row_acc(i) >row_acc(i-1) && row_acc(i)>=row_acc(i+1)  % `&&` 表示逻辑与
    peaks(peakindex)=row_acc(i);
    peakindex = peakindex+1;
end

if h>1                            % 条件判断
    indx(p)=1;
else
    indx(p)=0;
end
```

**逻辑运算符**：
- `&&`：逻辑与
- `||`：逻辑或
- `~`：逻辑非
- `==`：等于
- `~=`：不等于
- `>`、`<`、`>=`、`<=`：比较运算符

---

## 6. 函数调用

### 6.1 内置函数

```matlab
clc                              % 清空命令窗口
clear all                        % 清除所有变量
close all                        % 关闭所有图形窗口

rand(1,SignalLen)               % 生成随机数
round(...)                       % 四舍五入
ifft(x)                          % 逆快速傅里叶变换
pwelch(ReData,window,noverlap,length(ReData),Fs)  % Welch功率谱估计
fftshift(Pxx1)                   % 将零频移到中心
log10(...)                       % 以10为底的对数
find(peaks_dou>H)                % 查找满足条件的索引
sum(indx)                        % 求和
```

### 6.2 工具箱函数

```matlab
awgn(TrData,SNR(index),'measured')  % 添加高斯白噪声（需要Communications Toolbox）
hanning(150)                        % 生成汉宁窗（需要Signal Processing Toolbox）
```

---

## 7. 复数运算

MATLAB 原生支持复数运算。

```matlab
sqrt(-1)                         % 虚数单位，等同于 i 或 j
x=ich1+qch1.*sqrt(-1);           % 创建复数数组
TrData=ich4+qch4.*sqrt(-1);      % 复数数据
```

---

## 8. 绘图

```matlab
figure;                          % 创建新图形窗口
plot(SNR,corr,'r--d');           % 绘制折线图
% 参数说明：
%   SNR: x轴数据
%   corr: y轴数据
%   'r--d': 红色虚线，数据点用菱形标记

xlabel('SNR');                   % x轴标签
ylabel('正确率');                 % y轴标签
legend('信号检测正确率');         % 图例
grid off;                         % 关闭网格
```

---

## 9. 特殊语法

### 9.1 冒号运算符

```matlab
1:10                             % 生成 [1, 2, 3, ..., 10]
-10:1:0                          % 从-10到0，步长为1
m-1:-1:2                         % 从m-1递减到2，步长为-1
```

### 9.2 转置运算符

```matlab
(linspace(0,1,length(Pxx_int)))'  % `'` 表示转置（共轭转置）
```

### 9.3 空数组

```matlab
ParaBitSig=[];                   % 创建空数组
```

---

## 10. 代码组织

### 10.1 代码节（Section）

```matlab
%% 信号1                          % 代码节标记
% 相关代码...

%% 两个信号叠加                   % 另一个代码节
% 相关代码...
```

代码节可以用于：
- 代码折叠和展开
- 快速导航
- 分段执行（按 `Ctrl+Enter`）

### 10.2 多行语句

MATLAB 中，一行可以写多个语句（用分号或逗号分隔），但通常不推荐：

```matlab
% 不推荐的做法
a=1; b=2; c=3;
```

---

## 11. 常见编程模式

### 11.1 预分配数组

```matlab
peaks = linspace(0,0,m);         % 预分配数组，提高性能
valleys = linspace(0,0,m);
```

### 11.2 计数器模式

```matlab
peakindex = 1;                   % 初始化计数器
for i = 2:m-1
    if 条件
        peaks(peakindex)=row_acc(i);
        peakindex = peakindex+1;  % 递增计数器
    end
end
```

### 11.3 累积计算

```matlab
for x=1:1:size(Pxx,1)
    if x==1
        Pxx_int(x)=(1/(size(Pxx,1)-1))*Pxx(x+1);
    else
        Pxx_int(x)=Pxx_int(x-1)+(1/(size(Pxx,1)-1))*Pxx(x);  % 累积
    end
end
```

---

## 12. 重要注意事项

1. **索引从1开始**：MATLAB 数组索引从1开始，不是0
2. **区分大小写**：变量名 `SNR` 和 `snr` 是不同的
3. **矩阵运算 vs 元素级运算**：
   - `*`：矩阵乘法
   - `.*`：元素级乘法
   - `/`：矩阵除法
   - `./`：元素级除法
4. **分号的作用**：
   - 语句末尾加分号：不显示结果
   - 不加分号：显示结果
5. **数组是列优先**：MATLAB 按列存储数组
6. **动态类型**：变量类型由赋值决定，无需声明

---

## 13. 调试技巧

1. **使用分号**：在循环和大量计算中，使用分号避免输出过多信息
2. **注释代码**：使用 `%` 注释掉不需要的代码，而不是删除
3. **使用断点**：在代码行左侧点击设置断点，便于调试
4. **检查变量**：在命令窗口输入变量名查看其值
5. **使用 `disp()`**：在代码中插入 `disp(变量名)` 输出中间结果

---

## 14. 示例：理解代码片段

让我们分析文件中的一段代码：

```matlab
for i=4:length(y_)-4
    n1(i)=3*(loc(i,2))-(loc(i-1,2))-(loc(i-2,2))-(loc(i-3,2));
    n2(i)=loc(i+3,2)+loc(i+2,2)+loc(i+1,2)-3*(loc(i,2));
    d(i)=n2(i)-n1(i);
end
```

**解释**：
- `for i=4:length(y_)-4`：循环从4开始，到 `length(y_)-4` 结束
- `loc(i,2)`：访问矩阵 `loc` 的第 `i` 行第2列
- `n1(i)`、`n2(i)`、`d(i)`：计算并存储结果
- 这是拐点检测算法的一部分

---

## 总结

MATLAB 语法特点：
- ✅ 简洁直观，适合数值计算
- ✅ 强大的矩阵运算能力
- ✅ 丰富的内置函数库
- ✅ 良好的可视化功能
- ⚠️ 注意索引从1开始
- ⚠️ 注意矩阵运算和元素级运算的区别

通过阅读和理解 `indentify_signal_exist.m` 文件，可以掌握 MATLAB 的基本语法和编程模式。

