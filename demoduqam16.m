% ============================================================================
% 文件名: demoduqam16.m
% 功能: 16QAM解调函数（最小欧氏距离判决）
% 描述:
%   对接收到的16QAM复数符号进行硬判决解调，输出比特流
%   使用最小欧氏距离准则进行符号判决
% 输入参数:
%   Rx_serial_complex_symbols - 接收到的串行复数符号序列（IQ平面散点）
% 输出参数:
%   demodu_bit_symbol - 解调后的比特序列（按左高位输出）
% ============================================================================
function [demodu_bit_symbol]=demoduqam16(Rx_serial_complex_symbols)

%------------------------------------------------------------------------------
% 【算法流程：最小欧氏距离判决】
%------------------------------------------------------------------------------
complex_symbols=reshape(Rx_serial_complex_symbols,length(Rx_serial_complex_symbols),1);
d=1;  % 最小距离（与qam16.m中的d保持一致）

% 16QAM星座映射表（与qam16.m中的映射表一致）
mapping=[-3*d 3*d;
	   -d  3*d;
        d  3*d;
	  3*d  3*d;
       -3*d  d;
	    -d   d;
	     d   d;
	    3*d  d;
 	  -3*d  -d; 
	    -d  -d; 
	    d   -d;
       3*d  -d;
	 -3*d  -3*d;
	   -d  -3*d;
	    d  -3*d;
	  3*d  -3*d];
  
complex_mapping=complex(mapping(:,1),mapping(:,2));

% 对每个接收符号进行最小距离判决
for i=1:length(Rx_serial_complex_symbols)
    % 计算接收符号与所有16个理想星座点的距离
    for j=1:16
        metrics(j)=abs(complex_symbols(i,1)-complex_mapping(j,1));
    end
    % 找到距离最近的星座点索引
    [min_metric decode_symbol(i)]=min(metrics); 
end

% 将星座点索引转换为比特流（左高位MSB格式）
decode_bit_symbol=de2bi((decode_symbol-1)','left-msb');
demodu_bit_symbol=reshape(decode_bit_symbol',1,length(Rx_serial_complex_symbols)*4);
