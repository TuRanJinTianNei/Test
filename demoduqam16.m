function [demodu_bit_symbol]=demoduqam16(Rx_serial_complex_symbols)

% 功能: 16QAM星座的最小欧氏距离判决解调，输出比特流
% 输入: Rx_serial_complex_symbols —— 串行复数符号序列（IQ平面散点）
% 输出: demodu_bit_symbol —— 对应的比特序列（按左高位输出）

% 将得到的串行16QAM数据星座逆映射为比特流

complex_symbols=reshape(Rx_serial_complex_symbols,length(Rx_serial_complex_symbols),1);
d=1;
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
      for i=1:length(Rx_serial_complex_symbols)
          for j=1:16
              metrics(j)=abs(complex_symbols(i,1)-complex_mapping(j,1));
          end
          % 将离某星座点最近的值赋给decode_symble(i)
          [min_metric decode_symbol(i)]=min(metrics); 
      end
      decode_bit_symbol=de2bi((decode_symbol-1)','left-msb');
      demodu_bit_symbol=reshape(decode_bit_symbol',1,length(Rx_serial_complex_symbols)*4);
      
