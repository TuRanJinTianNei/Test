function [complex_qam_data]=qam16(bitdata)

% the function is aim to modulate bitdata to 16QAM complex signal

x1=reshape(bitdata,4,length(bitdata)/4)';
d=1;  %min distance of symble 
for i=1:length(bitdata)/4
    for j=1:4
        x1(i,j)=x1(i,j)*(2^(4-j)); % 16进制数对位转化成10进制
    end
    source(i,1)=1+sum(x1(i,:)); % 将16进制数转化为十进制后，再加一对位
end

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
	         d  -d;
           3*d  -d;
	    -3*d  -3*d;
	      -d  -3*d;
	       d  -3*d;
	    3*d  -3*d];
for i=1:length(bitdata)/4
    qam_data(i,:)=mapping(source(i),:); % data mapping
end
complex_qam_data=complex(qam_data(:,1),qam_data(:,2));

    

