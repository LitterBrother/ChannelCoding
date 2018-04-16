function Beta=m_sequence_interleaver(init_x,feedback) 
%生成m序列交织器索引的函数，通过m序列的生成机制，遍历状态寄存器的所有状态
%状态是一个10位二进制的序列，将其转化为相应的10进制，作为索引输出
%输入init_x为初始状态，feedback为反馈系数，均用0,1序列表示
for i=1:2^(length(init_x))-1
    s=0;
    Beta(i)=bin_to_dec(init_x);  %将0,1序列所表示的状态转化为十进制整数存入Beta中
    for k=1:10
        s=s+init_x(10-k+1)*feedback(k);  %将有反馈线的寄存器数值相加
    end
    a=rem(s,2);  %模二加
    init_x=trans(init_x,a);  %反馈送回首位
end

function output=trans(x,a)  %移位函数，将a送入序列x的首位并将x右移一位
for i=1:length(x)-1
    x(length(x)-i+1)=x(length(x)-i);
end
    x(1)=a;
output=x;  


function m=bin_to_dec(x)  %二进制序列转化为相应的十进制
for i=1:length(x)
    s(i)=x(i)*2^(length(x)-i);
end
m=sum(s); 
    