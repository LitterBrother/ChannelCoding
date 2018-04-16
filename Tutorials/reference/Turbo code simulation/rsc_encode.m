function y=rsc_encode(g,x,terminated)
%rsc编码子程序
%terminated为尾比特处理标志，如terminated>0，有m个尾比特，编码至x最后一个比特到达最后一个寄存器；
%如terminated<0，没有尾比特，编码至x最后一个比特进入寄存器
[n,K]=size(g);
m=K-1;
if terminated>0  %由terminated决定编码输出
  L_info=length(x);
  L_total=L_info+m;
else
  L_total=length(x);
  L_info=L_total-m;
end  
state=zeros(1,m);  %初始化状态向量
for i=1:L_total  %产生码字
   if terminated<0|(terminated>0&i<=L_info)
      d_k=x(1,i);
   elseif terminated>0&i>L_info
      d_k=rem(g(1,2:K)*state',2);  %尾比特处理
   end
 
   a_k=rem(g(1,:)*[d_k state]',2);  %a_k是编码器的第一个寄存器输入
   [output_bits,state]=encode_bit(g,a_k,state);  %进行编码
   output_bits(1,1)=d_k;  %编码比特的第一位是信息位
   y(n*(i-1)+1:n*i)=output_bits;  %编码比特：信息位，校验位1，校验位2，……，校验位n-1，信息位，……
end