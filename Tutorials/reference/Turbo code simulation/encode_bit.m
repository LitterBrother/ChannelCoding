function [output,state]=encode_bit(g,input,state)
%得到编码比特和状态转移的子程序
[n,k]=size(g);
m=k-1;
for i=1:n  %得到下一个输出比特
   output(i)=g(i,1)*input;
   for j=2:k
      output(i)=xor(output(i),g(i,j)*state(j-1));
   end;
end
state=[input,state(1:m-1)];