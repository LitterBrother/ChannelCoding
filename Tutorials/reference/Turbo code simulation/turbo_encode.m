function x_code_msg=turbo_encode(x_msg,g,Alpha,puncture)
%Turbo码编码函数，Alpha为交织映射表，输出en_output为0,1行向量
[n,k]=size(g); 
m=k-1;
L_info=length(x_msg); 
L_total=L_info+m;  
input=x_msg;
output1=rsc_encode(g,input,1);  %rsc编码
y(1,:)=output1(1:2:2*L_total);  %输出信息位
y(2,:)=output1(2:2:2*L_total);  %输出校验位
for i=1:L_total
   input1(1,i)=y(1,Alpha(i)); 
end  %对输入信息进行交织
output2=rsc_encode(g,input1(1,1:L_total),-1);  %交织后的rsc编码
y(3,:)=output2(2:2:2*L_total);
if puncture>0  %不删余		
   for i=1:L_total
       for j=1:3
           en_output(1,3*(i-1)+j)=y(j,i);
       end
   end
else  %删余		
   for i=1:L_total
       en_output(1,n*(i-1)+1)=y(1,i);
       if rem(i,2)  %rsc1的奇数位
           en_output(1,n*i)=y(2,i);
       else  %rsc2的偶数位
           en_output(1,n*i) = y(3,i);
       end 
    end  
end
x_code_msg=en_output;