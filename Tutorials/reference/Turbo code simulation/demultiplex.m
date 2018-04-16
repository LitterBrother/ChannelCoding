function subr=demultiplex(r,alpha,puncture);  %解复用函数
L_total=length(r)/(2+puncture);
if puncture==1  %不删余       
  for i=1:L_total  
      x_sys(i)=r(3*(i-1)+1);
      for j=1:2
          subr(j,2*i)=r(3*(i-1)+1+j);
      end
   end
else  %删余
   for i=1:L_total
       x_sys(i)=r(2*(i-1)+1);
       for j=1:2
          subr(j,2*i)=0;
       end   
       if rem(i,2)>0
          subr(1,2*i)=r(2*i);
       else
          subr(2,2*i)=r(2*i);
       end      
   end
end       

for j=1:L_total
   subr(1,2*(j-1)+1)=x_sys(j);
   subr(2,2*(j-1)+1)=x_sys(alpha(j));
end    
