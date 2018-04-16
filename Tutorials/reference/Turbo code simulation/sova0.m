function L_all=sova(rec_s,g,L_a,ind_dec)
%软输出维特比译码算法
L_total=length(L_a);
[n,K]=size(g); 
m=K-1;
nstates=2^m;
Infty=1e10;
delta=30;
[next_out,next_state,last_out,last_state]=trellis(g);
for t=1:L_total+1
   for state=1:nstates
      path_metric(state,t)=-Infty;
   end
end
path_metric(1,1)=0;
for t=1:L_total
   y=rec_s(2*t-1:2*t);
   for state=1:nstates
      sym0=last_out(state,1:2);
      sym1=last_out(state,3:4);
      state0=last_state(state,1);
      state1=last_state(state,2);
      Mk0=y*sym0'-L_a(t)/2+path_metric(state0,t);
      Mk1=y*sym1'+L_a(t)/2+path_metric(state1,t);
      
      if Mk0>Mk1
         path_metric(state,t+1)=Mk0;
         Mdiff(state,t+1)=Mk0-Mk1;
         prev_bit(state, t+1)=0;
      else
         path_metric(state,t+1)=Mk1;
         Mdiff(state,t+1)=Mk1-Mk0;
         prev_bit(state,t+1)=1;
      end
   end
end
if ind_dec==1
   mlstate(L_total+1)=1;
else
   mlstate(L_total+1)=find(path_metric(:,L_total+1)==max(path_metric(:,L_total+1)) );
end

for t=L_total:-1:1
   est(t)=prev_bit(mlstate(t+1),t+1);
   mlstate(t)=last_state(mlstate(t+1),est(t)+1);
end

for t=1:L_total
   llr=Infty;
   for i=0:delta
      if t+i<L_total+1
         bit=1-est(t+i);
         temp_state=last_state(mlstate(t+i+1),bit+1);
         for j=i-1:-1:0
            bit=prev_bit(temp_state,t+j+1);
            temp_state=last_state(temp_state,bit+1);
         end
         if bit~=est(t) 
            llr=min(llr,Mdiff(mlstate(t+i+1),t+i+1));
         end
      end
   end
   L_all(t)=(2*est(t)-1)*llr;
end    