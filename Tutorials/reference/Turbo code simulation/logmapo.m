function L_all=logmapo(rec_s,g,L_a,ind_dec)
%Log-Map译码算法
L_total=length(rec_s)/2;
[n,K]=size(g); 
m=K-1;
nstates=2^m;
[next_out,next_state,last_out,last_state]=trellis(g);
Infty=1e10;

Alpha(1,1)=0; 
Alpha(1,2:nstates)=-Infty*ones(1,nstates-1);  %初始化Alpha

if ind_dec==1
   Beta(L_total,1)=0;
   Beta(L_total,2:nstates)=-Infty*ones(1,nstates-1); 
elseif ind_dec==2
   Beta(L_total,1:nstates)=zeros(1,nstates);
else
   fprintf('ind_dec is limited to 1 and 2!\n');
end  %初始化Beta

for k=2:L_total+1
    for state2=1:nstates
      gamma=-Infty*ones(1,nstates);
      gamma(last_state(state2,1))=(-rec_s(2*k-3)+rec_s(2*k-2)*last_out(state2,2))....
           -log(1+exp(L_a(k-1)));
      gamma(last_state(state2,2))=(rec_s(2*k-3)+rec_s(2*k-2)*last_out(state2,4))....
           +L_a(k-1)-log(1+exp(L_a(k-1)));

      if(sum(exp(gamma+Alpha(k-1,:)))<1e-300)
         Alpha(k,state2)=-Infty;
      else
         Alpha(k,state2)=log(sum(exp(gamma+Alpha(k-1,:))));  
      end   
    end
    tempmax(k)=max(Alpha(k,:));
    Alpha(k,:)=Alpha(k,:)-tempmax(k);
end   

for k=L_total-1:-1:1
  for state1=1:nstates
     gamma=-Infty*ones(1,nstates);
     gamma(next_state(state1,1))=(-rec_s(2*k+1)+rec_s(2*k+2)*next_out(state1,2))....
           -log(1+exp(L_a(k+1)));
     gamma(next_state(state1,2))=(rec_s(2*k+1)+rec_s(2*k+2)*next_out(state1,4))....
           +L_a(k+1)-log(1+exp(L_a(k+1)));
     if(sum(exp(gamma+Beta(k+1,:)))<1e-300)
        Beta(k,state1)=-Infty;
     else
        Beta(k,state1)=log(sum(exp(gamma+Beta(k+1,:))));
     end   
  end
  Beta(k,:)=Beta(k,:)-tempmax(k+1);
end

for k=1:L_total
  for state2=1:nstates
     gamma0=(-rec_s(2*k-1)+rec_s(2*k)*last_out(state2,2))....
           -log(1+exp(L_a(k)));
     gamma1=(rec_s(2*k-1)+rec_s(2*k)*last_out(state2,4))...
           +L_a(k)-log(1+exp(L_a(k)));
     temp0(state2)=exp(gamma0+Alpha(k,last_state(state2,1))+Beta(k,state2));
     temp1(state2)=exp(gamma1+Alpha(k,last_state(state2,2))+Beta(k,state2));
  end
  L_all(k)=log(sum(temp1))-log(sum(temp0));
end