function [next_out,next_state,last_out,last_state]=trellis(g)
[n,K]=size(g);
m=K-1;
max_state=2^m;
for state=1:max_state
   state_vector=bin_state(state-1,m);
   d_k=0;
   a_k=rem(g(1,:)*[0 state_vector]',2);
   [out_0,state_0]=encode_bit(g,a_k,state_vector);
   out_0(1)=0;
   d_k=1;
   a_k=rem(g(1,:)*[1 state_vector]',2);
   [out_1,state_1]=encode_bit(g,a_k,state_vector);
   out_1(1)=1;
   next_out(state,:)=2*[out_0 out_1]-1;
   next_state(state,:)=[(int_state(state_0)+1) (int_state(state_1)+1)];
end

last_state=zeros(max_state,2);
for bit=0:1
   for state=1:max_state
      last_state(next_state(state,bit+1),bit+1)=state;
      last_out(next_state(state, bit+1),bit*2+1:bit*2+2) ...
         =next_out(state,bit*2+1:bit*2+2);
   end 
end