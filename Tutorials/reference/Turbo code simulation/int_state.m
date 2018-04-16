function int_state=int_state(state)
[dummy,m]=size(state);
for i=1:m
   vect(i)=2^(m-i);
end
int_state=state*vect';