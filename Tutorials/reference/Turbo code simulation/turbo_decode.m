function x_decode_msg=turbo_decode(rec,g,Alpha,puncture,niter,method)
%Turbo�����뺯����recΪ�������У���+1��-1���м���������õ������У�niterΪ�������������x_decode_msgΪ0,1������
[n,k]=size(g); 
m=k-1;
L=length(rec);

yk=demultiplex(rec,Alpha,puncture);
L_e(1:1:L/(2+puncture))=zeros(1,L/(2+puncture));
for iter=1:niter  %������1
    L_a(Alpha)=L_e;
    if method==0
        L_all=logmapo(yk(1,:),g,L_a,1);
    else
        L_all=sova0(yk(1,:),g,L_a,1);
    end
    L_e=L_all-2*yk(1,1:2:2*L/(2+puncture))-L_a;
    L_a=L_e(Alpha);  %������2
    if method==0
        L_all=logmapo(yk(2,:),g,L_a,2);
    else
        L_all=sova0(yk(2,:),g,L_a,2);
    end
    L_e=L_all-2*yk(2,1:2:2*L/(2+puncture))-L_a;
end

x_decode_msg(Alpha)=(sign(L_all)+1)/2;
L=length(x_decode_msg);
tmp=x_decode_msg(1:1:L-m);
x_decode_msg=tmp;