function x_decode_msg=turbo_decode(rec,g,Alpha,puncture,niter,method)
%Turbo码译码函数，rec为接收序列，即+1，-1序列加入噪声后得到的序列，niter为迭代次数，输出x_decode_msg为0,1行向量
[n,k]=size(g); 
m=k-1;
L=length(rec);

yk=demultiplex(rec,Alpha,puncture);
L_e(1:1:L/(2+puncture))=zeros(1,L/(2+puncture));
for iter=1:niter  %解码器1
    L_a(Alpha)=L_e;
    if method==0
        L_all=logmapo(yk(1,:),g,L_a,1);
    else
        L_all=sova0(yk(1,:),g,L_a,1);
    end
    L_e=L_all-2*yk(1,1:2:2*L/(2+puncture))-L_a;
    L_a=L_e(Alpha);  %解码器2
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