function Alpha=berrou_interleaver
%����Berrou��֯��ӳ�����еĺ���
x=1:1024;
m=1;
q=1;
in=reshape(x,32,32)';  %����д��
for k=32:-1:1
    for n=32:-1:1
        y(m)=in(n,k);  %���ж���
        m=m+1;
    end
end
for p=2:1024
    Alpha(q)=y(p);  %���ɽ�֯ӳ��
    q=q+1;
end


    