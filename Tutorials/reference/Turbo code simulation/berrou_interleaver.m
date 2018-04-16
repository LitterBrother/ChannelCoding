function Alpha=berrou_interleaver
%生成Berrou交织器映射序列的函数
x=1:1024;
m=1;
q=1;
in=reshape(x,32,32)';  %按行写入
for k=32:-1:1
    for n=32:-1:1
        y(m)=in(n,k);  %按列读出
        m=m+1;
    end
end
for p=2:1024
    Alpha(q)=y(p);  %生成交织映射
    q=q+1;
end


    