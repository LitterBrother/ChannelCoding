function x_uncode_msg=uncode(x_msg,EbN0)  %δ����ʱ��Ӳ�о����
Eb=1;
Liner_EbN0=10^(EbN0/10);
x_bpsk_msg=2*x_msg-1;
variance=0.5*(Eb/Liner_EbN0);  %�����������δ�����룬����r=1
rec=x_bpsk_msg+sqrt(variance)*randn(1,length(x_msg));
for i=1:length(x_msg)  %Ӳ�о�
    if rec(i)<=0
        x_uncode_msg(i)=0;
    else 
        x_uncode_msg(i)=1;
    end
end
