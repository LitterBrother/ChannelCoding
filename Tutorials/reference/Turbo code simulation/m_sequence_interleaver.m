function Beta=m_sequence_interleaver(init_x,feedback) 
%����m���н�֯�������ĺ�����ͨ��m���е����ɻ��ƣ�����״̬�Ĵ���������״̬
%״̬��һ��10λ�����Ƶ����У�����ת��Ϊ��Ӧ��10���ƣ���Ϊ�������
%����init_xΪ��ʼ״̬��feedbackΪ����ϵ��������0,1���б�ʾ
for i=1:2^(length(init_x))-1
    s=0;
    Beta(i)=bin_to_dec(init_x);  %��0,1��������ʾ��״̬ת��Ϊʮ������������Beta��
    for k=1:10
        s=s+init_x(10-k+1)*feedback(k);  %���з����ߵļĴ�����ֵ���
    end
    a=rem(s,2);  %ģ����
    init_x=trans(init_x,a);  %�����ͻ���λ
end

function output=trans(x,a)  %��λ��������a��������x����λ����x����һλ
for i=1:length(x)-1
    x(length(x)-i+1)=x(length(x)-i);
end
    x(1)=a;
output=x;  


function m=bin_to_dec(x)  %����������ת��Ϊ��Ӧ��ʮ����
for i=1:length(x)
    s(i)=x(i)*2^(length(x)-i);
end
m=sum(s); 
    