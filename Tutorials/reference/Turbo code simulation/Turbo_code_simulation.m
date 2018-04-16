%Turbo������������
clear all
niter=3;  %��������Ϊ5��
L_frame=1021;  %֡��
n_frame=10;  %֡��
start=0;
step=0.5;
finish=3;  %EbN0��0dB��3dB������Ϊ0.5
g=[1,1,1;1,0,1];  %���ɾ���
[n,k]=size(g);
m=k-1;
method=0;  %���뷽ʽ��0Ϊlogmap���룬1Ϊsova����
puncture=0;  %�Ƿ�ɾ�࣬0ɾ�࣬����1/2��1��ɾ�࣬����1/3
r=1/(puncture+2);  %����
init_x2=[0 0 0 1 1 1 0 1 0 0];  %��ѡm���з������ĳ�ʼ״̬
feedback2=[1 1 1 0 1 1 0 0 0 1 1];  %��ѡm���з���������ϵ��
Alpha=berrou_interleaver;
Beta=m_sequence_interleaver(init_x2,feedback2);  %��ѡm���н�֯ӳ������
Eb=1;
plot_pe1=[];  %����Berrou��֯����Turbo�����ʱ��������
plot_pe2=[];  %������ѡm���н�֯����Turbo�����ʱ��������
plot_pe3=[];  %δ����ʱ��������
Q=1;  %��������
axis_EbN0=start:step:finish;  %������Ϊ�����

for EbN0=start:step:finish
    Liner_EbN0=10^(EbN0/10);
    pe_number1=0;
    pe_number2=0;  
    pe_number3=0;  %��ʼ���������Ϊ0
    variance=0.5*(Eb/Liner_EbN0)/r;  %��������
    for i=1:1:n_frame
        x_msg=randint(1,L_frame,2);  %��randint�������0,1��������
        x_code_msg1=turbo_encode(x_msg,g,Alpha,puncture);
        x_code_msg2=turbo_encode(x_msg,g,Beta,puncture);  %Turbo�����
        x_bpsk_msg1=2*x_code_msg1-1; 
        x_bpsk_msg2=2*x_code_msg2-1;  %bpsk����
        rec1=x_bpsk_msg1+sqrt(variance)*randn(1,length(x_code_msg1));  
        rec2=x_bpsk_msg2+sqrt(variance)*randn(1,length(x_code_msg2));   %ͨ��AWGN�ŵ�����
        x_decode_msg1=turbo_decode(rec1,g,Alpha,puncture,niter,method);  
        x_decode_msg2=turbo_decode(rec2,g,Beta,puncture,niter,method);  %Turbo�����
        x_uncode_msg=uncode(x_msg,EbN0);  %δ������ʱ��Ӳ�о����
        pe_number1=pe_number1+sum(x_msg~=x_decode_msg1);
        pe_number2=pe_number2+sum(x_msg~=x_decode_msg2);
        pe_number3=pe_number3+sum(x_msg~=x_uncode_msg);  %�ֱ�ͳ��������
        
        current_time=fix(clock);
        fprintf('i am working %g      %g��  %g��  %g��  %gʱ  %g��  %g��\n\n',Q,current_time(1),current_time(2),current_time(3),current_time(4),current_time(5),current_time(6))
        Q=Q+1;
        fprintf('\n\n')  %ѭ�����
        
    end
    pe1=pe_number1/(L_frame*n_frame);  
    pe2=pe_number2/(L_frame*n_frame); 
    pe3=pe_number3/(L_frame*n_frame);  %����������
    plot_pe1=[plot_pe1,pe1];
    plot_pe2=[plot_pe2,pe2];
    plot_pe3=[plot_pe3,pe3];
end

semilogy(axis_EbN0,plot_pe1,'b*-',axis_EbN0,plot_pe2,'bs-',axis_EbN0,plot_pe3,'bo-')  %��ͼ
legend('���齻֯��','��ѡm���н�֯��','�ޱ������');
xlabel('Eb/N0 �źŹ���/��������(db)')
ylabel('BER')
title('������')
grid on
        