%Turbo码仿真的主程序
clear all
niter=3;  %迭代次数为5次
L_frame=1021;  %帧长
n_frame=10;  %帧数
start=0;
step=0.5;
finish=3;  %EbN0从0dB到3dB，步长为0.5
g=[1,1,1;1,0,1];  %生成矩阵
[n,k]=size(g);
m=k-1;
method=0;  %译码方式，0为logmap译码，1为sova译码
puncture=0;  %是否删余，0删余，码率1/2；1不删余，码率1/3
r=1/(puncture+2);  %码率
init_x2=[0 0 0 1 1 1 0 1 0 0];  %优选m序列发生器的初始状态
feedback2=[1 1 1 0 1 1 0 0 0 1 1];  %优选m序列发生器反馈系数
Alpha=berrou_interleaver;
Beta=m_sequence_interleaver(init_x2,feedback2);  %优选m序列交织映射索引
Eb=1;
plot_pe1=[];  %采用Berrou交织器的Turbo码编码时的误码率
plot_pe2=[];  %采用优选m序列交织器的Turbo码编码时的误码率
plot_pe3=[];  %未编码时的误码率
Q=1;  %工作次数
axis_EbN0=start:step:finish;  %横坐标为信噪比

for EbN0=start:step:finish
    Liner_EbN0=10^(EbN0/10);
    pe_number1=0;
    pe_number2=0;  
    pe_number3=0;  %初始化误码个数为0
    variance=0.5*(Eb/Liner_EbN0)/r;  %噪声方差
    for i=1:1:n_frame
        x_msg=randint(1,L_frame,2);  %由randint随机产生0,1输入序列
        x_code_msg1=turbo_encode(x_msg,g,Alpha,puncture);
        x_code_msg2=turbo_encode(x_msg,g,Beta,puncture);  %Turbo码编码
        x_bpsk_msg1=2*x_code_msg1-1; 
        x_bpsk_msg2=2*x_code_msg2-1;  %bpsk调制
        rec1=x_bpsk_msg1+sqrt(variance)*randn(1,length(x_code_msg1));  
        rec2=x_bpsk_msg2+sqrt(variance)*randn(1,length(x_code_msg2));   %通过AWGN信道传输
        x_decode_msg1=turbo_decode(rec1,g,Alpha,puncture,niter,method);  
        x_decode_msg2=turbo_decode(rec2,g,Beta,puncture,niter,method);  %Turbo码解码
        x_uncode_msg=uncode(x_msg,EbN0);  %未经编码时的硬判决输出
        pe_number1=pe_number1+sum(x_msg~=x_decode_msg1);
        pe_number2=pe_number2+sum(x_msg~=x_decode_msg2);
        pe_number3=pe_number3+sum(x_msg~=x_uncode_msg);  %分别统计误码数
        
        current_time=fix(clock);
        fprintf('i am working %g      %g年  %g月  %g日  %g时  %g分  %g秒\n\n',Q,current_time(1),current_time(2),current_time(3),current_time(4),current_time(5),current_time(6))
        Q=Q+1;
        fprintf('\n\n')  %循环标记
        
    end
    pe1=pe_number1/(L_frame*n_frame);  
    pe2=pe_number2/(L_frame*n_frame); 
    pe3=pe_number3/(L_frame*n_frame);  %计算误码率
    plot_pe1=[plot_pe1,pe1];
    plot_pe2=[plot_pe2,pe2];
    plot_pe3=[plot_pe3,pe3];
end

semilogy(axis_EbN0,plot_pe1,'b*-',axis_EbN0,plot_pe2,'bs-',axis_EbN0,plot_pe3,'bo-')  %绘图
legend('分组交织器','优选m序列交织器','无编码输出');
xlabel('Eb/N0 信号功率/噪声功率(db)')
ylabel('BER')
title('误码率')
grid on
        