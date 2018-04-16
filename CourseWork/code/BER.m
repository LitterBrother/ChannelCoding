close all;clc;
[A]=xlsread('(7,5)_HardViterbi.csv');
[B]=xlsread('(7,5)_SoftViterbi.csv');
[C]=xlsread('(7,5)_uncode.csv');
[D]=xlsread('(7,5)_bcjr.csv');
hard_viterbi=semilogy(A(:,1),A(:,2),'-kd','MarkerSize',10);hold on;
soft_viterbi=semilogy(B(:,1),B(:,2),'-k^','MarkerSize',10);hold on;
uncode=semilogy(C(:,1),C(:,2),'-ks','MarkerSize',10);hold on;
bcjr=semilogy(D(:,1),D(:,2),'-ko','MarkerSize',10);
set(gca,'FontSize',15,'XTick',0:1:12);
grid on;
xlabel('E_b/N_0(dB)');
ylabel('BER');
title('BER performance of (7,5)_8 conv. code over AWGN channel using BPSK');
legend([uncode,hard_viterbi,soft_viterbi,bcjr],'uncoded','Hard decision Viterbi','Soft decision Viterbi','BCJR');
axis([0 12 1e-6 1]);
% Constraint Length
% figure;
% hard_viterbi=semilogy(A(:,1),A(:,2),'-kd','MarkerSize',10);hold on;
% hard_viterbi_15_13=semilogy(E(:,1),E(:,2),'-ks','MarkerSize',10);hold on;
% xlabel('E_b/N_0(dB)');
% ylabel('BER');
% title('BER performance of rate 1/2 conv. code over AWGN channel using BPSK');
% legend([hard_viterbi,hard_viterbi_15_13],'K=3','K=4');
% axis([0 12 1e-6 1]);

% BCJR Decoding Algorithm
% different constraint length
figure;
[A]=xlsread('(7,5)_bcjr.csv');
[B]=xlsread('(15,13)_bcjr.csv');
[C]=xlsread('(23,35)_bcjr.csv');
bcjr_4state=semilogy(A(:,1),A(:,2),'-ks','MarkerSize',10);hold on;
bcjr_8state=semilogy(B(:,1),B(:,2),'-kd','MarkerSize',10);hold on;
bcjr_16state=semilogy(C(:,1),C(:,2),'-k^','MarkerSize',10);
set(gca,'FontSize',15,'XTick',0:1:12);
grid on;
xlabel('E_b/N_0(dB)');
ylabel('BER');
title('BER performance of different conv. code over AWGN channel using BPSK');
legend([bcjr_4state,bcjr_8state,bcjr_16state],'4-state (7,5)_8 conv. code','8-state (15,13)_8 conv. code','16-state (23,35)_8 conv. code');
axis([0 10 1e-6 1]);
% BCJR Decoding Algorithm
% different rate
figure;
[A]=xlsread('(7,5)_bcjr.csv');
[B]=xlsread('(7,7,5)_bcjr.csv');
bcjr_rate_1_2=semilogy(A(:,1),A(:,2),'-ks','MarkerSize',10);hold on;
bcjr_rate_1_3=semilogy(B(:,1),B(:,2),'-kd','MarkerSize',10);hold on;
set(gca,'FontSize',15,'XTick',0:1:12);
grid on;
xlabel('E_b/N_0(dB)');
ylabel('BER');
title('BER performance of different rate conv. code over AWGN channel using BPSK');
legend([bcjr_rate_1_2,bcjr_rate_1_3],'rate 1/2 (7,5)_8 conv. code','rate 1/3 (7,7,5)_8 conv. code');
axis([0 10 1e-6 1]);
