close all;clc;
%% (7,5) convolutional code
% load('7_6_HardViterbi.mat','Eb_N0_dB','simBer_HardViterbi');
% load('7_6_SoftViterbi.mat','Eb_N0_dB','simBer_SoftViterbi');
A = load('(7,5)_HardViterbi.csv');
B = load('(7,5)_SoftViterbi.csv');
C = load('(7,5)_uncode.csv');
D = load('(7,5)_bcjr.csv');
hard_viterbi=semilogy(A(:,1),A(:,2),'-bd','MarkerSize',12,'LineWidth',2);hold on;
soft_viterbi=semilogy(B(:,1),B(:,2),'-r^','MarkerSize',12,'LineWidth',2);hold on;
% hard_viterbi=semilogy(Eb_N0_dB,simBer_HardViterbi,'-kd','MarkerSize',12);hold on;
% soft_viterbi=semilogy(Eb_N0_dB,simBer_SoftViterbi,'-k^','MarkerSize',12);hold on;
uncode=semilogy(C(:,1),C(:,2),'-ks','MarkerSize',12,'LineWidth',2);hold on;
bcjr=semilogy(D(:,1),D(:,2),'-mo','MarkerSize',12,'LineWidth',2);
grid on;
xlabel('$E_b/N_0$(dB)','interpreter','latex');
ylabel('BER');
set(gca,'FontName','Times New Roman','FontSize',12,'xtick',0:1:12,'ytick',[1e-5 1e-4 1e-3 1e-2 1e-1 1]);
title('BER performances');
legend([uncode,hard_viterbi,soft_viterbi,bcjr],'uncoded','Hard decision Viterbi','Soft decision Viterbi','BCJR');
axis([0 10 1e-5 1]);
FigTool(1);

%% BCJR Decoding Algorithm
% different constraint length
figure;
[A] = load('(7,5)_bcjr.csv');
[B] = load('(15,13)_bcjr.csv');
[C] = load('(23,35)_bcjr.csv');
bcjr_4state = semilogy(A(:,1),A(:,2),'-rs','MarkerSize',12,'LineWidth',2);hold on;
bcjr_8state = semilogy(B(:,1),B(:,2),'-bd','MarkerSize',12,'LineWidth',2);hold on;
bcjr_16state = semilogy(C(:,1),C(:,2),'-g^','MarkerSize',12,'LineWidth',2);
set(gca,'FontSize',15,'XTick',0:1:12);
grid on;
xlabel('$E_b/N_0$(dB)','interpreter','latex');
ylabel('BER');
title('BER performances of Soft Decision Viterbi');
legend([bcjr_4state,bcjr_8state,bcjr_16state],'4-state (7,5)_8 conv. code','8-state (15,13)_8 conv. code','16-state (23,35)_8 conv. code');
axis([0 10 1e-5 1]);
set(gca,'FontName','Times New Roman','FontSize',15,'XTick',0:1:12,'ytick',[1e-5,1e-4,1e-3,1e-2,1e-1,1]);
FigTool(1);
%% BCJR Decoding Algorithm
% different polynomial
figure;
[A] = load('(7,5)_SoftViterbi.csv');
[B] = load('7_6_SoftViterbi.mat','Eb_N0_dB','simBer_SoftViterbi');
poly1 = semilogy(A(:,1),A(:,2),'-rs','MarkerSize',12,'LineWidth',2);hold on;
poly2 = semilogy(Eb_N0_dB,simBer_SoftViterbi,'-bd','MarkerSize',12,'LineWidth',2);hold on;
set(gca,'FontName','Times New Roman','FontSize',15,'XTick',0:1:12);

grid on;
xlabel('$E_b/N_0$(dB)','interpreter','latex');
ylabel('BER');
title('BER performances of Soft Decision Viterbi');
legend([poly1, poly2],'(7,5)_8','(7,6)_8');
axis([0 10 1e-5 1]);
set(gca,'FontSize',15,'XTick',0:1:12,'ytick',[1e-5,1e-4,1e-3,1e-2,1e-1,1]);
FigTool(1);