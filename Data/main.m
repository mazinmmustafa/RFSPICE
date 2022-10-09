close all; clear; clc;
%%
GHz     =   1E9;
Data    =   load('Test.dat');
%%
j       =   sqrt(-1);
freq    =   Data(:,1);
S11     =   Data(:,2)+j*Data(:,3);
S12     =   Data(:,4)+j*Data(:,5);
S21     =   Data(:,6)+j*Data(:,7);
S22     =   Data(:,8)+j*Data(:,9);
%%
figure()
hold on
plot(freq/GHz,20*log10(abs(S11)),'-','LineWidth',1)
plot(freq/GHz,20*log10(abs(S12)),'-','LineWidth',1)
plot(freq/GHz,20*log10(abs(S21)),'--','LineWidth',1)
plot(freq/GHz,20*log10(abs(S22)),'--','LineWidth',1)
hold off
xlabel('Freq [GHz]','Interpret','Latex','FontSize',12)
ylabel('$|S|$ [dB]','Interpret','Latex','FontSize',12)
set(gca,'TickLabel','Latex','FontSize',12)
legend({'$S_{11}$','$S_{12}$','$S_{21}$','$S_{22}$'},...
    'Interpreter','Latex','FontSize',12,'NumColumns',2,...
    'Location','SouthEast')
xlim([1 2])
ylim([-30 0])
pbaspect([2 1 1])
%%