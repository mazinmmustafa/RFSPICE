close all; clear; clc;
%% Definitions
GHz     =   1E9;
%%
Data    =   csvread('ResultsData.csv');
freq    =   Data(:,1);
[Ns,M]  =   size(Data);
S       =   zeros(Ns, (M-1)/2);
for i=1:Ns
    count = 0;
    for ii=2:2:M-1
        count       =   count+1;
        S(i,count) 	=   Data(i,ii)+1i*Data(i,ii+1);
    end
end
%% Reflection S11
figure()
hold on
plot(freq'/GHz,20*log10(abs(S(:,1))),'-','LineWidth',1)
hold off
xlabel('Freq [GHz]','Interpret','Latex','FontSize',12)
ylabel('$|S_{11}|$ [dB]','Interpret','Latex','FontSize',12)
set(gca,'TickLabel','Latex','FontSize',12)
pbaspect([2 1 1])
xlim([1 2])
ylim([-30 0])
%%
exportgraphics(gcf,'Figure1.pdf','ContentType','vector');
%% Transmission S21
figure()
hold on
plot(freq'/GHz,20*log10(abs(S(:,3))),'-','LineWidth',1)
hold off
xlabel('Freq [GHz]','Interpret','Latex','FontSize',12)
ylabel('$|S_{21}|$ [dB]','Interpret','Latex','FontSize',12)
set(gca,'TickLabel','Latex','FontSize',12)
pbaspect([2 1 1])
xlim([1 2])
ylim([-30 0])
%%
exportgraphics(gcf,'Figure2.pdf','ContentType','vector');
%%