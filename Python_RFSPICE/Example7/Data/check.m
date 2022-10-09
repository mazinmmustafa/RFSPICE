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
%%
%% Reflection S11
figure()
hold on
plot(freq'/GHz,20*log10(abs(S(:,1))),'-','LineWidth',1)
hold off
xlabel('Freq [GHz]','Interpret','Latex','FontSize',12)
ylabel('$|S_{11}|$ [dB]','Interpret','Latex','FontSize',12)
set(gca,'TickLabel','Latex','FontSize',12)
pbaspect([2 1 1])
xlim([0 2])
ylim([-60 0])
%%
exportgraphics(gcf,'Figure1.pdf','ContentType','vector');
%% Reflection S22
figure()
hold on
plot(freq'/GHz,20*log10(abs(S(:,6))),'-','LineWidth',1)
hold off
xlabel('Freq [GHz]','Interpret','Latex','FontSize',12)
ylabel('$|S_{22}|$ [dB]','Interpret','Latex','FontSize',12)
set(gca,'TickLabel','Latex','FontSize',12)
pbaspect([2 1 1])
xlim([0 2])
ylim([-60 0])
%%
exportgraphics(gcf,'Figure2.pdf','ContentType','vector');
%% Reflection S33
figure()
hold on
plot(freq'/GHz,20*log10(abs(S(:,11))),'-','LineWidth',1)
hold off
xlabel('Freq [GHz]','Interpret','Latex','FontSize',12)
ylabel('$|S_{33}|$ [dB]','Interpret','Latex','FontSize',12)
set(gca,'TickLabel','Latex','FontSize',12)
pbaspect([2 1 1])
xlim([0 2])
ylim([-60 0])
%%
exportgraphics(gcf,'Figure3.pdf','ContentType','vector');
%% Reflection S44
figure()
hold on
plot(freq'/GHz,20*log10(abs(S(:,16))),'-','LineWidth',1)
hold off
xlabel('Freq [GHz]','Interpret','Latex','FontSize',12)
ylabel('$|S_{44}|$ [dB]','Interpret','Latex','FontSize',12)
set(gca,'TickLabel','Latex','FontSize',12)
pbaspect([2 1 1])
xlim([0 2])
ylim([-60 0])
%%
exportgraphics(gcf,'Figure4.pdf','ContentType','vector');
%% Next 21
figure()
hold on
plot(freq'/GHz,20*log10(abs(S(:,2))),'-','LineWidth',1)
hold off
xlabel('Freq [GHz]','Interpret','Latex','FontSize',12)
ylabel('$|S_{21}|$ [dB]','Interpret','Latex','FontSize',12)
set(gca,'TickLabel','Latex','FontSize',12)
pbaspect([2 1 1])
xlim([0 2])
ylim([-60 0])
%%
exportgraphics(gcf,'Figure5.pdf','ContentType','vector');
%% Next S43
figure()
hold on
plot(freq'/GHz,20*log10(abs(S(:,12))),'-','LineWidth',1)
hold off
xlabel('Freq [GHz]','Interpret','Latex','FontSize',12)
ylabel('$|S_{43}|$ [dB]','Interpret','Latex','FontSize',12)
set(gca,'TickLabel','Latex','FontSize',12)
pbaspect([2 1 1])
xlim([0 2])
ylim([-60 0])
%%
exportgraphics(gcf,'Figure6.pdf','ContentType','vector');
%% Insertion S31
figure()
hold on
plot(freq'/GHz,20*log10(abs(S(:,3))),'-','LineWidth',1)
hold off
xlabel('Freq [GHz]','Interpret','Latex','FontSize',12)
ylabel('$|S_{31}|$ [dB]','Interpret','Latex','FontSize',12)
set(gca,'TickLabel','Latex','FontSize',12)
pbaspect([2 1 1])
xlim([0 2])
ylim([-60 0])
%%
exportgraphics(gcf,'Figure7.pdf','ContentType','vector');
%% Insertion 42
figure()
hold on
plot(freq'/GHz,20*log10(abs(S(:,8))),'-','LineWidth',1)
hold off
xlabel('Freq [GHz]','Interpret','Latex','FontSize',12)
ylabel('$|S_{42}|$ [dB]','Interpret','Latex','FontSize',12)
set(gca,'TickLabel','Latex','FontSize',12)
pbaspect([2 1 1])
xlim([0 2])
ylim([-60 0])
%%
exportgraphics(gcf,'Figure8.pdf','ContentType','vector');
%% Fext S41
figure()
hold on
plot(freq'/GHz,20*log10(abs(S(:,4))),'-','LineWidth',1)
hold off
xlabel('Freq [GHz]','Interpret','Latex','FontSize',12)
ylabel('$|S_{41}|$ [dB]','Interpret','Latex','FontSize',12)
set(gca,'TickLabel','Latex','FontSize',12)
pbaspect([2 1 1])
xlim([0 2])
ylim([-60 0])
%%
exportgraphics(gcf,'Figure9.pdf','ContentType','vector');
%% Fext S32
figure()
hold on
plot(freq'/GHz,20*log10(abs(S(:,7))),'-','LineWidth',1)
hold off
xlabel('Freq [GHz]','Interpret','Latex','FontSize',12)
ylabel('$|S_{32}|$ [dB]','Interpret','Latex','FontSize',12)
set(gca,'TickLabel','Latex','FontSize',12)
pbaspect([2 1 1])
xlim([0 2])
ylim([-60 0])
%%
exportgraphics(gcf,'Figure10.pdf','ContentType','vector');
%%

