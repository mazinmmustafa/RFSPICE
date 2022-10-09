close all; clear; clc;
%% Definitions
GHz     =   1E9;
%%
Data    =   load('data.dat');
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
VS  	=   zeros(Ns,1);
for i=1:Ns
    SM      =   [ S(i,1) S(i,2) ; S(i,3) S(i,4) ]; 
    Z0      =   50*eye(2);
    Z       =   (eye(2)+SM)*inv(eye(2)-SM)*Z0;
    Zp      =   [ 50 0 ; 0 100 ];
    Vp      =   [ 1.0 ; 0 ];
    V_      =   inv(eye(2)+Zp*inv(Z))*Vp;
    VS(i)   =   V_(2);
end
%% 
figure()
hold on
plot(freq',abs(VS),'-','LineWidth',1)
hold off
xlabel('Freq [Hz]','Interpret','Latex','FontSize',12)
ylabel('$|V|$ [V]','Interpret','Latex','FontSize',12)
set(gca,'TickLabel','Latex','FontSize',12)
set(gca,'XScale','log','YScale','log')
pbaspect([2 1 1])
xlim([1E6 0.5E9])
ylim([1E-3 1E0])
%%
exportgraphics(gcf,'Figure1.pdf','ContentType','vector');
%%
figure()
hold on
plot(freq',(180/pi)*angle(VS),'-','LineWidth',1)
hold off
xlabel('Freq [Hz]','Interpret','Latex','FontSize',12)
ylabel('$\angle V$ [deg]','Interpret','Latex','FontSize',12)
set(gca,'TickLabel','Latex','FontSize',12)
set(gca,'XScale','log')
pbaspect([2 1 1])
xlim([1E6 0.5E9])
ylim([-1 +1]*200)
%%
exportgraphics(gcf,'Figure2.pdf','ContentType','vector');
%%