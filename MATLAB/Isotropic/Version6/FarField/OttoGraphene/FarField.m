close all; clear; clc;
%%
Ns          =   1e3;
[~,~,~,um,~,~,~,~,~]=Units();
%% Definitions
J           =   1;  
theta0      =   90;
phi0        =   0;
x_          =   0*um;
y_          =   0*um;
z_          =   5*um;
phi         =   0;
Data_min  	=   -25;
%%
theta       =   linspace(0,2*pi,Ns);
[E_theta,E_phi]=FarFieldThetaCut(theta,phi,theta0,phi0,x_,y_,z_,J);
%% Plot
figure()
theta       =   linspace(0,2*pi,Ns);
Data        =   abs(E_theta)./max(abs(E_theta)); 
Data        =   20*log10(Data);
Data(Data<Data_min)=Data_min;
PolardB(theta,Data,[Data_min 0],6,'k',1)
title('$20\log_{10}|E_{\theta}|$ vs $\theta$','Interpret','Latex')
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\FarFieldTheta' -pdf -transparent
%% Plot
figure()
theta       =   linspace(0,2*pi,Ns);
Data        =   abs(E_phi)./max(abs(E_phi)); 
Data        =   20*log10(Data);
Data(Data<Data_min)=Data_min;
PolardB(theta,Data,[Data_min 0],6,'k',1)
title('$20\log_{10}|E_{\phi}|$ vs $\theta$','Interpret','Latex')
%% Export to pdf
% export_fig 'C:\Users\Mazin\Desktop\FarFieldPhi' -pdf -transparent
%%



















%%