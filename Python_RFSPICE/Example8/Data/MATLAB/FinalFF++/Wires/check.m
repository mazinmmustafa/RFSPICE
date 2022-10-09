close all; clear; clc;
%%
GHz         =   1E9;
mm          =   1E-3;
%%
c0          =   299792458;
mu0         =   4*pi*1E-7;
eps0        =   1/(mu0*c0^2);
eta0        =   sqrt(mu0/eps0);
%%
sigma       =   58E6;
r           =   1*mm;
S           =   10*mm;
freq        =   1*GHz;
%%
Z0          =   (eta0/pi)*acosh(S/(2*r));
Rs          =   sqrt(pi*freq*mu0/sigma);
R_          =   (Rs/(pi*r))*(S/r)/sqrt((S/r)^2-4);
C_          =   1/(c0*Z0);
L_          =   Z0/c0;
%%
fprintf('R = %0.2f [Ohm/m]\n',R_);
fprintf('L = %0.2f [nH/m]\n',L_/1E-9);
fprintf('C = %0.2f [pF/m]\n',C_/1E-12);
%%