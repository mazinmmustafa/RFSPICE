close all; clear; clc;
%% Definitions
GHz     =   1E9;
MHz     =   1E6;
m       =   1;
%%
Ns      =   1601;
L1    	=   2.0*m;
L2    	=   1.0*m;
L3    	=   1.0*m;
fmin    =   1.0*MHz;
fmax    =   0.5*GHz;
%%
%% Load L & C matricies
L1_  	=   load('L1.dat');
C1_     =   load('C1.dat');
L2_  	=   load('L2.dat');
C2_     =   load('C2.dat');
L3_  	=   load('L3.dat');
C3_     =   load('C3.dat');
%%
freq    =   logspace(log10(fmin),log10(fmax),Ns);
file1   =   fopen('Device1.csv','w');
file2   =   fopen('Device2.csv','w');
file3   =   fopen('Device3.csv','w');
for i=1:Ns
    SIdeal_1	=   CTL(freq(i),L1,0*L1_,0*L1_,L1_,C1_);
    SIdeal_2	=   CTL(freq(i),L2,0*L2_,0*L2_,L2_,C2_);
    SIdeal_3	=   CTL(freq(i),L3,0*L3_,0*L3_,L3_,C3_);
    [M,~]       =   size(SIdeal_1);
    fprintf(file1,"%21.14E,",freq(i));
    for ii=1:M
        for iii=1:M
            fprintf(file1,"%21.14E,%21.14E",...
                real(SIdeal_1(ii,iii)),imag(SIdeal_1(ii,iii)));
            if ~(ii==M && iii==M)
                fprintf(file1,",");
            end
        end
    end
    fprintf(file1,"\n");
    %
    [M,~]       =   size(SIdeal_2);
    fprintf(file2,"%21.14E,",freq(i));
    for ii=1:M
        for iii=1:M
            fprintf(file2,"%21.14E,%21.14E",...
                real(SIdeal_2(ii,iii)),imag(SIdeal_2(ii,iii)));
            if ~(ii==M && iii==M)
                fprintf(file2,",");
            end
        end
    end
    fprintf(file2,"\n");
    %
    [M,~]       =   size(SIdeal_3);
    fprintf(file3,"%21.14E,",freq(i));
    for ii=1:M
        for iii=1:M
            fprintf(file3,"%21.14E,%21.14E",...
                real(SIdeal_3(ii,iii)),imag(SIdeal_3(ii,iii)));
            if ~(ii==M && iii==M)
                fprintf(file3,",");
            end
        end
    end
    fprintf(file3,"\n");
    T1          =   S2T(SIdeal_1);
    T2          =   S2T(SIdeal_2);
    T3          =   S2T(SIdeal_3);
end
fclose(file1);
fclose(file2);
fclose(file3);
%% 
function[T]=TChain(TA,TB)
[MA,NA]   =   size(TA);
assert((MA==NA)&&(mod(MA,2)==0));
[MB,NB]   =   size(TB);
assert((MB==NB)&&(mod(MB,2)==0));
assert(MA==NB);
N       =   NA;
TA11   	=   TA(1:N/2,1:N/2);
TA12   	=   TA(1:N/2,1+N/2:N);
TA21   	=   TA(1+N/2:N,1:N/2);
TA22   	=   TA(1+N/2:N,1+N/2:N);
TB11   	=   TB(1:N/2,1:N/2);
TB12   	=   TB(1:N/2,1+N/2:N);
TB21   	=   TB(1+N/2:N,1:N/2);
TB22   	=   TB(1+N/2:N,1+N/2:N);
T       =   [ TA11*TB11+TA12*TB21 TA11*TB12+TA12+TB22 ; ...
              TA21*TB11+TA22*TB21 TA21*TB12+TA22*TB22 ];
end
%%
function[T]=S2T(S)
[M,N]   =   size(S);
assert((M==N)&&(mod(M,2)==0));
S11     =   S(1:N/2,1:N/2);
S12     =   S(1:N/2,1+N/2:N);
S21     =   S(1+N/2:N,1:N/2);
S22     =   S(1+N/2:N,1+N/2:N);
T       =   [ inv(S21) -inv(S21)*S22 ; ...
              S11*inv(S21) S12-S11*inv(S21)*S22 ];
end
%%
function[S]=T2S(T)
[M,N]   =   size(T);
assert((M==N)&&(mod(M,2)==0));
T11     =   T(1:N/2,1:N/2);
T12     =   T(1:N/2,1+N/2:N);
T21     =   T(1+N/2:N,1:N/2);
T22     =   T(1+N/2:N,1+N/2:N);
S       =   [ T21*inv(T11) T22-T21*inv(T11)*T12 ; ...
              inv(T11) -inv(T11)*T12 ];
end
%% 
function[S]=CTL(freq,L,R_,G_,L_,C_)
%%
j       =   sqrt(-1);
omega   =   2*pi*freq;
Z_      =   R_+j*omega*L_;
Y_      =   G_+j*omega*C_;
[N,M]   =   size(L_);
assert(N==M);
gamma   =   sqrtm(Z_*Y_);
Z0      =   inv(gamma)*Z_;
Z1      =   50*eye(N);
Z2      =   50*eye(N);
I       =   eye(N);
%%
Gamma_l =   inv(Z2*inv(Z0)+I)*(Z2*inv(Z0)-I);
Gamma_0 =   Gamma_l*expm(-2*L*gamma);
A       =   I+Gamma_0;
B       =   Z1*inv(Z0)*(I-Gamma_0);
C       =   I+Gamma_l;
D       =   Z2*inv(Z0)*(I-Gamma_l);
S11     =   (A-B)*inv(A+B);
S21     =   (C+D)*expm(-gamma*L)*inv(A+B);
S       =   [S11 S21 ; S21 S11];
end
%%
