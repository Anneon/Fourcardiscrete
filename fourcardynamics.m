close all; clc
clear
%syms  k1 k2 k3 k4
%syms alpha
A = [0 1 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0;
     0 0 0 1 0 0 0 0;
     0 0 0 0 0 0 0 0;
     0 0 0 0 0 1 0 0;
     0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 1;
     0 0 0 0 0 0 0 0];
alpha = 0.01;
S = 1;
B = [0 0 0 0;
     1 0 0 0;
     0 0 0 0;
     0 1 0 0;
     0 0 0 0;
     0 0 1 0;
     0 0 0 0;
     0 0 0 1];
C = [1 -1 0 0 -1 0 0 0;
    0 -1 1 0 0 1 -1 0;
    0 0 1 -1 0 0 1 -1;
    0 0 0 0 -1 0 0 0];
D = zeros(size(C,1),size(B,2));
states = {'x1', 'h1','x2','h2','x3','h3','x4','h4'};
inputs = {'ux1', 'uh1','ux2','uh2','ux3','uh3','ux4','uh4'};
outputs = {'p1','p2','p3','p4'};
sysc = ss(A,B,C,D);
step(sysc);
k = [-1 0 0 0;
     0 -1 0 0;
     0 0 -1 0;
     0 0 0 -1];
UA = ((A+B*alpha*k*C)+B*k*C);
UB = -(B*alpha*k*C);
eig(UA);
%Discrete Time Model
%sampling time
dT = 0.05; %(1/0.05) = 20Hz sampling frequency
sys_d = c2d(sysc,dT,'Zoh'); %zero-order hold

Aol = sys_d.A;
Bol = sys_d.B;
Col = sys_d.C;
Dol = sys_d.D;

k = [-1 0 0 0;
     0 -1 0 0;
     0 0 -1 0;
     0 0 0 -1];
Acl = ((A+B*alpha*k*C)+B*k*C);
Bcl = -(B*alpha*k*C*S);
K=Acl+Bcl;
%Kol = sys_d.K
%Observer design

%parameters G and H
G = eye(4); %because 4 states so 4 disturbances (p = n)
H = zeros(4,4); %4 outputs so 0 matrix of 4 x 4

%covariance matrices, process Q, measurement R
Qcov = diag(0.15*ones(1,4));%Q is 4x4
Rcov = diag(0.05*ones(1,4));%R is 4x4
%Bnew = [Bol*G]
sys_kf = ss(Aol,Bol*G,Col,Dol*H,dT);
%Obtain Kalman gain L and P,assume W and v are uncorrelated
[kest,L,P] = kalman(sys_kf,Qcov,Rcov,0);
%check the value of L that matlab returns
%Compare L_bar to L(Should be equal)
L_gain = (Aol*P*Col')/(Col*P*Col'+Rcov);
Error = norm(abs(L_gain-L));
%Assess the stability
Acb = Aol-L*Col;
eig(Acb);
