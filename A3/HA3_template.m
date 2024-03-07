clc
close all
clear all

%% Q3c
% solve the LP using the formulation shown in (5)
A=[0.4889  0.2939;
   1.0347 -0.7873;
   0.7269  0.8884; 
  -0.3034 -1.1471];
b= [-1.0689;
    -0.8095;
    -2.9443;
     1.4384];
A_in = [A -ones(4,1) ; -A -ones(4,1)];
b_in = [ b;
        -b];
c=[0 ; 0 ; 1]';
sol = linprog(c,A_in,b_in);
z_3c = sol;

%% Q3e
% solve the dual problem derived in Q3d
f = [b;-b];
A_eq = A_in';
b_eq = -[0 0 1]';
lb = zeros(length(f),1);
ub = [];
sol =  linprog(f,[],[],A_eq,b_eq,lb,ub);
mu_3e = sol;

%% Q3f
% Solve the primal problem by using the solution of the dual problem
% obtained in Q3e
active_const= find(mu_3e > 0);
x_primal = A_in(active_const,:)\b_in(active_const);

z_3f = x_primal;

%% Q4a
% solve the QP
A = 0.4;
B = 1;
x0 = 1.5;
H = 0.5 * eye(4);
ub = [5 , 0.5 , 2 , 2];
%ub = [inf , 0.5 , 2 , 2];
lb = [2.5 ,-0.5 , -2 ,-2];
%%lb = [-inf ,-0.5 , -2 ,-2];
A_eq = [1 0 -1 0;  %% x1 = Ax0 + u0 ==> x1- u0 = Ax0
       -A 1 0 -B]; %% x2 = Ax1 + u1 ==> x2 -Ax1 -u1 = 0 
b_eq = [A*x0;0];
[sol,fval,~,~,lambda] = quadprog(2*H,[],[],[],A_eq,b_eq,lb,ub);


lambda_eq = lambda.eqlin;
lambda_ineq = [lambda.upper;lambda.lower];
A_in = [eye(4);-eye(4)];
b_in = [5 , 0.5 , 2 , 2 , - 2.5 , 0.5 , 2 ,2]';

%KKT
% eq constraint h(x*)=0
eq_cond = A_eq * sol - b_eq;
% inq constraint g(x*) = 0
inq_cond = A_in * sol - b_in; 
% mu >=0
mu_cond = lambda_ineq ;
% mu*g(x*) = 0
mu_g_cond = lambda_ineq'*(A_in*sol-b_in) ;
% df + mu dg(x*) + lamda dh(x*) = 0
df = sol;
dg = A_in';
dh = A_eq';
gradient_cond = df + dg *lambda_ineq + dh * lambda_eq ;

x_4a = sol(1:2);
u_4a = sol(3:4);
