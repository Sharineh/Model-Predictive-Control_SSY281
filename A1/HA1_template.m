clear;clc;close all;
% Use this template for Home Assignment 1.
% You may add as many variables as you want but DO NOT CHANGE THE NAME OF
% THE VARIABLES WHICH ARE ALREADY IN THE TEMPLATE.
% Also, DO NOT REUSE THE VARIABLES WHICH ARE ALREADY IN THE TEMPLATE FOR
% SOMETHING ELSE.
% The grading will be partially automatized, so it is very important that
% the names of the variables are what the grader script expects to see.
% All the variables whose name must not be changed have a suffix related to
% the number of the question, so just pay extra attention to them.
% Keep also in mind that the variables are expected to have NUMERICAL
% VALUES and not strings of characters.

%% Q1a (discretization of a state-space model)
a = -0.9421;
b = 82.7231;
c = 14.2306;
p = -3.7808;
q = 4.9952;
r = 57.1120;
h = 0.1;
% Contious system matrices
A = [0 1 0 0; b 0 0 a; 0 0 0 1; q 0 0 p];
B = [0; c; 0; r;];
C = [1 0 0 0];
% Discrete Model
fun = @(t) expm(A.*t)*B;
A_1a = expm(A*h);
B_1a = integral(fun, 0, h, 'ArrayValued', true);
C_1a = C;

eig_1a = eig(A_1a);

%% Q1b
% Delayed model
tau = 0.8*h;
sys_1b = ss(A,B,C,0,'InputDelay',tau);
sysd_1b = c2d(sys_1b,h);

Aa_1b = sysd_1b.A;
Ba_1b = sysd_1b.B;
Ca_1b = sysd_1b.C;
eig_1b = eig(Aa_1b);

%% Q2a (Dynamic Programming solution of the LQ problem)
h = 0.01;
A = [1.0041 0.0100 0 0;
    0.8281 1.0041 0 -0.0093;
    0.0002 0.0000 1 0.0098;
    0.0491 0.0002 0 0.9629];

B = [0.0007;
    0.1398;
    0.0028;
    0.5605];

C = [1 0 0 0;
    0 0 1 0];
Q = eye(4);
P_f = 10*eye(4);
R = 1;
% finite horizon problem
i = 0;
P = P_f;
flag = true;
while flag
    i = i+ 1;
    K_d = -(R + B' * P * B)^(-1) * B' * P * A;
    P = (Q + A' * P * A) - (A' * P * B * -K_d);
    if max(abs(eig(A + B*K_d ))) < 1
        flag = false;
    end
end
N_2a = i;
K_2a = K_d;

%% Q2b

% infinite horizon problem
[P_r, K_r ,~] = idare(A, B , Q , R);

P = P_f;

n_2b = 1;
flag = true;
while flag
    P_old = P;
    K_d = -(R + B' * P * B)^(-1) * B' * P * A;
    P = (Q + A' * P * A) - (A' * P * B * -K_d);
    if  norm(P - P_old) < 0.1
        flag = false;
    end
    n_2b = n_2b+1;
end

P_inf_2b = P_r;
P_inf_DP_2b = P;
N_2b = n_2b;

%% Q2c
P = P_r;
n_2c = 0;
flag = true;
while flag
    P_old = P;
    K_d = -(R + B' * P * B)^(-1) * B' * P * A;
    P = (Q + A' * P * A) - (A' * P * B * -K_d);
    n_2c = n_2c+1;
    if  norm(P - P_old) < 0.1
        flag = false;
    end
end

N_2c = n_2c;
K_2c = K_d;

%% Q3a (batch solution of the LQ problem)
A = [1.0041 0.0100 0 0;
    0.8281 1.0041 0 -0.0093;
    0.0002 0.0000 1 0.0098;
    0.0491 0.0002 0 0.9629];

B = [0.0007;
    0.1398;
    0.0028;
    0.5605];

C = [1 0 0 0;
    0 0 1 0];
Q = eye(4);
P_f = 10*eye(4);
R = 1;

gamma = B;
omega = A;
R_bar = 1;

m = 0;
flag = true;
while flag
    j = 0;
    m = m+1;
    Q_bar = blkdiag(kron(eye(m-1),Q),P_f);
    R_bar = kron(eye(m),R);
    gamma = expand_gamma(gamma, m, A, B);
    for i = 1 : m
        omega (i+j:i+j+3,:) = A^i ;
        j=j+3;
    end

    K_vec = -(gamma' * Q_bar * gamma + R_bar)^(-1) * gamma' * Q_bar * omega;
    K_b= K_vec(1,:);
    if max(abs(eig(A + B* K_b ))) < 1
        flag = false;
    end
end

N_3a = size(gamma,2);
K0_3a = K_b;

%% Q4 (Receding horizon control)
R_1 = 1   ; N_1 = 40;
R_2 = 1   ; N_2 = 80;
R_3 = 0.1 ; N_3 = 40;
R_4 = 0.1 ; N_4 = 80;
tf = 2000;
x_0 = [pi/38 ; 0 ; 0 ; 0];
n = size(A,2);
x1_sim = zeros(n,tf) ; x1_sim(:,1) = x_0 ; u1_sim = zeros(tf,1);
x2_sim = zeros(n,tf) ; x2_sim(:,1) = x_0 ; u2_sim = zeros(tf,1);
x3_sim = zeros(n,tf) ; x3_sim(:,1) = x_0 ; u3_sim = zeros(tf,1);
x4_sim = zeros(n,tf) ; x4_sim(:,1) = x_0 ; u4_sim = zeros(tf,1);

for t = 1 :tf
    % controller 1
    [u1,x1] = URHC(A,B,N_1,Q,R_1,P_f,x1_sim(:,t));
    % controller 2
    [u2,x2] = URHC(A,B,N_2,Q,R_2,P_f,x2_sim(:,t));
    % controller 3
    [u3,x3] = URHC(A,B,N_3,Q,R_3,P_f,x3_sim(:,t));
    % controller 4
    [u4,x4] = URHC(A,B,N_4,Q,R_4,P_f,x4_sim(:,t));
    % Simulation data
    x1_sim(:,t+1) = x1; u1_sim(t) = u1(1);
    x2_sim(:,t+1) = x2; u2_sim(t) = u2(1);
    x3_sim(:,t+1) = x3; u3_sim(t) = u3(1);
    x4_sim(:,t+1) = x4; u4_sim(t) = u4(1);
end
x_sim = {x1_sim , x2_sim , x3_sim , x4_sim};
u_sim = {u1_sim,u2_sim,u3_sim,u4_sim};
t = 100;
for i = 1:4
    subplot(2,2,i); hold on; grid on
    for j = 1 : 4
        if i == 4
            plot(0:30,u_sim{j}(1:31),'Linewidth',2)
        end
        if i < 4
            if i == 3
                t = tf;
            end
            plot(0:t,x_sim{j}(i,1:t+1),'Linewidth',2)
        end
    end
    if i < 4
        title(sprintf('X%d for different params',i))
        legend('N=40 & R=1','N=80 & R=1','N=40 & R=0.1','N=80 & R=0.1')
    else
        title('control action U for different params')
        legend('N=40 & R=1','N=80 & R=1','N=40 & R=0.1','N=80 & R=0.1')
    end
end
sgtitle("URHC")
%% Q5 (constrained receding horizon control)
m = size(B,2);
n = length(A);
x2_max = 1; 
u_max  = 8;
x1_sim_c = zeros(n,tf) ; x1_sim_c(:,1) = x_0 ; u1_sim_c = zeros(tf,1);
x2_sim_c = zeros(n,tf) ; x2_sim_c(:,1) = x_0 ; u2_sim_c = zeros(tf,1);
x3_sim_c = zeros(n,tf) ; x3_sim_c(:,1) = x_0 ; u3_sim_c = zeros(tf,1);
x4_sim_c = zeros(n,tf) ; x4_sim_c(:,1) = x_0 ; u4_sim_c = zeros(tf,1);
% same plot as for Q4 but for the constrained problem
for t = 1 :tf
    % controller 1 % CRHC(A,B,N,Q,R,Pf,x0,x_max,u_max,n)
    [z1,~] = CRHC(A,B,N_1,Q,R_1,P_f,x1_sim_c(:,t),x2_max,u_max,n);
    % controller 2
    [z2,~] = CRHC(A,B,N_2,Q,R_2,P_f,x2_sim_c(:,t),x2_max,u_max,n);
    % controller 3
    [z3,~] = CRHC(A,B,N_3,Q,R_3,P_f,x3_sim_c(:,t),x2_max,u_max,n);
    % controller 4
    [z4,~] = CRHC(A,B,N_4,Q,R_4,P_f,x4_sim_c(:,t),x2_max,u_max,n);
    % Simulation data
    x1_sim_c(:,t+1) = z1(1:n); u1_sim_c(t) = z1(n*N_1+1);
    x2_sim_c(:,t+1) = z2(1:n); u2_sim_c(t) = z2(n*N_2+1);
    x3_sim_c(:,t+1) = z3(1:n); u3_sim_c(t) = z3(n*N_3+1);
    x4_sim_c(:,t+1) = z4(1:n); u4_sim_c(t) = z4(n*N_4+1);
end
x_sim_c = {x1_sim_c , x2_sim_c , x3_sim_c , x4_sim_c};
u_sim_c = {u1_sim_c , u2_sim_c , u3_sim_c , u4_sim_c};
%%
t = 100;
for i = 1:4
    subplot(2,2,i); hold on; grid on
    for j = 1 : 4
        if i == 4
            plot(0:30,u_sim_c{j}(1:31),'Linewidth',2)
        end
        if i < 4
            if i == 3
                t = tf;
            end
            plot(0:t,x_sim_c{j}(i,1:t+1),'Linewidth',2)
        end
    end
    if i < 4
        title(sprintf('X%d for different params',i))
        legend('N=40 & R=1','N=80 & R=1','N=40 & R=0.1','N=80 & R=0.1')
    else
        title('control action U for different params')
        legend('N=40 & R=1','N=80 & R=1','N=40 & R=0.1','N=80 & R=0.1')
    end
end
sgtitle("CRHC")
%%
function g = expand_gamma(gamma, m, A , B)
g = gamma;
j = 1;
step = 4;
rows_to_change = (m - 1)*step + 1: m*step;
for i=m :-1: 1
    if i == 1
        res = B;
    else
        res = A^(i-1) * B;
    end
    g(rows_to_change, j) = res;
    j = j + 1;
end
end
%%
function [u,x] = URHC(A,B,N,Q,R,Pf,x0)
P = Pf;
for i = 1 : N
    K_d = -(R + B' * P * B)^(-1) * B' * P * A;
    P = (Q + A' * P * A) - (A' * P * B * -K_d);
end
u = K_d * x0;
x = A * x0 + B * u;
end
%%
function [Z,VN] = CRHC(A,B,N,Q,R,Pf,x0,x_max,u_max,n)
Qbar = blkdiag(kron(eye(N-1),Q),Pf);
Rbar = kron(eye(N),R);
H    = blkdiag(Qbar,Rbar);
f    = [];

% Equality constraints
I    = eye(n);
Aeq1 = kron(eye(N),I)+kron(diag(ones(N-1,1),-1),-A);
Aeq2 = kron(eye(N),-B);
Aeq  = [Aeq1 Aeq2];
beq  = [A*x0;zeros(n*(N-1),1)];

% Inequality constraints
F      = kron([eye(N); zeros(N)], [0 1 0 0;0 -1 0 0]);
G      = kron([zeros(N); eye(N)], [1; -1]);
h      = [x_max*ones(2*N,1); u_max*ones(2*N,1)];

Ain  = [F G];
bin  = h;

% Solve QP
[Z,VN,~,~,~] = quadprog(2*H,f,Ain,bin,Aeq,beq);

end 