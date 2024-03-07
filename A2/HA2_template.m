clear;clc;close all;
% Use this template for Home Assignment 2.
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

%% Q1a
A = [1.0041 0.0100 0 0;
     0.8281 1.0041 0 -0.0093;
     0.0002 0.0000 1 0.0098;
     0.0491 0.0002 0 0.9629];
B = [0.0007 0.01;
     0.1398 1;
     0.0028 0;
     0.5605 0];
C= [1 0 0 0;
    0 0 1 0];
ys = [pi/18;-pi];
b_eq = [zeros(4,1);ys];
A_eq = [eye(size(A))-A , -B;
        C              , zeros(2)];
sol = linsolve(A_eq,b_eq);

xs_1a = sol(1:4);
us_1a = sol(5:end);


%% Q1b
B = [0.01 ; 1 ; 0 ;0];
A_eq = [eye(size(A))-A, -B];
b_eq = zeros(4,1);
H = zeros(5);
H(1,1) = 2;
H(3,3) = 2;
f=[-2*ys(1);0;-2*ys(2);0;0];
sol = quadprog(H,f,[],[],A_eq,b_eq);

xs_1b = sol(1:4);
us_1b = sol(end);

%% Q1c
B = [0.0007 0.01;
     0.1398 1;
     0.0028 0;
     0.5605 0];
C = [1 0 0 0];
ys = pi/18;
A_eq = [eye(size(A))-A , -B ;C 0 0];
b_eq = [zeros(4,1);ys];
setPoints = linsolve(A_eq,b_eq);
H = eye(6);
H(5,5) = 2;
H(6,6) = 2;
f=[];
sol = quadprog(H,f,[],[],A_eq,b_eq);
xs_1c = sol(1:4);
us_1c = sol(5:end);

%% Q2a
B = [0.0007;
     0.1398;
     0.0028;
     0.5605];
C = [1 0 0 0;
     0 0 1 0];
Cp = [0;1];
Bp = B;
nd1 = 1; Bd1 = zeros(4,1); Cd1 = [1;1];
nd2 = 2; Bd2 = zeros(4,2); Cd2 = eye(2);
nd3 = 2; Bd3 = [zeros(4,1) Bp]; Cd3 = eye(2);

Ae_1_2a = [A Bd1; zeros(1,4) 1];
Be_1_2a = [B;0];
Ce_1_2a = [C Cd1];
detc_1 = rank([eye(size(A))-A -Bd1; C Cd1]) == length(B) + nd1 ;

Ae_2_2a = [A Bd2; zeros(2,5) ones(2,1)];
Be_2_2a = [B;zeros(2,1)];
Ce_2_2a = [C Cd2];
detc_2 = rank([eye(size(A))-A -Bd2; C Cd2]) == length(B) + nd2 ;

Ae_3_2a = [A Bd3; zeros(2,4) eye(2)];
Be_3_2a = [B;zeros(2,1)];
Ce_3_2a = [C Cd3];
detc_3 = rank([eye(size(A))-A -Bd3; C Cd3]) == length(B) + nd3 ;

%% Q2b
% Case 1
Q1 = eye(5); % process noise
R1 = eye(2); % measurement noise
[P1,L1,~,~] = idare(Ae_1_2a',Ce_1_2a',Q1,R1);
K1 = P1*Ce_1_2a'/(Ce_1_2a*P1*Ce_1_2a' + R1);
% Case 3
Q3 = eye(6); % process noise
R3 = eye(2); % measurement noise
[P3,L3,~,~] = idare(Ae_3_2a',Ce_3_2a',Q3,R3);
K3 = P3*Ce_3_2a'/(Ce_3_2a*P3*Ce_3_2a' + R3);
Le_1_2b = K1;
Le_2_2b = 0;
Le_3_2b = K3;

%% Q2c
H = [0 1];

Mss_1_2c = [eye(size(A))-A -B ; H*C 0] \ [Bd1 ; -H*Cd1];
Mss_2_2c = 0;
Mss_3_2c = [eye(size(A))-A -B ; H*C 0] \ [Bd3 ; -H*Cd3];

%% Q2d
N  = 50; M  = 40; R  = 0.1; tf = 1000; n= size(A,1);
A1 = Ae_1_2a; A3 = Ae_3_2a; B1 = Be_1_2a; B3 = Be_3_2a; C1 = Ce_1_2a; C3 = Ce_3_2a;
Q = diag([5 , 2 , 0.5 ,0.1]);
P_f = Q;
x0 = [pi/36 ; 0 ; 0 ; 0];


%%
% Model state vector
x_s1 = zeros(n,tf)  ; x_s1(:,1) = x0 ; u1 = zeros(tf,1);
x_s3 = zeros(n,tf)  ; x_s3(:,1) = x0 ; u3 = zeros(tf,1);
% Observer state vector
x_hat1 = zeros(n+nd1,tf); x_hat1(:,1) = zeros(n+nd1,1);
x_hat3 = zeros(n+nd3,tf); x_hat3(:,1) = zeros(n+nd3,1);
% Output vector
ys = zeros(n,tf);
y_meas1 = zeros(size(C,1),tf);
y_meas3 = zeros(size(C,1),tf);
% disturbance function
p = @(k) (k > 50) * 0.2;

for k = 1 :tf
    % System 1
    d_hat1 = x_hat1(end,k);
    var1 = Mss_1_2c * d_hat1;
    xs_1 = var1(1:4); us_1 = var1(5);
    delta_x1 = x_hat1(1:n,k)- xs_1;

    [z1,~] = CRHC(A,B,N,M,Q,R,P_f,delta_x1,n);

    delta_u1 = z1(n*N+1);

    u1(k) = delta_u1 + us_1;

    y_meas1(:,k) = C * x_s1(:,k);

    x_hat1_corr = x_hat1(:,k) + K1*(y_meas1(:,k) - C1 * x_hat1(:,k));
    x_hat1(:,k+1) = A1*x_hat1_corr + B1*u1(k);

    x_s1(:,k+1) = A * x_s1(:,k) + B * u1(k) + Bp*p(k);

    % System 3

    d_hat3 = x_hat3(end-1:end,k);
    var3 = Mss_3_2c * d_hat3;
    xs_3 = var3(1:4); us_3 = var3(5);
    delta_x3 = x_hat3(1:n,k)- xs_3;

    [z3,~] = CRHC(A,B,N,M,Q,R,P_f,delta_x3,n);

    delta_u3 = z3(n*N+1);

    u3(k) = delta_u3 + us_3;
    y_meas3(:,k) = C * x_s3(:,k);

    x_hat3_corr = x_hat3(:,k) + K3*(y_meas3(:,k) - C3 * x_hat3(:,k));
    x_hat3(:,k+1) = A3*x_hat3_corr + B3*u3(k);

    x_s3(:,k+1) = A * x_s3(:,k) + B * u3(k) + Bp*p(k);
end
%%
% provide plots showing how the states evolve for all the observable systems
figure;
for i = 1:5
    subplot(5,1,i); hold on; grid on
    t=tf;
        if i == 4
            p1 = plot(0:t-1,u1(1:t),'Linewidth',2);
            p2 = yline(0,'r--');
            p3 = plot(50,u1(50),"+k",MarkerSize=18);
            title("Control input u")
            legend([p1,p2,p3],"Control input","Set point",'Disturbances entry')
        elseif i < 4
            if i == 3
                t =tf; 
            end
            p1 = plot(0:t-1,x_s1(i,1:t),'Linewidth',2);
            p2 = yline(0,'r--');
            p3 = plot(50,x_s1(i,50),"+k",MarkerSize=18);
            title(sprintf("x%d",i))
            legend([p1 p2 p3],"State trajectory","Set point",'Disturbances entry')
        else
            p1 = plot(0:t-1,x_hat1(5,1:t),'Linewidth',2);
            p2 = yline(0,'r--');
            p3 = plot(50,x_hat1(5,50),"+k",MarkerSize=18);
            title(sprintf("Disturbance state estimation x%d",i))
            legend([p1 p2 p3],"State trajectory","Set point",'Disturbances entry')
        end
end
sgtitle("Model 1")

figure;
for i = 1:6
    subplot(6,1,i); hold on; grid on
    t=tf;
        if i == 4
            t = tf;
            p1 = plot(0:t-1,u3(1:t),'Linewidth',2);
            p2 = yline(0,'r--');
            p3 = plot(50,u3(50),'+k',MarkerSize=18);
            title("Control input u")
            legend([p1 p2 p3],"control input","Set point",'Disturbances entry')
        elseif i < 4
            if i ==3
                t= tf;
            end
            p1 = plot(0:t-1,x_s3(i,1:t),'Linewidth',2);
            p2 = yline(0,'r--');
            p3 = plot(50,x_s3(i,50),'+k','MarkerSize',18);
            title(sprintf("x%d",i))
            legend([p1 p2 p3],"State trajectory","Set point",'Disturbances entry')
        else
            p1 = plot(0:t-1,x_hat3(i,1:t),'Linewidth',2);
            p2 = yline(0,'r--');
            p3 = plot(50,x_hat3(i,50),"+k",MarkerSize=18);
            title(sprintf("Disturbance state estimation x%d",i))
            legend([p1 p2 p3],"State trajectory","Set point",'Disturbances entry')
        end
end

sgtitle("Model 3")

%%
function [Z,VN] = CRHC(A,B,N,M,Q,R,Pf,x0,n)
Qbar = blkdiag(kron(eye(N-1),Q),Pf);
Rbar = kron(eye(M),R);
H    = blkdiag(Qbar,Rbar);
f    = [];

% Equality constraints
I    = eye(n);
Aeq1 = kron(eye(N),I)+kron(diag(ones(N-1,1),-1),-A);
Aeq2 = kron(eye(M),-B);
Aeq_aux  = repmat(Aeq2(end-3 : end,:),N-M,1);
Aeq2 = [Aeq2 ; Aeq_aux];
Aeq  = [Aeq1 Aeq2];
beq  = [A*x0;zeros(n*(N-1),1)];

% Solve QP
[Z,VN,~,~,~] = quadprog(2*H,f,[],[],Aeq,beq);

end 