clc
close all
clear

%% Q1a
% find terminal weight to guarantee asymptotic stability and plot X_N
A = [1.2 1;
    0 1];
B = [0 ; 1];
Q = eye(2); R = 100; N = 4; X_f = Polyhedron([0 0]);

sys = LTISystem('A',A,'B',B);
sys.x.min = [-15 ; -15];
sys.x.max = [15 ; 15 ];
sys.u.min = -1;
sys.u.max = 1;

sys.x.penalty = QuadFunction(Q);
sys.u.penalty = QuadFunction(R);
sys.x.with('terminalSet');
sys.x.terminalSet = X_f;

X = Polyhedron('lb',sys.x.min,'ub',sys.x.max);
U = Polyhedron('lb',sys.u.min,'ub',sys.u.max);

[P,~,~] = idare(A,B,Q,R);
ctrblty = rank(ctrb(A,B));

X_N = sys.reachableSet('X',X_f,'U',U,'direction','backward','N',N);

Pf_1a = P;
X_N_1a = X_N;

P = QuadFunction(Pf_1a);
sys.x.with('terminalPenalty')
sys.x.terminalPenalty = P;

hold on
plot(X_N,'alpha',0.7)
plot(0,0,'Marker','o','MarkerEdgeColor','K','MarkerFaceColor','auto')
xlabel('$x_1$',Interpreter='latex',FontSize=18)
ylabel('$x_2$',Interpreter='latex',FontSize=18)
legend({'$\mathcal{X}_N$','$\mathcal{X}_f = \mathbf{0}$'},Interpreter="latex",FontSize=18)

%% Q1b
% design 3 different MPC controllers, simulate them and plot the results
% (both for the simulated states and for the actual ones)
x0 = [7 ; -4];
N = {10,15,20};
mpc_ctrl = cell(3,1);
Cloop = cell(3,1);
data  = cell(3,1);
Oloop = cell(3,1);
for i = 1 : 3
    mpc_ctrl{i} = MPCController(sys,N{i});
    Cloop{i} = ClosedLoop(mpc_ctrl{i},sys);
    data{i} = Cloop{i}.simulate(x0,N{i});
    [~,~,Oloop{i}] = mpc_ctrl{i}.evaluate(x0);
    figure
    hold on
    plot(0:N{i},data{i}.X(1,:),"LineWidth",2)
    plot(0:N{i},data{i}.X(2,:),"LineWidth",2)
    plot(0:N{i}-1,data{i}.U,'Color','k',"LineWidth",2)
    plot(0:N{i},Oloop{i}.X(1,:),'LineStyle','--',"LineWidth",2)
    plot(0:N{i},Oloop{i}.X(2,:),'LineStyle','--',"LineWidth",2)
    plot(0:N{i}-1,Oloop{i}.U,'LineStyle','--','Color','k',"LineWidth",2)
    legend({'Actual $x_1$','Actual $x_2$','Actual $u$','Predicted $x_1$','Predicted $x_2$','Predicted $u$'},Interpreter="latex",FontSize=16)
    title(sprintf('MPC controller with horizon N=%d and objective value V=%.2f',N{i},Oloop{i}.cost),"FontSize",17)
    grid on

end


%% Q1c
% find the explicit MPC solution to reach the target set and plot the
% state-space partitions
N = 20;
Ainq = [eye(2); -eye(2)];
binq = 0.01*ones(4,1);
sys.x.with('initialSet')
sys.x.initialSet = X;
X_f = Polyhedron('A',Ainq,'b',binq);
sys.x.with('terminalSet')
sys.x.terminalSet = X_f;

sys.x.with('terminalPenalty')
[P_f,~ , ~]=idare(A,B,Q,R);
sys.x.terminalPenalty = QuadFunction(P_f);

ctrller = MPCController(sys,N);
ex_ctrl = ctrller.toExplicit();

ex_ctrl.partition.plot('wire',false)
hold on
%X.plot('wire', true, 'linestyle', '-.', 'linewidth', 2)
xlabel('$x_1$',Interpreter='latex',FontSize=18)
ylabel('$x_2$',Interpreter='latex',FontSize=18)
title('State-space partitions','Interpreter','latex','FontSize',18)


%% Q1d
% choose a target set Xf such that persistent feasibility is guaranteed for
% all initial state x0 in C_inf
c_inf = sys.invariantSet();
X_f = sys.reachableSet('X',c_inf,'direction','forward') ;
X_f = X_f.intersect(c_inf).minHRep();

Xf_1d = X_f;

c_inf.plot('alpha',0.7)
hold on
X_f.plot('color','b','alpha',0.4)
legend({'$\mathcal{C}_{\infty}$','$\mathcal{X}_f$'},Interpreter="latex",FontSize=18)
xlabel('$x_1$',Interpreter='latex',FontSize=18)
ylabel('$x_2$',Interpreter='latex',FontSize=18)

%% Q2a
% find the matrices for the discretized system
% Model parameters
clc, clear,close all
Ls     = 1.0;
ds     = 0.02;
Js     = 0;
JM     = 0.5;
BetaM  = 0.1;
R      = 20;
KT     = 10;
roh    = 20;
K_teta = 1280.2;
JL     = 50*JM;
BetaL  = 25;

% Sampling interval
h = 0.1;

A = [0                ,  1        ,  0                  ,     0;
    -K_teta/JL        , -BetaL/JL ,  K_teta/(roh*JL)    ,     0;
     0                ,  0        ,  0                  ,     1;
     K_teta/(roh*JM)  ,  0        , -K_teta/(roh^2*JM)  ,    -(BetaM+(KT^2/R)/JM)];

B = [0 ; 
     0 ;
     0 ;
     KT/(R*JM)];

C = [K_teta , 0 , -K_teta/roh , 0];

model = ss(A,B,C,[]);
model = c2d(model,h);



A_2a = model.A;
B_2a = model.B;
C_2a = model.C;

%% Q2b
% design a minimum time (i.e. minimum N) controller that brings the system
% to standstill and plot the predicted states and the output
x0 = [0 ; 2.5 ; 0 ; 75];
sys_2b = LTISystem('A',A_2a,'B',B_2a,'C',C_2a);
sys_2b.u.min = -200;
sys_2b.u.max = 200;
U = Polyhedron('lb',sys_2b.u.min,'ub',sys_2b.u.max);

% x2 = x4 = 0
X_f = Polyhedron('A',[0 1 0 0 ; 0 0 0 1],'b',[0;0]) ;
sys_2b.x.with('terminalSet')
sys_2b.x.terminalSet = X_f;

XN = X_f;
N = 0;
while ~XN.contains(x0)
XN = sys_2b.reachableSet('X',XN,'U',U,'direction','backward');
N=N+1;
end

N_2b = N;
% Predicted system states and output
sys_2b.x.with('penalty')
sys_2b.x.penalty = QuadFunction(eye(4));

sys_2b.x.with('terminalPenalty')
sys_2b.x.terminalPenalty = QuadFunction(eye(4));

sys_2b.u.with('penalty')
sys_2b.u.penalty = QuadFunction(1);

mpc_ctrl = MPCController(sys_2b,N);
[~,~,openloop] = mpc_ctrl.evaluate(x0);
y = C_2a * openloop.X;
minTime = 0:h:h*N;
% plots
subplot(3,1,1)
for i = 1 :4
hold on 
plot(minTime,openloop.X(i,:),'LineWidth',2)
end
legend({'$x_1$','$x_2$','$x_3$','$x_4$'},'Interpreter','latex','FontSize',18)
title('Predicted states','FontSize',18)
xticks(unique(round(get(gca, 'xTick'),1)));
grid on

subplot(3,1,2)
hold on
stairs(minTime(1:end),[openloop.U openloop.U(end)],'LineWidth',2,'Color','#290')
legend({'Input u'},'FontSize',18)
xticks(unique(round(get(gca, 'xTick'),1)));
title('Predicted input','FontSize',18)
grid on

subplot(3,1,3)
plot(minTime,y,'LineWidth',2)
legend({'Output y'},'FontSize',18)
xticks(unique(round(get(gca, 'xTick'),1)));
title('Predicted output','FontSize',18)
grid on

%% Q2c
% design a minimum time controller that brings the system to standstill,
% while having a constrained torsional torque, and plot the predicted states
% and the output
x0 = [0 ; 2.5 ; 0 ; 75];

sys_2c = LTISystem('A',A_2a,'B',B_2a,'C',C_2a);

sys_2c.u.min = -200;
sys_2c.u.max =  200;

U = Polyhedron('lb',-200,'ub',200);
X = Polyhedron('A',[C_2a;-C_2a],'b',[150 ; 150]);

sys_2c.x.with('setConstraint')
sys_2c.x.setConstraint = X;

sys_2c.x.with('penalty')
sys_2c.x.penalty = QuadFunction(eye(4));

sys_2c.x.with('terminalPenalty')
sys_2c.x.terminalPenalty = QuadFunction(eye(4));

sys_2c.u.with('penalty')
sys_2c.u.penalty = QuadFunction(1);

X_f = Polyhedron('A',[0 1 0 0 ; 0 0 0 1],'b',[0;0]) ;
sys_2c.x.with('terminalSet')
sys_2c.x.terminalSet = X_f;


XN = X_f;
N = 0;
while ~XN.contains(x0)
XN = sys_2c.reachableSet('X',XN,'U',U,'direction','backward');
XN = XN.intersect(X);
N  = N+1;
end

N_2c = N;

mpc_ctrl = MPCController(sys_2c,N);
[~,~,openloop] = mpc_ctrl.evaluate(x0);
y = C_2a * openloop.X;
minTime = 0:h:h*N;
%plots
subplot(3,1,1)
for i = 1 :4
hold on 
plot(minTime,openloop.X(i,:),'LineWidth',2)
end
legend({'$x_1$','$x_2$','$x_3$','$x_4$'},'Interpreter','latex','FontSize',18)
title('Predicted states','FontSize',18)
xticks(unique(round(get(gca, 'xTick'),1)));
grid on

subplot(3,1,2)
hold on
stairs(minTime,[openloop.U openloop.U(end)],'LineWidth',2,'Color','#290')
legend({'Input u'},'FontSize',18)
xticks(unique(round(get(gca, 'xTick'),1)));
title('Predicted input','FontSize',18)
grid on

subplot(3,1,3)
plot(minTime,y,'LineWidth',2)
legend({'Output y'},'FontSize',18)
xticks(unique(round(get(gca, 'xTick'),1)));
title('Predicted output','FontSize',18)
grid on
%% Q2d
% design a minimum time controller that brings the system close to
% standstill and plot the predicted states and output

% then, perform a closed loop simulation with a controller having the N you
% just found and plot the states and outputs for 2s
x0 = [0 ; 2.5 ; 0 ; 75];

sys_2d = LTISystem('A',A_2a,'B',B_2a,'C',C_2a);

sys_2d.u.min = -200;
sys_2d.u.max =  200;

U = Polyhedron('lb',-200,'ub',200);
X = Polyhedron('A',[C_2a;-C_2a],'b',[150 ; 150]);

sys_2d.x.with('setConstraint')
sys_2d.x.setConstraint = X;

sys_2d.x.with('penalty')
sys_2d.x.penalty = QuadFunction(eye(4));

sys_2d.x.with('terminalPenalty')
sys_2d.x.terminalPenalty = QuadFunction(eye(4));

sys_2d.u.with('penalty')
sys_2d.u.penalty = QuadFunction(1);


c_inf = sys_2d.invariantSet();
X_target = Polyhedron('lb',[-10 ; -0.01 ; -10 ; -0.01],'ub',[10 ; 0.01 ; 10 ; 0.01]);
X_f = c_inf.intersect(X_target);

sys_2d.x.with('terminalSet')
sys_2d.x.terminalSet = X_f;

XN = X_f;
N = 0;
while ~XN.contains(x0)
XN = sys_2d.reachableSet('X',XN,'U',U,'direction','backward');
XN = XN.intersect(X);
N  = N+1;
end

N_2d = 'minimum horizon length that allows the system to reach almost standstill';


