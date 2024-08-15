%% Matlab skript for the Paper 2024:
%
% - reproduce the simulations
% - reproduce optimisation results
% - reproduce all figures.png (converted to .eps using LibreOfficeDraw)
% - test average computation time
%

%% Numerical Simulations
%
% Create figures form chapter: Numerical simulations
% 
close all;
[Output1, Output2, time]=fun_1(1,0,'sol',1e4,1e-6),
%
%figure(1); set(gcf,'color','w'); f = gcf; exportgraphics(f,'Figure_6.png');
%figure(2); set(gcf,'color','w'); f = gcf; exportgraphics(f,'Figure_7.png');
%figure(3); set(gcf,'color','w'); f = gcf; exportgraphics(f,'Figure_8_ERD.png');
%

%% Average simulation time:
%
% calculate average computation time of the simulation above
%
T=zeros(1,100); 
for i=1:length(T)
[~, ~, time]=fun_1(1,0.1,'sol',1e4,1e-6);
T(i)=time;
end
mean(T) %= 0.9600
%

%% Optimisation - case1 (epsilon constraint method)
%
% Create Paretofront for single SWRO
%
load("DATA_case1.mat")
% 
% load("DATA_case1.mat")
% save DATA_case1.mat X1 Y1 X1_P Y1_P time1 time1_P
% 
% X1  -  set of operating conditions for epsilon constraint method (2x200)
% Y1  -  system output using X2 (5x40)
%
% X1_P  -  set of operating conditions from paretosearch (2x200)
% Y1_P  -  system output using X2_P (4x200)
%
% time1    - time for epsilon constraint method (40points)
% time1_P  - time for paretosearch solver (200 points)
%
startTime=datetime("now");
A= [-1 1]; b= -0; Aeq=[]; beq=[]; lb = [30;30]; ub = [70;70];
option_mesh = 1e4; option_BVP = 1e-6; option_data = 1;
foptions = optimoptions('fmincon','Display','off','Algorithm','interior-point', 'StepTolerance', 1e-12, 'OptimalityTolerance',1e-4, 'MaxFunEvals',100);
rng default
%
FWmin=linspace(0.0065, 0.3861, 200); 
X0=[linspace(60,69.99,200);linspace(55,50,200)];
%
parfor i=1:200
    epsilon=1;
    minFW=-FWmin(i);
    x0=X0(:,i);
    while epsilon > .1
        [x_neu, minSEC] = fmincon(@(x)fun_1(x,option_data,'SEC',option_mesh,option_BVP),x0   ,A,b,Aeq,beq,lb,ub,@(x)nonlcon(x,'FW' ,option_data,option_mesh,option_BVP,-minFW ),foptions);
        [x0, minFW]     = fmincon(@(x)fun_1(x,option_data,'FW' ,option_mesh,option_BVP),x_neu,A,b,Aeq,beq,lb,ub,@(x)nonlcon(x,'SEC',option_data,option_mesh,option_BVP,-minSEC),foptions);
        y_neu = fun_1(x_neu,option_data,'Pareto',option_mesh,option_BVP),
        y_0   = fun_1(x0,option_data,'Pareto',option_mesh,option_BVP),
        epsilon = min(norm(x0-x_neu), norm(y_neu-y_0));
    end
    X1(:,i)=x0;
    Y1(:,i)=fun_1([X1(:,i)],option_data,'sol',option_mesh,option_BVP);
    i,
end
endTime=datetime("now");      % computation time homePC: 02:22:08 
time1 = endTime - startTime;  % computation time uni-PC: 04:55:45
save DATA_case1.mat X1 Y1 X1_P Y1_P time1 time1_P
%%scatter
%f=figure(1); f.Position = [1200 500 800 500];
scatter(Y1(1,:),Y1(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
legend('Single SWRO unit','Location', 'Northeast');ylabel('FW [m^3/h]','FontSize',16);xlabel('SEC_{net} [kWh/m^3]','FontSize',16);

%% Optimisation - case1 (paretosearch)
%
load("DATA_case1.mat")
% 
% load("DATA_case1.mat")
% save DATA_case1.mat X1 Y1 X1_P Y1_P time1 time1_P
% 
% X1  -  set of operating conditions for epsilon constraint method (2x200)
% Y1  -  system output using X2 (5x40)
%
% X1_P  -  set of operating conditions from paretosearch (2x200)
% Y1_P  -  system output using X2_P (4x200)

startTime=datetime("now");
A= [-1 1]; b= -0; Aeq=[]; beq=[]; lb = [30;30]; ub = [70;70];
option_mesh = 1e4; option_BVP = 1e-6; option_data = 1;
rng default
%
%X0=[linspace(X1(1,1),X1(1,end),200); linspace(X1(2,1),X1(2,end),200)]';
X_init=[X1(1,1:4:end); X1(2,1:4:end)];
X0=repmat(X_init, 1, 4)';
options = optimoptions('paretosearch','ParetoSetSize',200, 'InitialPoints',X0,'Display','iter', 'MaxFunctionEvaluations',10000, 'ParetoSetChangeTolerance',1e-4,'UseParallel', true);
% paretosearch:
[X, Y] = paretosearch(@(x)fun_1(x,2,'Pareto',option_mesh,option_BVP),2,A,b,Aeq,beq,lb,ub,@(x)nonlcon(x,'default'),options); 
%
parfor i=1:200
    X1_P(:,i)=X(i,:)
    Y1_P(:,i)=fun_1([X1_P(:,i)],option_data,'sol',option_mesh,option_BVP);
end
endTime=datetime("now");
time1_P = endTime - startTime;
save DATA_case1.mat X1 Y1 X1_P Y1_P time1 time1_P
% scatter
f=figure(1); f.Position = [1200 500 800 500];
scatter(Y1(1,:),Y1(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
scatter(Y1_P(1,:),Y1_P(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','m'); hold on
grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
legend('epsilon constraint method', 'paretosearch','Location', 'Northeast');ylabel('FW [m^3/h]','FontSize',16);xlabel('SEC_{net} [kWh/m^3]','FontSize',16);

















%% Optimisation - case2 (epsilon constraint method)
%
% Create Paretofront for SWRO+ERD 
%
load("DATA_case2.mat")
% 
% load("DATA_case2.mat")
% save DATA_case2.mat X2 Y2 X2_P Y2_P time2 time2_P
% 
% X2  -  set of operating conditions for epsilon constraint method (2x40)
% Y2  -  system output using X2 (4x40)
%
% X2_P  -  set of operating conditions from paretosearch (2x200)
% Y2_P  -  system output using X2_P (4x200)
%
% time2    - time for epsilon constraint method (40points)
% time2_P  - time for paretosearch solver (200 points)
% 
startTime=datetime("now");
A= [-1 1]; b= -0; Aeq=[]; beq=[]; lb = [30;30]; ub = [70;70];
option_mesh = 1e4; option_BVP = 1e-6; option_data = 2;
foptions = optimoptions('fmincon','Display','off','Algorithm','interior-point', 'StepTolerance', 1e-12, 'OptimalityTolerance',1e-4, 'MaxFunEvals',100);
rng default
%
FWmin=linspace(0.0015, 0.3769, 200);
%X0=X1;
X0=[linspace(37,70,16);linspace(36,50,16)];
%
parfor i=1:200
    epsilon=1;
    minFW=-FWmin(i);
    x0=X0(:,i);
    while epsilon > .001
        [x_neu, minSEC] = fmincon(@(x)fun_1(x,option_data,'SEC',option_mesh,option_BVP),x0   ,A,b,Aeq,beq,lb,ub,@(x)nonlcon(x,'FW' ,option_data,option_mesh,option_BVP,-minFW ),foptions);
        [x0, minFW]     = fmincon(@(x)fun_1(x,option_data,'FW' ,option_mesh,option_BVP),x_neu,A,b,Aeq,beq,lb,ub,@(x)nonlcon(x,'SEC',option_data,option_mesh,option_BVP,-minSEC),foptions);
        y_neu = fun_1(x_neu,option_data,'Pareto',option_mesh,option_BVP),
        y_0   = fun_1(x0,option_data,'Pareto',option_mesh,option_BVP),
        epsilon = max(norm(x0-x_neu), norm(y_neu-y_0));
    end
    X2(:,i)=x0; 
    Y2(:,i)=fun_1([X2(:,i)],option_data,'sol',option_mesh,option_BVP);
    i,
end
endTime=datetime("now");
time2 = endTime - startTime;
save DATA_case2.mat X2 Y2 X2_P Y2_P time2 time2_P
%%scatter
scatter(Y2(1,:),Y2(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
legend('Single SWRO unit','SWRO with ERD','Location', 'Northeast');ylabel('FW [m^3/h]','FontSize',16);xlabel('SEC_{net} [kWh/m^3]','FontSize',16);

%% test 
X0=X1;
%X0=[linspace(37,70,16);linspace(36,50,16)];
A=zeros(1,200);B=A;
parfor i=1:200
    Y0=fun_1([X0(:,i)],2,'sol',1e4,1e-4);
    A(i)=Y0(1);
    B(i)=Y0(2);
end
scatter(A,B,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','g'); hold on



















%% Optimisation - case2 (with paretosearch)
load("DATA_case2.mat")
%
startTime=datetime("now");
A= [-1 1]; b= -.1; Aeq=[]; beq=[]; lb = [30;30]; ub = [70;70];
option_mesh = 1e4; option_BVP = 1e-4; option_data = 2;
%
rng default
% take 200 random data points
% run paretosearch 
%
% endTime
% runTime ...
% scater plot togerther with other 40 points
% define all variables as floating points and see that the


%% Optimisation-2 (red)
close all;clear all;clc
c3=datetime("now");
x0=[69.996767760342830 49.753887785404920];
FWmin=linspace(0.3862, 0.0066, 40);
for i=2:39
clc,
i,
close all; A= [-1 1]; b= -.1; Aeq=[]; beq=[]; lb = [30;30]; ub = [70;70];
foptions = optimoptions('fmincon','Display','iter-detailed','Algorithm','interior-point', 'StepTolerance', 1e-16, 'OptimalityTolerance',1e-4, 'MaxFunEvals',100);
option_mesh = 1e4; option_BVP = 1e-4; option_data = 2;
[x_opt, f_opt]=fmincon(@(x)fun_1(x,2,'SEC',option_mesh,option_BVP),x0,A,b,Aeq,beq,lb,ub,@(x)nonlcon(x,'FW',2,option_mesh,option_BVP, FWmin(i)),foptions);
end
c4=datetime("now"); %von 2 bis 39
save DATA_T2.mat c3 c4
(c4-c3)/38, %00:01:54

%% hybrid paretosearch - old
close all; clear all; clc
rng default % For reproducibility
A= []; b= []; Aeq=[]; beq=[]; lb = [30;10;10;1.05;1.5]; ub = [70;20;20;2;1.5];
% initial points:
load('DATA_Approx.mat');
X0 = X_approx'; %size(X0)= (93,5)
%remove infeasible poitns
z=length(X_approx(1,:));
for i=1:z
    k=0;
    for j=1:4
    if X_approx(j,z+1-i) > ub(j); k=1;end
    if X_approx(j,z+1-i) < lb(j); k=1;end
    end
    if k==1;X0(z+1-i,:)=[]; end
end
%After removing points:
%size(X0)= (61,5)
%
options = optimoptions('paretosearch','ParetoSetSize',200, 'InitialPoints',X0,'Display','iter', 'MaxFunctionEvaluations',5000);
% x = paretosearch(fun,nvars,A,b,Aeq,beq,lb,ub,nonlcon,options)
[X_paretosearch_3, Y_paretosearch_3] = paretosearch(@(x)fun_1(x,3,'Pareto',1e4,1e-4),5,A,b,Aeq,beq,lb,ub,@(x)nonlcon(x,'default'),options); 

%%
load DATA_Pareto.mat
scatter(SEC_Pareto_1,FW_Pareto_1);hold on
scatter(SEC_Pareto_2,FW_Pareto_2)

%% red - Paretosearch
%
%  UseCompletePoll =true --> paretosearch in parallel
%  UseParallel =true --> paralell paretosearch
%
%
close all; clear all; clc
rng default % For reproducibility
%
A= [-1 1]; b= 0; Aeq=[]; beq=[]; lb = [30;30]; ub = [70;70];
%A= []; b= []; Aeq=[]; beq=[]; lb = [30;10;10;1.05;1.5]; ub = [70;20;20;2;1.5];
%
% initial points:
load DATA_Pareto.mat
X0 = [X_Pareto_2';X_Pareto_2';X_Pareto_2';X_Pareto_2';X_Pareto_2'];
%
options = optimoptions('paretosearch','ParetoSetSize',200, 'InitialPoints',X0,'Display','iter', 'MaxFunctionEvaluations',5000);
c5=datetime("now");
% x = paretosearch(fun,nvars,A,b,Aeq,beq,lb,ub,nonlcon,options)
[X, Y] = paretosearch(@(x)fun_1(x,2,'Pareto',1e4,1e-4),2,A,b,Aeq,beq,lb,ub,@(x)nonlcon(x,'default'),options); 
c6=datetime("now");
scatter(-Y(:,1),-Y(:,2));hold on
scatter(SEC_Pareto_2,FW_Pareto_2)
save DATA_T3 c5 c6
%c6-c5 = 00:56:57

%% blue - Paretosearch
%
close all; clear all; clc
rng default % For reproducibility
%
A= [-1 1]; b= 0; Aeq=[]; beq=[]; lb = [30;30]; ub = [70;70];
%A= []; b= []; Aeq=[]; beq=[]; lb = [30;10;10;1.05;1.5]; ub = [70;20;20;2;1.5];
%
% initial points:
load DATA_Pareto.mat
X0 = [X_Pareto_1';X_Pareto_1';X_Pareto_1';X_Pareto_1';X_Pareto_1'];
%
options = optimoptions('paretosearch','ParetoSetSize',200, 'InitialPoints',X0,'Display','iter', 'MaxFunctionEvaluations',5000);
c7=datetime("now");
% x = paretosearch(fun,nvars,A,b,Aeq,beq,lb,ub,nonlcon,options)
[X, Y] = paretosearch(@(x)fun_1(x,1,'Pareto',1e4,1e-4),2,A,b,Aeq,beq,lb,ub,@(x)nonlcon(x,'default'),options); 
c8=datetime("now");
scatter(-Y(:,1),-Y(:,2));hold on
scatter(SEC_Pareto_1,FW_Pareto_1)
save DATA_T4 c7 c8
%c8-c7 = 01:05:44
