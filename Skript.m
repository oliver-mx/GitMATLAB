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

%% Overview of computation times:
%
%                               |  home PC   |  uni PC  |
%                               |  i12900k   |
%                               |  16 cores  | 12 cores |
%                               |  3.2 GHz   |
%                               |  64GB RAM  |                      
%
% case 1 - epsilon constraint   |  02:22:08  | 04:55:45 |
% case 1 - paretosearch         |  00:52:50  |    |
%
% case 2 - epsilon constraint   |  06:39:07  |    |
% case 2 - paretosearch         |  00:17:24  |    |
% 
% case 3 - paretosearch         |    |    |
%

%% plot of the Pareto fronts
close all;
load("DATA_case1.mat")
load("DATA_case2.mat")
f=figure(1); f.Position = [1200 500 800 500];
scatter(Y1(1,:),Y1(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
scatter(Y1_P(1,:),Y1_P(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','m'); hold on
scatter(Y2(1,:),Y2(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
scatter(Y2_P(1,:),Y2_P(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','c'); hold on

grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
legend('Case1: \epsilon-constraint method', 'Case1: paretosearch','Case2: \epsilon-constraint method', 'Case2: paretosearch','Location', 'Northeast');
ylabel('FW [m^3/h]','FontSize',16);xlabel('SEC_{net} [kWh/m^3]','FontSize',16);

%% Optimization - case1 (epsilon constraint method)
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
% Y1_P  -  system output using X2_P (5x200)
%
% time1    - time for epsilon constraint method
% time1_P  - time for paretosearch solver
%
startTime=datetime("now");
A= [-1 1]; b= -0; Aeq=[]; beq=[]; lb = [30;30]; ub = [70;70];
option_mesh = 1e4; option_BVP = 1e-6; option_data = 1;
foptions = optimoptions('fmincon','Display','off','Algorithm','interior-point', 'StepTolerance', 1e-12, 'OptimalityTolerance',1e-4, 'MaxFunEvals',100);
rng default
FWmin=linspace(0.0065, 0.3861, 200); 
X0=[linspace(60,69.99,200);linspace(55,50,200)];
%
parfor i=1:200
    epsilon=1;
    minFW=-FWmin(i);
    x0=X0(:,i);
    while epsilon > .01
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

%% Optimisation - case1 (paretosearch)
%
load("DATA_case1.mat")
% 
startTime=datetime("now");
A= [-1 1]; b= -0; Aeq=[]; beq=[]; lb = [30;30]; ub = [70;70];
option_mesh = 1e4; option_BVP = 1e-6; option_data = 1;
rng default
X_init=[X1(1,1:20:end); X1(2,1:20:end)];
X0=repmat(X_init, 1, 20)';
options = optimoptions('paretosearch','ParetoSetSize',200, 'InitialPoints',X0,'Display','iter', 'MaxFunctionEvaluations',10000, 'ParetoSetChangeTolerance',1e-5,'UseParallel', true);
X = paretosearch(@(x)fun_1(x,option_data,'Pareto',option_mesh,option_BVP),2,A,b,Aeq,beq,lb,ub,@(x)nonlcon(x,'default'),options); 
%
parfor i=1:200
    X1_P(:,i)=X(i,:)
    Y1_P(:,i)=fun_1([X1_P(:,i)],option_data,'sol',option_mesh,option_BVP);
end
endTime=datetime("now");        % computation time homePC: 00:52:50
time1_P = endTime - startTime;  %
save DATA_case1.mat X1 Y1 X1_P Y1_P time1 time1_P

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
FWmin=linspace(0.0015, 0.3770, 200);
X0=[[linspace(36.4,69.9,160) linspace(69.9,69.999,40)];[linspace(36.2,60.41,160) linspace(60.41,49.31,40)]];
%
parfor i=1:200
    epsilon=1;
    minFW=-FWmin(i);
    x0=X0(:,i);
    while epsilon > .01
        [x_neu, minSEC] = fmincon(@(x)fun_1(x,option_data,'SEC',option_mesh,option_BVP),x0   ,A,b,Aeq,beq,lb,ub,@(x)nonlcon(x,'FW' ,option_data,option_mesh,option_BVP,-minFW ),foptions);
        [x0, minFW]     = fmincon(@(x)fun_1(x,option_data,'FW' ,option_mesh,option_BVP),x_neu,A,b,Aeq,beq,lb,ub,@(x)nonlcon(x,'SEC',option_data,option_mesh,option_BVP,-minSEC),foptions);
        y_neu = fun_1(x_neu,option_data,'Pareto',option_mesh,option_BVP),
        y_0   = fun_1(x0,option_data,'Pareto',option_mesh,option_BVP),
        epsilon = min(norm(x0-x_neu), norm(y_neu-y_0));
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
legend('Case2: Epsilon constraint method','Location', 'Northeast');ylabel('FW [m^3/h]','FontSize',16);xlabel('SEC_{net} [kWh/m^3]','FontSize',16);

%% Optimisation - case2 (paretosearch)
%
load("DATA_case2.mat")
% 
startTime=datetime("now");
A= [-1 1]; b= -0; Aeq=[]; beq=[]; lb = [30;30]; ub = [70;70];
option_mesh = 1e4; option_BVP = 1e-6; option_data = 2;
rng default
X_init=[X2(1,1:10:end); X2(2,1:10:end)];
X0=repmat(X_init, 1, 10)';
options = optimoptions('paretosearch','ParetoSetSize',200, 'InitialPoints',X0,'Display','iter', 'MaxFunctionEvaluations',10000, 'ParetoSetChangeTolerance',1e-6,'UseParallel', true);
X = paretosearch(@(x)fun_1(x,option_data,'Pareto',option_mesh,option_BVP),2,A,b,Aeq,beq,lb,ub,@(x)nonlcon(x,'default'),options); 
%
parfor i=1:200
    X2_P(:,i)=X(i,:)
    Y2_P(:,i)=fun_1([X2_P(:,i)],option_data,'sol',option_mesh,option_BVP);
end
endTime=datetime("now");        % computation time homePC: 00:52:50
time2_P = endTime - startTime;  %
save DATA_case2.mat X2 Y2 X2_P Y2_P time2 time2_P
%%scatter
f=figure(1); f.Position = [1200 500 800 500];
scatter(Y2_P(1,:),Y2_P(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','c'); hold on
grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
legend('Case2: Paretosearch','Location', 'Northeast');ylabel('FW [m^3/h]','FontSize',16);xlabel('SEC_{net} [kWh/m^3]','FontSize',16);



















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

