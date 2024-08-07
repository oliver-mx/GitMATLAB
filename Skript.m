%% Matlab file for the Paper 2024:
% - reproduce the simulations
% - reproduce optimisation results
% - reproduce all figures.png (converted to .eps using LibreOfficeDraw)
% - test average computation time

%% Numerical Simulations
% figures for SWRO, PRO and the ERD
% original .m file from: C:>User>bav1839>Documents>MATLAB>Counter-current-solver
close all;
[Output1, Output2, time]=fun_1(1,0,'sol',1e6,1e-6),
%figure(1); set(gcf,'color','w'); f = gcf; exportgraphics(f,'Figure_6.png');
%figure(2); set(gcf,'color','w'); f = gcf; exportgraphics(f,'Figure_7.png');
%figure(3); set(gcf,'color','w'); f = gcf; exportgraphics(f,'Figure_8_ERD.png');
%

%% Average simulation time:
% manually disable figures first
T=zeros(1,100);
for i=1:length(T)
[Sim_3b_output1, Sim_3b_output2, time]=fun_1(1,.32,'sol',1e6,1e-6);
T(i)=time;
end
mean(T)
%5.9319

%% Optimisation-1(blue)
c1=datetime("now");
x0=[69.9 51.5];
FWmin=linspace(0.3770, 0.0015, 40);
for i=2:39
clc,
i, 
close all; A= [-1 1]; b= -.1; Aeq=[]; beq=[]; lb = [30;30]; ub = [70;70];
foptions = optimoptions('fmincon','Display','iter-detailed','Algorithm','interior-point', 'StepTolerance', 1e-16, 'OptimalityTolerance',1e-4, 'MaxFunEvals',100);
option_mesh = 1e4; option_BVP = 1e-4; option_data = 1;
[x_opt, f_opt]=fmincon(@(x)fun_1(x,1,'SEC',option_mesh,option_BVP),x0,A,b,Aeq,beq,lb,ub,@(x)nonlcon(x,'FW',1,option_mesh,option_BVP, FWmin(i)),foptions);
end
c2=datetime("now"); %von 2 bis 39
save DATA_T1.mat c1 c2 % 00:02:30

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
