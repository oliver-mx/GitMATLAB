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
%                               |  home PC   |  uni PC    |
%                               |  i12900k   |  i135000k  |
%                               |  16 cores  |  14 cores  |
%                               |  3.2 GHz   |  2.5 GHz   |
%                               |  64GB RAM  |  16GB RAM  |                  
%
% case 1 - epsilon constraint   |  02:22:08  |  04:55:45  |
% case 1 - paretosearch         |  00:52:50! |  01:07:44  |
%
% case 2 - epsilon constraint   |  06:39:07  |  !  |
% case 2 - paretosearch         |  00:17:24! |  00:46:28  |
% 
% case 3 - paretosearch         |  !  |  !  |
%

%% plot of the Pareto fronts
close all;
load("DATA_case1.mat")
f=figure(1); f.Position = [1200 500 800 500];
scatter(Y1(1,:),Y1(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
scatter(Y1_P(1,:),Y1_P(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','m'); hold on
scatter(Y1_R(1),Y1_R(2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','y'); hold on
load("DATA_case2.mat")
scatter(Y2(1,:),Y2(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
scatter(Y2_P(1,:),Y2_P(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','c'); hold on
scatter(Y2_R(1),Y2_R(2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','y'); hold on
load("DATA_case3.mat")
%scatter(Y3_R(1),Y3_R(2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','y'); hold on
scatter(-2.2348,0.3440,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','y'); hold on %/bei ~ 10/3 Längenverhältnis



grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
legend('Case1: \epsilon-constraint method', 'Case1: paretosearch','','Case2: \epsilon-constraint method', 'Case2: paretosearch','', '','Location', 'Northeast');
ylabel('FW [m^3/h]','FontSize',16);xlabel('SEC_{net} [kWh/m^3]','FontSize',16);

%% Case1: epsilon-constraint method
%
load("DATA_case1.mat")
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
endTime=datetime("now");     
time1 = endTime - startTime;
save DATA_case1.mat X1 Y1 X1_P Y1_P X1_R Y1_R time1 time1_P

%% Case1: Max Revenue 
%
load("DATA_case1.mat")
[a1, b1]=max(Y1(3,:));
A= [-1 1]; b= -0; Aeq=[]; beq=[]; lb = [30;30]; ub = [70;70];
option_mesh = 1e4; option_BVP = 1e-6; option_data = 1; 
foptions = optimoptions('fmincon','Display','iter','Algorithm','interior-point', 'StepTolerance', 1e-12, 'OptimalityTolerance',1e-4, 'MaxFunEvals',100);
rng default
X1_R= fmincon(@(x)fun_1(x,option_data,'Rev',option_mesh,option_BVP),X1(:,b1),A,b,Aeq,beq,lb,ub,@(x)nonlcon(x,'default'),foptions);
Y1_R=fun_1(X1_R,option_data,'sol',option_mesh,option_BVP);
save DATA_case1.mat X1 Y1 X1_P Y1_P X1_R Y1_R time1 time1_P

%% Case1: paretosearch
%
load("DATA_case1.mat")
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
endTime=datetime("now");       
time1_P = endTime - startTime;
save DATA_case1.mat X1 Y1 X1_P Y1_P X1_R Y1_R time1 time1_P

%% Case2: epsilon-constraint method
%
load("DATA_case2.mat")
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
        epsilon = max(norm(x0-x_neu), norm(y_neu-y_0));
    end
    X2(:,i)=x0; 
    Y2(:,i)=fun_1([X2(:,i)],option_data,'sol',option_mesh,option_BVP);
    i,
end
endTime=datetime("now");
time2 = endTime - startTime;
save DATA_case2.mat X2 Y2 X2_P Y2_P X2_R Y2_R time2 time2_P

%% Case2: Max Revenue 
%
load("DATA_case2.mat")
[a2, b2]=max(Y2(3,:));
A= [-1 1]; b= -0; Aeq=[]; beq=[]; lb = [30;30]; ub = [70;70];
option_mesh = 1e4; option_BVP = 1e-6; option_data = 2; 
foptions = optimoptions('fmincon','Display','iter','Algorithm','interior-point', 'StepTolerance', 1e-12, 'OptimalityTolerance',1e-4, 'MaxFunEvals',100);
rng default
X2_R=fmincon(@(x)fun_1(x,option_data,'Rev',option_mesh,option_BVP),X2(:,b2),A,b,Aeq,beq,lb,ub,@(x)nonlcon(x,'default'),foptions);
Y2_R=fun_1(X2_R,option_data,'sol',option_mesh,option_BVP);
save DATA_case2.mat X2 Y2 X2_P Y2_P X2_R Y2_R time2 time2_P

%% Case2: paretosearch
%
load("DATA_case2.mat")
% 
startTime=datetime("now");
A= [-1 1]; b= -0; Aeq=[]; beq=[]; lb = [30;30]; ub = [70;70];
option_mesh = 1e4; option_BVP = 1e-6; option_data = 2;
rng default
X_init=[X2(1,1:5:end); X2(2,1:5:end)];
X0=repmat(X_init, 1, 5)';
options = optimoptions('paretosearch','ParetoSetSize',200, 'InitialPoints',X0,'Display','iter', 'MaxFunctionEvaluations',10000, 'ParetoSetChangeTolerance',1e-6,'UseParallel', true);
X = paretosearch(@(x)fun_1(x,option_data,'Pareto',option_mesh,option_BVP),2,A,b,Aeq,beq,lb,ub,@(x)nonlcon(x,'default'),options); 
%
parfor i=1:200
    X2_P(:,i)=X(i,:)
    Y2_P(:,i)=fun_1([X2_P(:,i)],option_data,'sol',option_mesh,option_BVP);
end
endTime=datetime("now");     
time2_P = endTime - startTime; 
save DATA_case2.mat X2 Y2 X2_P Y2_P X2_R Y2_R time2 time2_P

%% Case3: Max Revenue 
%
%load("DATA_case3.mat")
A= [0 0 1 -1 0]; b= -.01; Aeq=[]; beq=[]; lb = [.5;30;10;10;1]; ub = [10;70;20;20;5];
option_mesh = 1e4; option_BVP = 1e-6;option_data = .3;
%foptions = optimoptions('fmincon','Display','iter','Algorithm','interior-point', 'StepTolerance', 1e-12, 'OptimalityTolerance',1e-4, 'MaxFunEvals',5000);
optsm = optimoptions("patternsearch",Algorithm="nups" ,Display="iter",PlotFcn="psplotbestf", ConstraintTolerance=1e-4, MaxIterations=1000, StepTolerance=1e-12, UseParallel=true);
rng default % ,Algorithm = {"classic"} | "nups" | "nups-gps" | "nups-mads"
%
%x0 = [ 1.3275   70.0000   11.9213   13.6448    2.0934]; %for nups-mads
x0 = [1.869313646668921 69.999179474445380 13.661473058241283 16.206591230548390 2.718807588032991]; %for nups-gps

%
%X3_R = fmincon(@(x)fun_1(x,option_data,'Rev',option_mesh,option_BVP),x0,A,b,Aeq,beq,lb,ub,@(x)nonlcon(x,'default'),foptions);
X3_R = patternsearch(@(x)fun_1(x,option_data,'Rev',option_mesh,option_BVP),x0,A,b,Aeq,beq,lb,ub,@(x)nonlcon(x,'default'),optsm);
Y3_R=fun_1(X3_R,option_data,'sol',option_mesh,option_BVP);
display(X3_R)
display(Y3_R)

% best X with nup-mads:
%  1.8693   69.9992   12.1019   14.7601    2.7188
% output:
% -2.4087    0.3406    0.4413   41.9332    4.4192

% best X with nup-mads:
% 1.8693   69.9992   13.6615   16.2066    2.7188
% output:
% -2.3508    0.3368    0.4461   42.5031    4.2593

%save DATA_case3 X3_R Y3_R

%% test locally around best point so far

a1=      1.1.*ones(1,2000) + .2.*rand(1,2000);
a2=       70.*ones(1,2000);
a4=       20*ones(1,2000);
a3= a4 - 1.1.*ones(1,2000) -  .3.*rand(1,2000);
a5=     1.0001.*ones(1,2000) + .02.*rand(1,2000);

X_test=[a1; a2; a3; a4; a5];

parfor i=1:2000
i,
Y_test(:,i)=fun_1(X_test(:,i), .3 ,'sol', 1e4 , 1e-6);
end

%%

scatter(Y_test(1,:),Y_test(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','g'); hold on

R=Y_test(3,:)
R = R(~isnan(R))';
R = R(~isinf(R))';
[a,b]=max(R)
c=find(Y_test(3,:)== a)


Y_test(:,c)'
beep,
X_test(:,c)

%
% testing with 2000 random pertubations
% 1.701334230850316 70 17.303900763215214 19.518212394164596 2.515543290319345
% -2.3210  0.3430  0.4594  41.8200  2.6930
%

%
% 2nd test
% 1.645379340650449 70 17.918231682265990 19.890375014904670 2.018569892899560
% -2.2546    0.3390    0.4654   42.4818    3.2361
%


%
% 3rd test
% 1.515502010177877 70 18.144811408171233 19.888390283989512 1.214990800883215
% -2.2550    0.3431    0.4709   41.9577    9.4891
%

%
% 4th test
% 1.581025550521746 70 18.086408450978364 20 1.109362628068349
% -2.2325    0.3427    0.4742   42.2457   40.7659
%

%
% 5th
% 1.414591373229294 70 18.465466683045230 20 1.051796875666936
% -2.2352    0.3437    0.4751   42.0992   47.1326
%

%
% 6th
% 1.217202844329564 70 18.869106506608680 20 1.019759965650840
% -2.2348    0.3440    0.4756   42.0134   51.5381
%



































%% hybrid system optimal length:
A= [0 0 1 -1 0]; b= -.01; Aeq=[]; beq=[]; lb = [.5;30;10;10;1]; ub = [10;70;20;20;5];
option_mesh = 1e4; option_BVP = 1e-6;option_data = .3;
%foptions = optimoptions('fmincon','Display','final','Algorithm','interior-point', 'StepTolerance', 1e-12, 'OptimalityTolerance',1e-4, 'MaxFunEvals',5000);
optsm = optimoptions("patternsearch",Algorithm="nups-mads",Display="iter",PlotFcn="psplotbestf", ConstraintTolerance=1e-4, MaxIterations=1000, StepTolerance=1e-12, UseParallel=true);
rng default %{"classic"} | "nups" | "nups-gps" | "nups-mads"
%
%x0 = [ 1.3275   70.0000   11.9213   13.6448    2.0934];
x0 = [2 69.9999 17.5 19.3 2.5];
%
%[x_opt, rev_opt] = fmincon(@(x)fun_1(x,option_data,'Rev',option_mesh,option_BVP),x0,A,b,Aeq,beq,lb(:,i),ub(:,i),@(x)nonlcon(x,'default'),foptions);
X3_R = patternsearch(@(x)fun_1(x,option_data,'Rev',option_mesh,option_BVP),x0,A,b,Aeq,beq,lb,ub,@(x)nonlcon(x,'default'),optsm);
Y3_R=fun_1(X3_R,option_data,'sol',option_mesh,option_BVP);
close all,
display(X3_R)

% best X with nups-gps:
%  1.3275   70.0000   11.9263   13.6482    2.0934
% output:
% -2.6397    0.3579    0.7268   38.0479    1.3155

















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

