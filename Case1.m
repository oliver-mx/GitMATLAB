%% Pareto Front for the simple SWRO system
% press [ctrl]+[enter] to run code sections
addpath('Input_DATA','Scaled_model','Unscaled_model','Output_DATA')

%% gather Initial data
clc;close all
load Output_DATA/Test2_Output % Simulation results
% single SWRO
X=[X0' Y1(:,1:2)];
Fmin=min(X(:,4));[Fmax, Fmax_index]=max(X(:,4));F=linspace(Fmin,Fmax,200);X1_init=zeros(2,200);X1_init(:,end)=X(Fmax_index,1:2);Y1_init=zeros(2,200);Y1_init(:,end)=X(Fmax_index,3:4);
for i= 1:199
    f1=find(X(:,4)>=F(i) & X(:,4)<F(i+1));X1=X(f1,:);[a,b]=max(X1(:,3));X1_init(:,i)=X1(b,1:2);Y1_init(:,i)=X1(b,3:4);
end
% SWRO+ERD
X=[X0' Y2(:,1:2)];
Fmin=min(X(:,4));[Fmax, Fmax_index]=max(X(:,4));F=linspace(Fmin,Fmax,200);X2_init=zeros(2,200);X2_init(:,end)=X(Fmax_index,1:2);Y2_init=zeros(2,200);Y2_init(:,end)=X(Fmax_index,3:4);
for i= 1:199
    f1=find(X(:,4)>=F(i) & X(:,4)<F(i+1));X1=X(f1,:);[a,b]=max(X1(:,3));X2_init(:,i)=X1(b,1:2);Y2_init(:,i)=X1(b,3:4);
end
% plot
f=figure(1);f.Position=[1000 727.6667 1207 510.0000];tiledlayout(1,2);nexttile
scatter(Y1_init(1,:),Y1_init(2,:),'m');hold on; scatter(Y2_init(1,:),Y2_init(2,:),'c');hold on;
xlim([-5.5 0]);ylim([0.1 1.45]);grid on;title('Paretosearch - initial data','FontSize',14);xlabel('SEC_{net} [kWh/m^3]','FontSize',12);ylabel('FW [m^3/h]','FontSize',12)
nexttile;
scatter(X1_init(1,:),X1_init(2,:),'m');hold on; scatter(X2_init(1,:),X2_init(2,:),'c');hold on;
xlim([29.4 70.6]);ylim([29.4 70.6]);grid on;title('operating pressures','FontSize',14);xlabel('P_d(0) [bar]','FontSize',12);ylabel('P_d(L) [bar]','FontSize',12)

%% Paretosearch - Case1
clc;close all
load("Output_DATA/DATA_Case_1.mat");clc
disp('___________________________________________________________')
disp('Starting Pareto front calculation for the single SWRO unit.')
startTime=datetime("now"); 
A= [-1 1; 1 -1]; b= [0; 3.4]; Aeq=[]; beq=[]; lb = [30;30]; ub = [70;70];
option_mesh = 1e4; option_BVP = 1e-4; option_data = 1;
% paretosearch
options = optimoptions('paretosearch','ParetoSetSize',200, 'InitialPoints',X1_init','Display','iter','MaxFunctionEvaluations',10000, 'ParetoSetChangeTolerance',1e-6,'UseParallel', true);
X = paretosearch(@(x)fun_unscaled(x,option_data,'Pareto',option_mesh,option_BVP),2,A,b,Aeq,beq,lb,ub,[],options);
% evaluate optimal points
X1_pareto=X';Y1_pareto=zeros(200,28);
parfor i=1:200
    Y1_pareto(i,:)=fun_unscaled(X(i,:),option_data,'sol',option_mesh,option_BVP);
end
endTime=datetime("now");
t1_pareto = endTime - startTime;
save Output_DATA/DATA_Case_1.mat X1_init Y1_init X1_pareto Y1_pareto t1_pareto X0 n
disp('Finished calculating the Pareto front for the single SWRO unit.')
disp('_______________________________________________________________')
system('git add .'); system('git commit -m "Case1-paretosearch"');system('git push https://github.com/oliver-mx/GitMATLAB.git');
%

%% Plot the pareto front 
clc;close all
load("Output_DATA/DATA_Case_1.mat");
% plot
f=figure(1);f.Position=[1000 727.6667 1207 510.0000];tiledlayout(1,2);nexttile
scatter(Y1_pareto(:,1),Y1_pareto(:,2),'red');hold on;
xlim([-5.5 0]);ylim([0.1 1.45]);grid on;title('Pareto front','FontSize',14);xlabel('SEC_{net} [kWh/m^3]','FontSize',12);ylabel('FW [m^3/h]','FontSize',12);legend('Case1: SWRO','Location', 'best');
nexttile;
scatter(X1_pareto(1,:),X1_pareto(2,:),'red');hold on;
xlim([29.4 70.6]);ylim([29.4 70.6]);grid on;title('Pareto front - operating pressures','FontSize',14);xlabel('P_d(0) [bar]','FontSize',12);ylabel('P_d(L) [bar]','FontSize',12);legend('Case1: SWRO','Location', 'best');




















%% Paretosearch - Case2
load Output_DATA/Test2_Output % data from simulation
load("Output_DATA/DATA_Case_2.mat");clc
disp('_____________________________________________________________')
disp('Starting Pareto front calculation for the SWRO unit with ERD.')
% Simulations:
% fun_unscaled([0; X0(:,i)],.2,'sol',1e4,1e-4);
% Now:
% % fun_unscaled([X0(:,i)],2,'sol',1e4,1e-4);
startTime=datetime("now");
% constraints
A= [-1 1; 1 -1]; b= [0; 3.4]; Aeq=[]; beq=[]; lb = [30;30]; ub = [70;70];
option_mesh = 1e4; option_BVP = 1e-4; option_data = 2;
%initial
X1_init=repmat(X0(:,I2),1,10);X2_inti=X2_init(:,1:200);
% paretosearch
options = optimoptions('paretosearch','ParetoSetSize',200, 'InitialPoints',X2_init','Display','iter','MaxFunctionEvaluations',10000, 'ParetoSetChangeTolerance',1e-6,'UseParallel', true);
X = paretosearch(@(x)fun_unscaled(x,option_data,'Pareto',option_mesh,option_BVP),2,A,b,Aeq,beq,lb,ub,[],options);
% evaluate optimal points
X2_pareto=X';Y2_pareto=zeros(200,28);
%
parfor i=1:200
    Y2_pareto(i,:)=fun_unscaled(X(i,:),option_data,'sol',option_mesh,option_BVP);
end
endTime=datetime("now");
t2_pareto = endTime - startTime;
save Output_DATA/DATA_Case_2.mat X2_init Y2_init X2_pareto Y2_pareto t2_pareto X0 n
disp('Finished calculating the Pareto front for the SWRO+ERD case.')
disp('_________________________________________________________________')
%
%system('git status');
%system('git add .');
%system('git commit -m "Paretosearch Case 1+2"');
%system('git push https://github.com/oliver-mx/GitMATLAB.git');
%system('shutdown /s /t 30');

%% plot in comparison with initial data
close all;
load("Output_DATA/Test2_Output.mat")
load("Output_DATA/DATA_Case_1.mat")
load("Output_DATA/DATA_Case_2.mat")
% Paretosearch
scatter(Y1_pareto(:,1),Y1_pareto(:,2),'r','filled');hold on; title('simple SWRO'); xlabel('SEC_{net} [kWh/m^3]'); ylabel('FW [m^3/h]');
%
I2=[2149 6654 3224 4819 5901 471 9614 3118 5718 5124 3653 4825 9337 8264 6711 1924 3250 4957 4439 3248 7393 6568 8272];
scatter(Y2(I2,1),Y2(I2,2),'k','filled');hold on; title('simple SWRO'); xlabel('SEC_{net} [kWh/m^3]'); ylabel('FW [m^3/h]');
% Paretosearch
scatter(Y2_pareto(:,1),Y2_pareto(:,2),'b','filled');hold on; title('simple SWRO'); xlabel('SEC_{net} [kWh/m^3]'); ylabel('FW [m^3/h]');
xlim([-6 0]);ylim([0.05 1.5]);view(2)















































%% find initial data
% based on simulation results from "Basic_test_2.m"
load Output_DATA/Basic_Test_Output
close all; clc;
f=figure(1); f.Position= [823.6667 167.6667 989.3333 483.3333];
Fmax=max(Y1(:,2));Fmin=min(Y1(:,2));
F=linspace(Fmin,Fmax,200);
Y1_init=zeros(200,2);X1_init=zeros(200,2);
for i=1:200
    index=find(Y1(:,2)>=[F(i)]);
    a=max(Y1(index,1));
    Y1_init(i,1)=a;
    b=find(Y1(:,1)==a);
    Y1_init(i,2)=Y1(b,2);
    X1_init(i,:)=X0(:,b);
end
%scatter3(Y1(:,1),Y1(:,2),1:n,'r','filled');hold on; title('simple SWRO'); xlabel('SEC_{net} [kWh/m^3]'); ylabel('FW [m^3/h]');
scatter3(Y1_init(:,1),Y1_init(:,2),1:200,'k','filled');hold on; 
xlim([-6 -1]);ylim([.1 1.4]);view(2)

%% Paretosearch - Case1
load("Output_DATA/DATA_Case_1.mat");clc
disp('___________________________________________________________')
disp('Starting Pareto front calculation for the single SWRO unit.')
%
startTime=datetime("now");
% constraints
A= [-1 1; 1 -1]; b= [0; 3.4]; Aeq=[]; beq=[]; lb = [30;30]; ub = [70;70];
option_mesh = 1e4; option_BVP = 1e-4; option_data = 1;
% paretosearch
options = optimoptions('paretosearch','ParetoSetSize',200, 'InitialPoints',X1_init','Display','iter','MaxFunctionEvaluations',10000, 'ParetoSetChangeTolerance',1e-6,'UseParallel', true);
X = paretosearch(@(x)fun_unscaled(x,option_data,'Pareto',option_mesh,option_BVP),2,A,b,Aeq,beq,lb,ub,[],options);
% evaluate optimal points
X1_pareto=X';Y1_pareto=zeros(200,18);
parfor i=1:200
    Y1_pareto(i,:)=fun_unscaled(X(i,:),option_data,'sol',option_mesh,option_BVP);
end
endTime=datetime("now");
t1_pareto = endTime - startTime;
save Output_DATA/DATA_Case_1.mat X1_init X1_pareto Y1_pareto t1_pareto
disp('Finished calculating the Pareto front for the single SWRO unit.')
disp('_______________________________________________________________')
%
%% Compare initial data to Pareto curve
scatter(Y1_init(:,1),Y1_init(:,2),'k','filled');hold on;xlabel('SEC_{net} [kWh/m^3]'); ylabel('FW [m^3/h]');
scatter(Y1_pareto(:,1),Y1_pareto(:,2),'r','filled');hold on;
xlim([-10 0]);ylim([0 2.2]);

































%% Paretosearch - Case2
load("Output_DATA/DATA_Case_2.mat");clc
disp('Starting Pareto front calculation for the SWRO unit with ERD.')
%
startTime=datetime("now");
% initial data
rng default; 
X_init = 31.*ones(1,200) + 39.*rand(2,200);
X_init(2,:) = X_init(1,:) - .1.*ones(1,200) - 3.3.*rand(1,200);
X_init(2,:) = max(X_init(2,:), 30*ones(1,200));
% constraints
A= [-1 1; 1 -1]; b= [0; 3.4]; Aeq=[]; beq=[]; lb = [30;30]; ub = [70;70];
option_mesh = 1e4; option_BVP = 1e-3; option_data = -.6; % <---IMPORTANT
% paretosearch
options = optimoptions('paretosearch','ParetoSetSize',200, 'InitialPoints',X_init','Display','iter','MaxFunctionEvaluations',10000, 'ParetoSetChangeTolerance',1e-4,'UseParallel', true);
X = paretosearch(@(x)fun_scaled(x,option_data,'Pareto',option_mesh,option_BVP),2,A,b,Aeq,beq,lb,ub,@(x)nonlcon(x,'default'),options);
% evaluate optimal points
X2c_pareto=X';Y2c_pareto=zeros(200,18);
parfor i=1:200
    Y2c_pareto(i,:)=fun_scaled(X(i,:),option_data,'sol',option_mesh,option_BVP);
end
endTime=datetime("now");
t2c_pareto = endTime - startTime;
save Output_DATA/DATA_Case_2.mat X2a_pareto X2b_pareto X2c_pareto Y2a_pareto Y2b_pareto Y2c_pareto t2a_pareto t2b_pareto t2c_pareto 

clc
disp('everything went sucsessful')
disp('please make sure to save all data')

%% plot the three Pareto curves:
close all;
load("Output_DATA/DATA_Case_1.mat");
f=figure(1); f.Position = [1200 500 800 500];
%
scatter(Y1_pareto(:,1),Y1_pareto(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
grid on; xlim([-4.4 -.5]); ylim([0 2.2]);view(2);
%%
scatter(Y2a_pareto(:,1),Y2a_pareto(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','#00008B'); hold on
scatter(Y2b_pareto(:,1),Y2b_pareto(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','#0096FF'); hold on
scatter(Y2c_pareto(:,1),Y2c_pareto(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
%
grid on; xlim([-4.4 -.5]); ylim([0 2.2]);view(2);
legend('Case1: ideal permeate flows', 'Case1: ICP', 'Case1: ICP+ECP','Case2: ideal permeate flows', 'Case2: ICP', 'Case2: ICP+ECP','Location', 'Northeast');
%legend('Case1: SWRO','','Case2: SWRO with ERD','','Case3: SWRO-PRO hybrid','Location', 'Northeast');
ylabel('FW [m^3/h]','FontSize',16);xlabel('SEC_{net} [kWh/m^3]','FontSize',16);

%% initialize the variables + DATA sets
%startTime=datetime("now");
%load("Output_DATA/DATA_Case_1.mat");clc
%X1a_pareto=X_init;X1b_pareto=X_init;X1c_pareto=X_init;
%Y1a_pareto=zeros(18,200);Y1b_pareto=zeros(18,200);Y1c_pareto=zeros(18,200);
%endTime=datetime("now");
%t1a_pareto = endTime - startTime;t1b_pareto=t1a_pareto;t1c_pareto=t1a_pareto;
%save Output_DATA/DATA_Case_1.mat X1a_pareto X1b_pareto X1c_pareto Y1a_pareto Y1b_pareto Y1c_pareto t1a_pareto t1b_pareto t1c_pareto 



























%% Epsilon constraint method 
%load("DATA_case1.mat") 
%startTime=datetime("now");
%A= [-1 1]; b= -0; Aeq=[]; beq=[]; lb = [30;30]; ub = [70;70];
%option_mesh = 1e4; option_BVP = 1e-6; option_data = 1;
%foptions = optimoptions('fmincon','Display','off','Algorithm','interior-point', 'StepTolerance', 1e-12, 'OptimalityTolerance',1e-4, 'MaxFunEvals',100);
%rng default
%FWmin=linspace(0.0065, 0.3861, 200); 
%X0=[linspace(60,69.99,200);linspace(55,50,200)];
%
parfor i=1:0
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
%endTime=datetime("now");     
%time1 = endTime - startTime;
%save DATA_case1.mat X1 Y1 X1_P Y1_P X1_R Y1_R time1 time1_P



%% plot of the Pareto fronts
close all;
load("DATA_case1.mat")
f=figure(1); f.Position = [1200 500 800 500];
scatter(Y1(1,:),Y1(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
scatter(Y1_P(1,:),Y1_P(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
scatter(Y1_R(1),Y1_R(2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','y'); hold on
load("DATA_case2.mat")
scatter(Y2(1,:),Y2(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
scatter(Y2_P(1,:),Y2_P(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
scatter(Y2_R(1),Y2_R(2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','y'); hold on
load("DATA_case3.mat")
scatter3(Y3_P(1,:),Y3_P(2,:),1:1:200,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.9290 0.6940 0.1250]); hold on
scatter(Y3_sqp(1),Y3_sqp(2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','y'); hold on
grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
legend('Case1: \epsilon-constraint method', 'Case1: paretosearch','','Case2: \epsilon-constraint method', 'Case2: paretosearch','','Case3: paretosearch', '','Location', 'Northeast');
%legend('Case1: SWRO','','Case2: SWRO with ERD','','Case3: SWRO-PRO hybrid','Location', 'Northeast');
ylabel('FW [m^3/h]','FontSize',16);xlabel('SEC_{net} [kWh/m^3]','FontSize',16);
