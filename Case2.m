%% Pareto Front for SWRO+ERD system
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

%% Paretosearch - Case2
clc;close all
load("Output_DATA/DATA_Case_2.mat");clc
disp('_____________________________________________________________')
disp('Starting Pareto front calculation for the SWRO unit with ERD.')
startTime=datetime("now"); 
A= [-1 1; 1 -1]; b= [0; 3.4]; Aeq=[]; beq=[]; lb = [30;30]; ub = [70;70];
option_mesh = 1e4; option_BVP = 1e-4; option_data = 2;
% paretosearch
options = optimoptions('paretosearch','ParetoSetSize',200, 'InitialPoints',X2_init','Display','iter','MaxFunctionEvaluations',10000, 'ParetoSetChangeTolerance',1e-6,'UseParallel', true);
X = paretosearch(@(x)fun_unscaled(x,option_data,'Pareto',option_mesh,option_BVP),2,A,b,Aeq,beq,lb,ub,[],options);
% evaluate optimal points
X2_pareto=X';Y2_pareto=zeros(200,28);
parfor i=1:200
    Y2_pareto(i,:)=fun_unscaled(X(i,:),option_data,'sol',option_mesh,option_BVP);
end
endTime=datetime("now");
t2_pareto = endTime - startTime;
save Output_DATA/DATA_Case_2.mat X2_init Y2_init X2_pareto Y2_pareto t2_pareto X0 n
disp('Finished calculating the Pareto front for the SWRO+ERD case.')
disp('_________________________________________________________________')
system('git add .'); system('git commit -m "Case2-paretosearch"');system('git push https://github.com/oliver-mx/GitMATLAB.git');
%

%% Plot the pareto front 
clc;close all
load("Output_DATA/DATA_Case_1.mat");load("Output_DATA/DATA_Case_2.mat");
% plot
f=figure(1);f.Position=[1000 727.6667 1207 510.0000];tiledlayout(1,2);nexttile
scatter(Y1_pareto(:,1),Y1_pareto(:,2),'red');hold on; scatter(Y2_pareto(:,1),Y2_pareto(:,2),'blue');hold on;
xlim([-5.5 0]);ylim([0.1 1.45]);grid on;title('Pareto front','FontSize',14);xlabel('SEC_{net} [kWh/m^3]','FontSize',12);ylabel('FW [m^3/h]','FontSize',12);legend('Case1: SWRO','Case2: SWRO+ERD','Location', 'best');
nexttile;
scatter(X1_pareto(1,:),X1_pareto(2,:),'red');hold on;scatter(X2_pareto(1,:),X2_pareto(2,:),'blue');hold on;
xlim([29.4 70.6]);ylim([29.4 70.6]);grid on;title('Pareto front - operating pressures','FontSize',14);xlabel('P_d(0) [bar]','FontSize',12);ylabel('P_d(L) [bar]','FontSize',12);legend('Case1: SWRO','Case2: SWRO+ERD','Location', 'best');

