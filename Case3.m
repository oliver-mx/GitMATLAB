%% Pareto Front for SWRO+ERD system
% press [ctrl]+[enter] to run code sections
addpath('Input_DATA','Scaled_model','Unscaled_model','Output_DATA')

%% gather Initial data
clc;close all
load("Test_3_DATA.mat");
% remove data with REC_PRO < 70%
for i=1:n
    if Y_hyI(5,i) < 70
        Y_hyI(:,i)=NaN(28,1);
    end
end
% Hybrid I
X=[X_hyI' Y_hyI(1:2,:)'];
Fmin=min(X(:,6));[Fmax, Fmax_index]=max(X(:,6));F=linspace(Fmin,Fmax,200);X3_init=zeros(4,200);X3_init(:,end)=X(Fmax_index,1:4);Y3_init=zeros(2,200);Y3_init(:,end)=X(Fmax_index,5:6);
for i= 1:199
    f1=find(X(:,6)>=F(i) & X(:,6)<F(i+1));
    X1=X(f1,:);
    [a,b]=max(X1(:,5));
    if b>=1
        X3_init(:,i)=X1(b,1:4);
        Y3_init(:,i)=X1(b,5:6);
    else
        X3_init(:,i)=X3_init(:,i-1);
        Y3_init(:,i)=Y3_init(:,i-1); 
    end
end
% plot
f=figure(1);f.Position=[1000 727.6667 1207 510.0000];tiledlayout(1,2);nexttile
scatter(Y1_init(1,:),Y1_init(2,:),'m');hold on; scatter(Y2_init(1,:),Y2_init(2,:),'c');hold on;
scatter(Y3_init(1,:),Y3_init(2,:),'k');hold on;
xlim([-5.5 0]);ylim([0.1 1.45]);grid on;title('Paretosearch - initial data','FontSize',14);xlabel('SEC_{net} [kWh/m^3]','FontSize',12);ylabel('FW [m^3/h]','FontSize',12)
nexttile;
scatter(X1_init(1,:),X1_init(2,:),'m');hold on; scatter(X2_init(1,:),X2_init(2,:),'c');hold on;
xlim([29.4 70.6]);ylim([29.4 70.6]);grid on;title('operating pressures','FontSize',14);xlabel('P_d(0) [bar]','FontSize',12);ylabel('P_d(L) [bar]','FontSize',12)

%% Paretosearch - Case3
clc;close all
load("Output_DATA/DATA_Case_3.mat");clc
disp('________________________________________________________________')
disp('Starting Pareto front calculation for the Hybrid I configuraton.')
startTime=datetime("now"); 
A= [-1 1 0 0; 1 -1 0 0]; b= [0; 3.4]; Aeq=[]; beq=[]; lb = [30;30;5;1]; ub = [70;70;20;5];
option_mesh = 1e4; option_BVP = 1e-4; option_data = 3;
% paretosearch
options = optimoptions('paretosearch','ParetoSetSize',200, 'InitialPoints',X3_init','Display','iter','MaxFunctionEvaluations',10000, 'ParetoSetChangeTolerance',1e-6,'UseParallel', true);
X = paretosearch(@(x)fun_scaled(x,option_data,'Pareto',option_mesh,option_BVP),4,A,b,Aeq,beq,lb,ub,[],options);
% evaluate optimal points
X3_pareto=X';Y3_pareto=zeros(200,28);
parfor i=1:200
    Y3_pareto(i,:)=fun_scaled(X(i,:),option_data,'sol',option_mesh,option_BVP);
end
endTime=datetime("now");
t3_pareto = endTime - startTime;
save Output_DATA/DATA_Case_3.mat X3_init Y3_init X3_pareto Y3_pareto t3_pareto X0 n
disp('Finished calculating the Pareto front for the Hybrid I case.')
disp('____________________________________________________________________')
system('git add .'); system('git commit -m "Case3-paretosearch"');system('git push https://github.com/oliver-mx/GitMATLAB.git');
%

%% Plot the pareto front 
clc;close all
load("Output_DATA/DATA_Case_1.mat");load("Output_DATA/DATA_Case_2.mat");load("Output_DATA/DATA_Case_3.mat");
% plot
f=figure(1);f.Position=[1000 727.6667 1207 510.0000];tiledlayout(1,2);nexttile
scatter(Y1_pareto(:,1),Y1_pareto(:,2),'red');hold on; scatter(Y2_pareto(:,1),Y2_pareto(:,2),'blue');hold on;scatter(Y3_pareto(:,1),Y3_pareto(:,2),'k');hold on;
xlim([-5.5 0]);ylim([0.1 1.45]);grid on;title('Pareto front','FontSize',14);xlabel('SEC_{net} [kWh/m^3]','FontSize',12);ylabel('FW [m^3/h]','FontSize',12);legend('Case1: SWRO','Case2: SWRO+ERD','Case3: Hybrid I','Location', 'best');
nexttile;
scatter(X1_pareto(1,:),X1_pareto(2,:),'red');hold on;scatter(X2_pareto(1,:),X2_pareto(2,:),'blue');hold on;
xlim([29.4 70.6]);ylim([29.4 70.6]);grid on;title('Pareto front - operating pressures','FontSize',14);xlabel('P_d(0) [bar]','FontSize',12);ylabel('P_d(L) [bar]','FontSize',12);legend('Case1: SWRO','Case2: SWRO+ERD','Location', 'best');

