%% Simulations for Hybrid I pareto front
% press [ctrl]+[enter] to run code sections
addpath('Input_DATA','Scaled_model','Unscaled_model','Output_DATA')

%% Pareto front plot
clc;close all
load("Output_DATA/DATA_Case_1.mat");load("Output_DATA/DATA_Case_2.mat");
% plot
f=figure(1);f.Position=[1000 727.6667 1207 510.0000];tiledlayout(1,2);nexttile
scatter(Y1_pareto(:,1),Y1_pareto(:,2),'red');hold on; scatter(Y2_pareto(:,1),Y2_pareto(:,2),'blue');hold on;
xlim([-5.5 0]);ylim([0.1 1.45]);grid on;title('Pareto front','FontSize',14);xlabel('SEC_{net} [kWh/m^3]','FontSize',12);ylabel('FW [m^3/h]','FontSize',12);legend('Case1: SWRO','Case2: SWRO+ERD','Location', 'best');
nexttile;
scatter(X1_pareto(1,:),X1_pareto(2,:),'red');hold on;scatter(X2_pareto(1,:),X2_pareto(2,:),'blue');hold on;
xlim([29.4 70.6]);ylim([29.4 70.6]);grid on;title('Pareto front - operating pressures','FontSize',14);xlabel('P_d(0) [bar]','FontSize',12);ylabel('P_d(L) [bar]','FontSize',12);legend('Case1: SWRO','Case2: SWRO+ERD','Location', 'best');

%% Simulate Hybrid I
load Output_DATA/DATA_Case_2.mat X2_pareto n
% define new input data
X_hyI=[X0; zeros(2,n)];Y_hyI=zeros(28,n);k=0;
P=[5 20];p=[1.3 1.25 1.2 1.15 1.1 1.08 1.05 1.03 1.01 1.008 1.005 1.003 1.001 1.0008 1.0005 1.0003 1.0001 1.00008 1.00005 1.00003 1.00001 1.000008 1.000005 1.000003 1.000001];
for i=1:200
    for j=1:2
        for t=1:25
            k=k+1;
            X_hyI(3,k)=P(j);
            X_hyI(4,k)=p(t);
        end
    end
end
% simulate the hybrid system
option_mesh = 1e4; option_BVP = 1e-4; option_data = 3;
parfor i=1:n
    Y_hyI(:,i)=fun_scaled(X_hyI(:,i),option_data,'sol',option_mesh,option_BVP);
end
save Test_3_DATA.mat X_hyI Y_hyI

%% Plot the results
clc;close all
load("Output_DATA/DATA_Case_1.mat");load("Output_DATA/DATA_Case_2.mat");load("Test_3_DATA.mat");
% remove data with REC_PRO < 70%
for i=1:n
    if Y_hyI(5,i) < 70
        Y_hyI(:,i)=NaN(28,1);
    end
end
% plot
f=figure(1);f.Position=[1000 727.6667 1207 510.0000];tiledlayout(1,2);nexttile
scatter(Y_hyI(1,:),Y_hyI(2,:),'k');hold on;
scatter(Y1_pareto(:,1),Y1_pareto(:,2),'red');hold on; scatter(Y2_pareto(:,1),Y2_pareto(:,2),'blue');hold on;
xlim([-5.5 0]);ylim([0.1 1.45]);grid on;title('Pareto front','FontSize',14);xlabel('SEC_{net} [kWh/m^3]','FontSize',12);ylabel('FW [m^3/h]','FontSize',12);legend('Hybrid I','Case1: SWRO','Case2: SWRO+ERD','Location', 'best');
nexttile;
scatter(X1_pareto(1,:),X1_pareto(2,:),'red');hold on;scatter(X2_pareto(1,:),X2_pareto(2,:),'blue');hold on;
xlim([29.4 70.6]);ylim([29.4 70.6]);grid on;title('Pareto front - operating pressures','FontSize',14);xlabel('P_d(0) [bar]','FontSize',12);ylabel('P_d(L) [bar]','FontSize',12);legend('Case1: SWRO','Case2: SWRO+ERD','Location', 'best');

% create skript Case3 and start paretosearch




























