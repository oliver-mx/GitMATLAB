%% Pareto Front for the simple SWRO system
% press [ctrl]+[enter] to run code sections
addpath('Input_DATA','Scaled_model','Unscaled_model','Output_DATA')

%% test random data
X0=zeros(2,10000);
for i=1:10000
    X0(1,i)= 31+39*rand();
    X0(2,i)= X0(1,i) - .1 - 3.3*rand();
end
Y1=zeros(10000,18);Y2=zeros(10000,18);Y3=zeros(10000,18);Y4=zeros(10000,18);Y5=zeros(10000,18);Y6=zeros(10000,18);
save TEST_if_it_works X0 Y1 Y2 Y3 Y4 Y5 Y6
%%
parfor i=1:10000
    Y1(i,:)=fun_scaled(X0(:,i),-.1,'sol',1e4,1e-3);
    Y2(i,:)=fun_scaled(X0(:,i),-.2,'sol',1e4,1e-3);
    Y3(i,:)=fun_scaled(X0(:,i),-.3,'sol',1e4,1e-3);
    Y4(i,:)=fun_scaled(X0(:,i),-.4,'sol',1e4,1e-3);
    Y5(i,:)=fun_scaled(X0(:,i),-.5,'sol',1e4,1e-3);
    Y6(i,:)=fun_scaled(X0(:,i),-.6,'sol',1e4,1e-3);
end
save TEST_if_it_works X0 Y1 Y2 Y3 Y4 Y5 Y6

%% evaluate the tests 
load("TEST_if_it_works.mat")
close all; clc
% see what we can display:
ev(Y1(1,:))
userNumber = input('Please enter a number of the set {1,2,3,...,18}: ');
% create meshgrid for X0
xq = linspace(min(X0(1,:)), max(X0(1,:)), 100); % Create a query grid for x
yq = linspace(min(X0(2,:)), max(X0(2,:)), 100); % Create a query grid for y
[X, Y] = meshgrid(xq, yq);
% Optional: Visualize the scattered points and the meshgrid
f=figure(1);
f.Position= [117.6667 334.3333 1.4200e+03 999.3334];
tiledlayout(2,3);
% same colorbar
all_vectors = [Y1(:,userNumber); Y2(:,userNumber); Y3(:,userNumber); Y4(:,userNumber); Y5(:,userNumber); Y6(:,userNumber)];
cmap = jet(256); % Use the 'jet' colormap with 256 colors
cmin = min(all_vectors(:));
cmax = max(all_vectors(:));
%
nexttile
scatter3(X0(1,:), 10*(X0(1,:)-X0(2,:)),Y1(:,userNumber),1.5,Y1(:,userNumber),'filled');hold on; title('simple SWRO (ideal)'); xlabel('Seawater inlet pressure'); ylabel('Pressure drop'); axis equal;view(2); grid on;colormap(cmap);colorbar
xlim([30.5 70]);ylim([1 33]);yticks(10:10:30);
nexttile
scatter3(X0(1,:), 10*(X0(1,:)-X0(2,:)),Y2(:,userNumber),1.5,Y2(:,userNumber),'filled');hold on; title('simple SWRO (with ICP)'); xlabel('Seawater inlet pressure'); ylabel('Pressure drop'); axis equal;view(2); grid on;colormap(cmap);colorbar
xlim([30.5 70]);ylim([1 33]);yticks(10:10:30);
nexttile
scatter3(X0(1,:), 10*(X0(1,:)-X0(2,:)),Y3(:,userNumber),1.5,Y3(:,userNumber),'filled');hold on; title('simple SWRO (ICP+ECP)'); xlabel('Seawater inlet pressure'); ylabel('Pressure drop'); axis equal;view(2); grid on;colormap(cmap);colorbar
xlim([30.5 70]);ylim([1 33]);yticks(10:10:30);
nexttile
scatter3(X0(1,:), 10*(X0(1,:)-X0(2,:)),Y4(:,userNumber),1,Y4(:,userNumber),'filled');hold on; title('SWRO+ERD (ideal)'); xlabel('Seawater inlet pressure'); ylabel('Pressure drop'); axis equal;view(2); grid on;colormap(cmap);colorbar
xlim([30.5 70]);ylim([1 33]);yticks(10:10:30);
nexttile
scatter3(X0(1,:), 10*(X0(1,:)-X0(2,:)),Y5(:,userNumber),1,Y5(:,userNumber),'filled');hold on; title('SWRO+ERD (with ICP)'); xlabel('Seawater inlet pressure'); ylabel('Pressure drop'); axis equal;view(2); grid on;colormap(cmap);colorbar
xlim([30.5 70]);ylim([1 33]);yticks(10:10:30);
nexttile
scatter3(X0(1,:), 10*(X0(1,:)-X0(2,:)),Y6(:,userNumber),1,Y6(:,userNumber),'filled');hold on; title('SWRO+ERD (ICP+ECP)'); xlabel('Seawater inlet pressure'); ylabel('Pressure drop'); axis equal;view(2); grid on;colormap(cmap);colorbar
xlim([30.5 70]);ylim([1 33]);yticks(10:10:30);



























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

%% Pareto Front calculation -- paretosearch
load("Output_DATA/DATA_Case_1.mat"); 
rng default; 
startTime=datetime("now");
A= [-1 1; 1 -1]; b= [0; 3.4]; Aeq=[]; beq=[]; lb = [30;30]; ub = [70;70];
option_mesh = 1e4; option_BVP = 1e-3; option_data = 1;
% initial data
X_init = [linspace(60.5e5,69.5e5,200);linspace(58.5e5,68.5e5,200)]';
% paretosearch
options = optimoptions('paretosearch','ParetoSetSize',200, 'InitialPoints',X_init,'Display','iter', 'MaxFunctionEvaluations',10000, 'ParetoSetChangeTolerance',1e-5,'UseParallel', true);
X = paretosearch(@(x)fun_scaled(x,option_data,'Pareto',option_mesh,option_BVP),2,A,b,Aeq,beq,lb,ub,@(x)nonlcon(x,'default'),options);
% evaluate optimal points
parfor i=1:200
    X1_pareto(:,i)=X(i,:)
    Y1_pareto(:,i)=fun_1([X1_pareto(:,i)],option_data,'sol',option_mesh,option_BVP);
end
endTime=datetime("now");       
t1_pareto = endTime - startTime;
save Output_DATA/DATA_Case_1.mat X1_pareto Y1_pareto t1_pareto X1_epsilon Y1_epsilon t1_epsilon X1_rev Y1_rev t1_rev

%%
%load("Output_DATA/DATA_Case_1.mat"); 
scatter(Y1_pareto(1,:),Y1_pareto(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
grid on;view(2);
%xlim([-5.501 -.5]); ylim([0 0.48])
legend('Case1: paretosearch','Location', 'Northeast');
ylabel('FW [m^3/h]','FontSize',16);xlabel('SEC_{net} [kWh/m^3]','FontSize',16);









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



