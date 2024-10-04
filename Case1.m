%% Pareto Front for the simple SWRO system
% press [ctrl]+[enter] to run code sections
addpath('Input_DATA','Scaled_model','Unscaled_model','Output_DATA')

%% test inital data
% random initial data
rng default
X_init = 31.*ones(1,200) + 39.*rand(2,200);
X_init(2,:) = X_init(1,:) - .1.*ones(1,200) - 3.3.*rand(1,200);
X_init(2,:) = max(X_init(2,:), 30*ones(1,200));
% test if initial data works
Y=zeros(200,2);
parfor i=1:200
    Y(i,:)=fun_scaled(X_init(:,i),-.6,'Pareto',1e4,1e-3);
end

%% Paretosearch - ideal SWRO (no ERD)
load("Output_DATA/DATA_Case_1.mat");clc
disp('Starting Pareto front calculation for the single SWRO unit.')
%
startTime=datetime("now");
% initial data
rng default; 
X_init = 31.*ones(1,200) + 39.*rand(2,200);
X_init(2,:) = X_init(1,:) - .1.*ones(1,200) - 3.3.*rand(1,200);
X_init(2,:) = max(X_init(2,:), 30*ones(1,200));
% constraints
A= [-1 1; 1 -1]; b= [0; 3.4]; Aeq=[]; beq=[]; lb = [30;30]; ub = [70;70];
option_mesh = 1e4; option_BVP = 1e-3; option_data = -.1; % <---IMPORTANT
% paretosearch
options = optimoptions('paretosearch','ParetoSetSize',200, 'InitialPoints',X_init','Display','iter','MaxFunctionEvaluations',10000, 'ParetoSetChangeTolerance',1e-5,'UseParallel', true);
X = paretosearch(@(x)fun_scaled(x,option_data,'Pareto',option_mesh,option_BVP),2,A,b,Aeq,beq,lb,ub,@(x)nonlcon(x,'default'),options);
% evaluate optimal points
X1a_pareto=X';Y1a_pareto=zeros(200,18);
parfor i=1:200
    Y1a_pareto(:,i)=fun_1(X(i,:),option_data,'sol',option_mesh,option_BVP);
end
endTime=datetime("now");
t1a_pareto = endTime - startTime;
save Output_DATA/DATA_Case_1.mat X1a_pareto X1b_pareto X1c_pareto Y1a_pareto Y1b_pareto Y1c_pareto t1a_pareto t1b_pareto t1c_pareto 
% push to github
system('git add .');
system('git commit -m "push1"');
system('git push https://github.com/oliver-mx/GitMATLAB.git');
% Paretosearch - ICP SWRO (no ERD)
load("Output_DATA/DATA_Case_1.mat");clc
disp('Starting Pareto front calculation for the single SWRO unit.')
%
startTime=datetime("now");
% initial data
rng default; 
X_init = 31.*ones(1,200) + 39.*rand(2,200);
X_init(2,:) = X_init(1,:) - .1.*ones(1,200) - 3.3.*rand(1,200);
X_init(2,:) = max(X_init(2,:), 30*ones(1,200));
% constraints
A= [-1 1; 1 -1]; b= [0; 3.4]; Aeq=[]; beq=[]; lb = [30;30]; ub = [70;70];
option_mesh = 1e4; option_BVP = 1e-3; option_data = -.2; % <---IMPORTANT
% paretosearch
options = optimoptions('paretosearch','ParetoSetSize',200, 'InitialPoints',X_init','Display','iter','MaxFunctionEvaluations',10000, 'ParetoSetChangeTolerance',1e-5,'UseParallel', true);
X = paretosearch(@(x)fun_scaled(x,option_data,'Pareto',option_mesh,option_BVP),2,A,b,Aeq,beq,lb,ub,@(x)nonlcon(x,'default'),options);
% evaluate optimal points
X1b_pareto=X';Y1b_pareto=zeros(200,18);
parfor i=1:200
    Y1b_pareto(:,i)=fun_1(X(i,:),option_data,'sol',option_mesh,option_BVP);
end
endTime=datetime("now");
t1b_pareto = endTime - startTime;
save Output_DATA/DATA_Case_1.mat X1a_pareto X1b_pareto X1c_pareto Y1a_pareto Y1b_pareto Y1c_pareto t1a_pareto t1b_pareto t1c_pareto 
% push to github
system('git add .');
system('git commit -m "push2"');
system('git push https://github.com/oliver-mx/GitMATLAB.git');
% Paretosearch - ICP+ECP SWRO (no ERD)
load("Output_DATA/DATA_Case_1.mat");clc
disp('Starting Pareto front calculation for the single SWRO unit.')
%
startTime=datetime("now");
% initial data
rng default; 
X_init = 31.*ones(1,200) + 39.*rand(2,200);
X_init(2,:) = X_init(1,:) - .1.*ones(1,200) - 3.3.*rand(1,200);
X_init(2,:) = max(X_init(2,:), 30*ones(1,200));
% constraints
A= [-1 1; 1 -1]; b= [0; 3.4]; Aeq=[]; beq=[]; lb = [30;30]; ub = [70;70];
option_mesh = 1e4; option_BVP = 1e-3; option_data = -.3; % <---IMPORTANT
% paretosearch
options = optimoptions('paretosearch','ParetoSetSize',200, 'InitialPoints',X_init','Display','iter','MaxFunctionEvaluations',10000, 'ParetoSetChangeTolerance',1e-5,'UseParallel', true);
X = paretosearch(@(x)fun_scaled(x,option_data,'Pareto',option_mesh,option_BVP),2,A,b,Aeq,beq,lb,ub,@(x)nonlcon(x,'default'),options);
% evaluate optimal points
X1c_pareto=X';Y1c_pareto=zeros(200,18);
parfor i=1:200
    Y1c_pareto(:,i)=fun_1(X(i,:),option_data,'sol',option_mesh,option_BVP);
end
endTime=datetime("now");
t1c_pareto = endTime - startTime;
save Output_DATA/DATA_Case_1.mat X1a_pareto X1b_pareto X1c_pareto Y1a_pareto Y1b_pareto Y1c_pareto t1a_pareto t1b_pareto t1c_pareto 
% push to github
system('git add .');
system('git commit -m "push3"');
system('git push https://github.com/oliver-mx/GitMATLAB.git');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Paretosearch - ideal SWRO with ERD
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
option_mesh = 1e4; option_BVP = 1e-3; option_data = -.4; % <---IMPORTANT
% paretosearch
options = optimoptions('paretosearch','ParetoSetSize',200, 'InitialPoints',X_init','Display','iter','MaxFunctionEvaluations',10000, 'ParetoSetChangeTolerance',1e-5,'UseParallel', true);
X = paretosearch(@(x)fun_scaled(x,option_data,'Pareto',option_mesh,option_BVP),2,A,b,Aeq,beq,lb,ub,@(x)nonlcon(x,'default'),options);
% evaluate optimal points
X2a_pareto=X';Y2a_pareto=zeros(200,18);
parfor i=1:200
    Y2a_pareto(:,i)=fun_1(X(i,:),option_data,'sol',option_mesh,option_BVP);
end
endTime=datetime("now");
t2a_pareto = endTime - startTime;
save Output_DATA/DATA_Case_2.mat X2a_pareto X2b_pareto X2c_pareto Y2a_pareto Y2b_pareto Y2c_pareto t2a_pareto t2b_pareto t2c_pareto 
% push to github
system('git add .');
system('git commit -m "push4"');
system('git push https://github.com/oliver-mx/GitMATLAB.git');
% Paretosearch - ICP SWRO with ERD
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
option_mesh = 1e4; option_BVP = 1e-3; option_data = -.5; % <---IMPORTANT
% paretosearch
options = optimoptions('paretosearch','ParetoSetSize',200, 'InitialPoints',X_init','Display','iter','MaxFunctionEvaluations',10000, 'ParetoSetChangeTolerance',1e-5,'UseParallel', true);
X = paretosearch(@(x)fun_scaled(x,option_data,'Pareto',option_mesh,option_BVP),2,A,b,Aeq,beq,lb,ub,@(x)nonlcon(x,'default'),options);
% evaluate optimal points
X2b_pareto=X';Y2b_pareto=zeros(200,18);
parfor i=1:200
    Y2b_pareto(:,i)=fun_1(X(i,:),option_data,'sol',option_mesh,option_BVP);
end
endTime=datetime("now");
t2b_pareto = endTime - startTime;
save Output_DATA/DATA_Case_2.mat X2a_pareto X2b_pareto X2c_pareto Y2a_pareto Y2b_pareto Y2c_pareto t2a_pareto t2b_pareto t2c_pareto 
% push to github
system('git add .');
system('git commit -m "push5"');
system('git push https://github.com/oliver-mx/GitMATLAB.git');
% Paretosearch - ICP+ECP SWRO with ERD
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
options = optimoptions('paretosearch','ParetoSetSize',200, 'InitialPoints',X_init','Display','iter','MaxFunctionEvaluations',10000, 'ParetoSetChangeTolerance',1e-5,'UseParallel', true);
X = paretosearch(@(x)fun_scaled(x,option_data,'Pareto',option_mesh,option_BVP),2,A,b,Aeq,beq,lb,ub,@(x)nonlcon(x,'default'),options);
% evaluate optimal points
X2c_pareto=X';Y2c_pareto=zeros(200,18);
parfor i=1:200
    Y2c_pareto(:,i)=fun_1(X(i,:),option_data,'sol',option_mesh,option_BVP);
end
endTime=datetime("now");
t2c_pareto = endTime - startTime;
save Output_DATA/DATA_Case_2.mat X2a_pareto X2b_pareto X2c_pareto Y2a_pareto Y2b_pareto Y2c_pareto t2a_pareto t2b_pareto t2c_pareto 

% push to github
system('git add .');
system('git commit -m "push6"');
system('git push https://github.com/oliver-mx/GitMATLAB.git');

% plot the three Pareto curves:
close all;
load("Output_DATA/DATA_Case_1.mat");
f=figure(1); f.Position = [1200 500 800 500];
%
scatter(Y1a_pareto(1,:),Y1a_pareto(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','#880808'); hold on
scatter(Y1b_pareto(1,:),Y1b_pareto(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','#EE4B2B'); hold on
scatter(Y1c_pareto(1,:),Y1c_pareto(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
%
scatter(Y2a_pareto(1,:),Y2a_pareto(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','#00008B'); hold on
scatter(Y2b_pareto(1,:),Y2b_pareto(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','#0096FF'); hold on
scatter(Y2c_pareto(1,:),Y2c_pareto(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
%
grid on; xlim([-15 -.5]); ylim([0 2]);view(2);
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
