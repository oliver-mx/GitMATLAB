%% Matlab skript for the Paper 2024:
%
% - reproduce the simulations
% - reproduce optimisation results
% - reproduce all figures.png (converted to .eps using LibreOfficeDraw)
% - test average computation time
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Numerical Simulations
%
% Create figures form chapter: Numerical simulations
% 
clc
load("DATA_case3.mat")
[output1, output2, time]=fun_1(X3_sqp, 3,'fig', 1e4 , 1e-6),
%
%figure(1); set(gcf,'color','w'); f = gcf; exportgraphics(f,'Figure_6.png');
%figure(2); set(gcf,'color','w'); f = gcf; exportgraphics(f,'Figure_7.png');
%figure(3); set(gcf,'color','w'); f = gcf; exportgraphics(f,'Figure_8.png');
%figure(4); set(gcf,'color','w'); f = gcf; exportgraphics(f,'Figure_9.png');
%

%% ERD test
output1=fun_1([70,50], 1,'fig', 1e4 , 1e-6),

%% Average simulation time:
%
% calculate average computation time of the simulation above
%
load("DATA_case3.mat")
T=zeros(1,100); 
for i=1:length(T)
[~, ~, time] = fun_1(X3_sqp, 3,'sol', 1e4 , 1e-6);
T(i)=time;
end
mean(T) %uni: 6.1215

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Overview of computation times:
%
%                               |  home PC   |  uni PC    |
%                               |  i12900k   |  i13500k   |
%                               |  16 cores  |  14 cores  |
%                               |  3.2 GHz   |  2.5 GHz   |
%                               |  64GB RAM  |  16GB RAM  |                  
%
% case 1 - epsilon constraint   |  02:22:08  |  04:55:45  |
% case 1 - paretosearch         |  00:52:50! |  01:07:44  |
%
% case 2 - epsilon constraint   |  06:39:07  |  18:29:29  |
% case 2 - paretosearch         |  00:17:24! |  00:46:28  |
% 
% case 3 - paretosearch         |  01:08:46! |  02:11:35  |
%

%% plot for paper
close all
revenue_pareto_plot([.8 1.4 2])
%
%figure(1); set(gcf,'color','w'); f = gcf; exportgraphics(f,'Figure_10.png');
%
%% another plot
close all
plot1()
%
%figure(1); set(gcf,'color','w'); f = gcf; exportgraphics(f,'Figure_11.png');
%

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Case1: epsilon-constraint method
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
%save DATA_case1.mat X1 Y1 X1_P Y1_P X1_R Y1_R time1 time1_P

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
%save DATA_case1.mat X1 Y1 X1_P Y1_P X1_R Y1_R time1 time1_P

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
%save DATA_case1.mat X1 Y1 X1_P Y1_P X1_R Y1_R time1 time1_P

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%save DATA_case2.mat X2 Y2 X2_P Y2_P X2_R Y2_R time2 time2_P

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
%save DATA_case2.mat X2 Y2 X2_P Y2_P X2_R Y2_R time2 time2_P

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
%save DATA_case2.mat X2 Y2 X2_P Y2_P X2_R Y2_R time2 time2_P

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Case3: Max Revenue 
%
load("DATA_case3.mat")
startTime=datetime("now");
A= []; b=[]; Aeq=[]; beq=[]; lb = [1.9;30;10;10;1]; ub = [2.1;70;20;20;5];
option_mesh = 1e4; option_BVP = 1e-6; option_data = .3;
%foptions = optimoptions('fmincon','Display','iter','Algorithm','interior-point', 'StepTolerance', 1e-16, 'OptimalityTolerance',1e-4, 'MaxFunEvals',15000);
foptions = optimoptions('fmincon','Display','iter','Algorithm','sqp', 'StepTolerance', 1e-16, 'OptimalityTolerance',1e-4, 'MaxFunEvals',15000);
rng default
%
X3_sqp = fmincon(@(x)fun_1(x,option_data,'Rev',option_mesh,option_BVP),X3_sqp_initial,A,b,Aeq,beq,lb,ub,@(x)nonlcon(x,'default'),foptions);
Y3_sqp = fun_1(X3_sqp,option_data,'sol',option_mesh,option_BVP);
display(X3_sqp)
display(Y3_sqp)
%
endTime=datetime("now");     
time3_sqp = endTime - startTime; 
%save DATA_case3 X3_sqp Y3_sqp  X3_P Y3_P time3_sqp time3_P X3_sqp_initial X3_P_initial

%% Case3: paretosearch
%
load("DATA_case3.mat")
%X3_P_initial=X3_P';
% 
startTime=datetime("now");
A= []; b=[]; Aeq=[]; beq=[]; lb = [X3_sqp(1);30;10;10;1]; ub = [X3_sqp(1);70;20;20;5];
option_mesh = 1e4; option_BVP = 1e-6; option_data = 3;
rng default
options = optimoptions('paretosearch','ParetoSetSize',200, 'InitialPoints',X3_P_initial,'Display','iter', 'MaxFunctionEvaluations',10000, 'ParetoSetChangeTolerance',1e-6,'UseParallel', true);
X = paretosearch(@(x)fun_1(x,option_data,'Pareto',option_mesh,option_BVP),5,A,b,Aeq,beq,lb,ub,@(x)nonlcon(x,'default'),options); 
%
parfor i=1:200
    X3_P(:,i)=X(i,:)
    Y3_P(:,i)=fun_1([X3_P(:,i)],option_data,'sol',option_mesh,option_BVP);
end
endTime=datetime("now");     
time3_P = endTime - startTime;
save DATA_case3 X3_sqp Y3_sqp  X3_P Y3_P time3_sqp time3_P X3_sqp_initial X3_P_initial

%system('git status');
%system('git add .');
%system('git commit -m "orange pareto set improvements"');
%system('git push https://github.com/oliver-mx/GitMATLAB.git');


%% find points
close all;
plot1()
% 

%%
clc;close all;rng default;load("DATA_case3.mat")
%
k=16*20;
X0=X3_P(:,118);
%
X_init=zeros(5,k);
for i=1:k
    X_init(:,i)=[X0(1) X0(2) X0(3)-.2+.4*rand() X0(4) X0(5)-.2+.4*rand()];
end
Y_init=zeros(5,k);
parfor i=1:k
    Y_init(:,i)=fun_1(X_init(:,i),3,'sol',1e4,1e-6);
end
%% plot zum vergleich
plot1()
scatter3(Y_init(1,:),Y_init(2,:),1:1:k,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
%
grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);

%% replace X and Y in X3_p:
load("DATA_case3.mat")
%
index_neu=[267];
%
index_alt=[154];
%
if length(index_neu)==length(index_alt)
    for i = 1:length(index_neu)
        X3_P(:,index_alt(i)) = X_init(:,index_neu(i));
        Y3_P(:,index_alt(i)) = Y_init(:,index_neu(i));
    end
else
    display('lengths are not equal !')
end
save DATA_case3 X3_sqp Y3_sqp  X3_P Y3_P time3_sqp time3_P X3_sqp_initial X3_P_initial
close all;
plot1()









%%
clc;
option_mesh = 1e4; option_BVP = 1e-6; option_data = 3;  
k=16*100;% initial guess
I1 = X3_sqp(1).*ones(1,k);
I2 = 63.*ones(1,k)+7.*rand(1,k);
I2(1:800)=70.*ones(1,800);
I4 = 17.*ones(1,k)+3.*rand(1,k);
I4(601:1000)=20.*ones(1,400);
I3 = I4- .7.*rand(1,k)- 6.*ones(1,k);
I5 = 1.5.*ones(1,k)+ 2*rand(1,k);
%
X_init=[I1;I2;I3;I4;I5];
Y_init=zeros(5,k);
parfor i=1:k
    Y_init(:,i)=fun_1(X_init(:,i),option_data,'sol',option_mesh,option_BVP);
end
% plots zum vergleich
scatter3(Y_init(1,:),Y_init(2,:),1:1:k,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
%
grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
beep 



%%
% remove points with negative REC
find(Y_Rt(5,:)<0)
R=Y_Rt(3,:);
R(3024)=0;R(10606)=0;R(14386)=0;
R = R(~isnan(R))';
R = R(~isinf(R))';
a=max(R);
c=find(Y_Rt(3,:)== a);
%
disp(X_Rt(:,c))
disp(Y_Rt(:,c)')






%% find best candidates






%% imporve fig 10 right
% lengths
I=[1.6 1.7 1.8 1.9 2.1 2.2 2.3 2.4 2.5];
% coressponding start vector:
X0=[1.6000   70.0000   17.6408   20.0000    1.5522;
    1.7000   70.0000   17.1616   20.0000    1.4811;
    1.8000   70.0000   16.8488   20.0000    1.9343;
    1.9000   70.0000   16.5807   20.0000    2.0394;
    2.1000   70.0000   15.6291   20.0000    2.7275;
    2.2000   70.0000   15.1129   20.0000    3.1129;
    2.3000   70.0000   14.5006   20.0000    3.8143;
    2.4000   70.0000   13.6224   20.0000    4.0084;
    2.5000   70.0000   13.6484   20.0000    4.6532]';
%
rng default
option_mesh = 1e4; option_BVP = 1e-6; option_data = 3;  
%
k=14*60;
X_init=zeros(5,9*k);
Y_init=zeros(5,9*k);
z=0;
for j=1:9
    A=X0(:,j)';
    for i=1:k
        z=z+1;
        X_init(:,z)=[A(1) A(2) A(3)-.2+.4*rand() A(4) A(5)-.2+.4*rand()];
    end
end
%%
parfor i2=1:k*9
        Y_init(:,i2)=fun_1(X_init(:,i2),option_data,'sol',option_mesh,option_BVP);
end

%%
system('git status');
system('git add .');
system('git commit -m "nach Urlaub"');
system('git push https://github.com/oliver-mx/GitMATLAB.git');

%%


%%
%%
clc
XX=X_init(:,1681:1:2238);
YY=Y_init(:,1681:1:2238);
%load Data_25.mat
%z=Y_25(3,:);
z=YY(3,:);
REV=z;
REV= REV(~isnan(REV));
a=max(REV);
c=find(z==a);
c=c(1);
%
%disp(X_25(:,c)')
%disp(Y_25(:,c)')
disp(XX(:,c)')
disp(YY(:,c)')


%%
z=1;
REV= Y_init(3,(z-1)*14*80:1:z*14*80);
REV= REV(~isnan(REV))';
a=max(REV);
c=find(Y_test(3,:)== a);
%
disp(X_test(c,:))
disp(Y_test(:,c)')

%%
load DATA_MaxRev.mat
y_rev(12)=.4825; %1.6 = 12
save DATA_MaxRev.mat x_length y_rev

%% plot functions:
close all
revenue_pareto_plot([.8 1.4 2])
function revenue_pareto_plot(I)
    fig = figure(1);
    fig.Position=[592.3333 516.3333 1.4077e+03 483.3333];
    tiledlayout(1, 2);
    %%%%%%%%
    nexttile;
    L = 0.5:0.1:2.5;
    load DATA_05.mat X_05 Y_05
    load DATA_06.mat X_06 Y_06
    load DATA_07.mat X_07 Y_07
    load DATA_08.mat X_08 Y_08
    load DATA_09.mat X_09 Y_09
    load DATA_10.mat X_10 Y_10
    load DATA_11.mat X_11 Y_11
    load DATA_12.mat X_12 Y_12
    load DATA_13.mat X_13 Y_13
    load DATA_14.mat X_14 Y_14
    load DATA_15.mat X_15 Y_15
    load DATA_16.mat X_16 Y_16
    load DATA_17.mat X_17 Y_17
    load DATA_18.mat X_18 Y_18
    load DATA_19.mat X_19 Y_19
    load DATA_20.mat X_20 Y_20
    load DATA_21.mat X_21 Y_21
    load DATA_22.mat X_22 Y_22
    load DATA_23.mat X_23 Y_23
    load DATA_24.mat X_24 Y_24
    load DATA_25.mat X_25 Y_25
    Color=summer(30);
    Color=flipud(Color());
    if any(I==.5);scatter(Y_05(1,:),Y_05(2,:),15,'MarkerEdgeColor','none','MarkerFaceColor',Color(1,:)); hold on; end
    if any(I==.6);scatter(Y_06(1,:),Y_06(2,:),15,'MarkerEdgeColor','none','MarkerFaceColor',Color(2,:)); hold on; end
    if any(I==.7);scatter(Y_07(1,:),Y_07(2,:),15,'MarkerEdgeColor','none','MarkerFaceColor',Color(3,:)); hold on; end
    if any(I==.8);scatter(Y_08(1,:),Y_08(2,:),15,'MarkerEdgeColor','none','MarkerFaceColor',Color(4,:)); hold on; end
    if any(I==.9);scatter(Y_09(1,:),Y_09(2,:),15,'MarkerEdgeColor','none','MarkerFaceColor',Color(5,:)); hold on; end
    if any(I==1);scatter(Y_10(1,:),Y_10(2,:),15,'MarkerEdgeColor','none','MarkerFaceColor',Color(6,:)); hold on; end
    if any(I==1.1);scatter(Y_11(1,:),Y_11(2,:),15,'MarkerEdgeColor','none','MarkerFaceColor',Color(7,:)); hold on; end
    if any(I==1.2);scatter(Y_12(1,:),Y_12(2,:),15,'MarkerEdgeColor','none','MarkerFaceColor',Color(8,:)); hold on; end
    if any(I==1.3);scatter(Y_13(1,:),Y_13(2,:),15,'MarkerEdgeColor','none','MarkerFaceColor',Color(9,:)); hold on; end
    if any(I==1.4);scatter(Y_14(1,:),Y_14(2,:),15,'MarkerEdgeColor','none','MarkerFaceColor',Color(10,:)); hold on; end
    if any(I==1.5);scatter(Y_15(1,:),Y_15(2,:),15,'MarkerEdgeColor','none','MarkerFaceColor',Color(11,:)); hold on; end
    if any(I==1.6);scatter(Y_16(1,:),Y_16(2,:),15,'MarkerEdgeColor','none','MarkerFaceColor',Color(12,:)); hold on; end
    if any(I==1.7);scatter(Y_17(1,:),Y_17(2,:),15,'MarkerEdgeColor','none','MarkerFaceColor',Color(13,:)); hold on; end
    if any(I==1.8);scatter(Y_18(1,:),Y_18(2,:),15,'MarkerEdgeColor','none','MarkerFaceColor',Color(14,:)); hold on; end
    if any(I==1.9);scatter(Y_19(1,:),Y_19(2,:),15,'MarkerEdgeColor','none','MarkerFaceColor',Color(15,:)); hold on; end
    if any(I==2);scatter(Y_20(1,:),Y_20(2,:),15,'MarkerEdgeColor','none','MarkerFaceColor',Color(16,:)); hold on; end
    if any(I==2.1);scatter(Y_21(1,:),Y_21(2,:),15,'MarkerEdgeColor','none','MarkerFaceColor',Color(17,:)); hold on; end
    if any(I==2.2);scatter(Y_22(1,:),Y_22(2,:),15,'MarkerEdgeColor','none','MarkerFaceColor',Color(18,:)); hold on; end
    if any(I==2.3);scatter(Y_23(1,:),Y_23(2,:),15,'MarkerEdgeColor','none','MarkerFaceColor',Color(19,:)); hold on; end
    if any(I==2.4);scatter(Y_24(1,:),Y_24(2,:),15,'MarkerEdgeColor','none','MarkerFaceColor',Color(20,:)); hold on; end
    if any(I==2.5);scatter(Y_25(1,:),Y_25(2,:),15,'MarkerEdgeColor','none','MarkerFaceColor',Color(21,:)); hold on; end 
    grid on;view(2);
    ylabel('FW [m^3/h]','FontSize',16);xlabel('SEC_{net} [kWh/m^3]','FontSize',16);
    %xlim([-5.501 -.5]); ylim([0 0.48])
    if length(I)==1;legend(['L^{RO}= 4m, L^{PRO}= ',num2str(I(1)),'m'],'Location', 'Northeast');end
    if length(I)==2;legend(['L^{RO}= 4m, L^{PRO}= ',num2str(I(1)),'m'],['L^{RO}= 4m, L^{PRO}= ',num2str(I(2)),'m'],'Location', 'SouthWest');end
    if length(I)==3;legend(['L^{RO}= 4m, L^{PRO}= ',num2str(I(1)),'m'],['L^{RO}= 4m, L^{PRO}= ',num2str(I(2)),'m'],['L^{RO}= 4m, L^{PRO}= ',num2str(I(3)),'m'],'Location', 'SouthWest');end
    if length(I)==4;legend(['L^{RO}= 4m, L^{PRO}= ',num2str(I(1)),'m'],['L^{RO}= 4m, L^{PRO}= ',num2str(I(2)),'m'],['L^{RO}= 4m, L^{PRO}= ',num2str(I(3)),'m'],['L^{RO}= 4m, L^{PRO}= ',num2str(I(4)),'m'],'Location', 'SouthWest');end
    if length(I)>4;legend(['L^{RO}= 4m, L^{PRO}= ',num2str(I(1)),'m'],['L^{RO}= 4m, L^{PRO}= ',num2str(I(2)),'m'],['L^{RO}= 4m, L^{PRO}= ',num2str(I(3)),'m'],['L^{RO}= 4m, L^{PRO}= ',num2str(I(4)),'m'],['L^{RO}= 4m, L^{PRO}= ',num2str(I(5)),'m'],'Location', 'SouthWest');end
    %
    %load DATA_case3.mat Y3_P
    %scatter3(Y3_P(1,:),Y3_P(2,:),1:1:200,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.9290 0.6940 0.1250]); hold on
    %
    %%%%%%%%
    nexttile;
    % Create a figure and axis
    Color=summer(30);
    Color=flipud(Color());
    a=zeros(1,21);
    s1=2;s2=.3;
    %
    load DATA_19.mat X_19 Y_19
    load DATA_22.mat X_22 Y_22
    load DATA_23.mat X_23 Y_23
    load DATA_24.mat X_24 Y_24
    load DATA_25.mat X_25 Y_25
    load DATA_case3.mat X3_sqp Y3_sqp
    %
    a(15)=max(s1.*Y_19(2,:) + s2.*Y_19(2,:).*Y_19(1,:));
    scatter(X_19(1,1),a(15),'MarkerEdgeColor','none','MarkerFaceColor',Color(19,:));hold on
    %
    a(18)=max(s1.*Y_22(2,:)+ s2.*Y_22(2,:).*Y_22(1,:));
    scatter(X_22(1,1),a(18),'MarkerEdgeColor','none','MarkerFaceColor',Color(22,:));hold on
    %
    a(19)=max(s1.*Y_23(2,:) + s2.*Y_23(2,:).*Y_23(1,:));
    scatter(X_23(1,1),a(19),'MarkerEdgeColor','none','MarkerFaceColor',Color(23,:));hold on
    %
    a(20)=max(s1.*Y_24(2,:) + s2.*Y_24(2,:).*Y_24(1,:));
    scatter(X_24(1,1),a(20),'MarkerEdgeColor','none','MarkerFaceColor',Color(24,:));hold on
    %
    a(21)=max(s1.*Y_25(2,:) + s2.*Y_25(2,:).*Y_25(1,:));
    scatter(X_25(1,1),a(21),'MarkerEdgeColor','none','MarkerFaceColor',Color(25,:));hold on
    %
    load DATA_MaxRev.mat x_length y_rev
    for i=1:21
        scatter(x_length(i),y_rev(i),'MarkerEdgeColor','none','MarkerFaceColor',Color(4+i,:));hold on
    end
    scatter(X3_sqp(1),Y3_sqp(3),'MarkerEdgeColor','none','MarkerFaceColor',[0.9290 0.6940 0.1250]);hold on
    %
    grid on; xlim([.46 2.54]); ylim([0.4551 .4951]);view(2);
    ylabel('REV [$/h]','FontSize',16);xlabel('L^{PRO} [m]','FontSize',16);
end
function plot1()
    load DATA_case1.mat Y1 Y1_P
    f=figure(1); f.Position = [1200 500 800 500];
    scatter(Y1(1,:),Y1(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
    scatter(Y1_P(1,:),Y1_P(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
    load DATA_case2.mat Y2 Y2_P
    scatter(Y2(1,:),Y2(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
    scatter(Y2_P(1,:),Y2_P(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
    load DATA_case3.mat Y3_P
    scatter3(Y3_P(1,:),Y3_P(2,:),1:1:200,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.9290 0.6940 0.1250]); hold on
    grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
    legend('Case1: SWRO','','Case2: SWRO with ERD','','Case3: SWRO-PRO hybrid','Location', 'Northeast');
    ylabel('FW [m^3/h]','FontSize',16);xlabel('SEC_{net} [kWh/m^3]','FontSize',16);
end