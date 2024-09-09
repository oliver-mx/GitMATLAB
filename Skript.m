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
scatter3(Y3_P(1,:),Y3_P(2,:), 1:1:200,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.9290 0.6940 0.1250]); hold on
scatter(-2.1330,0.3588,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','y'); hold on
grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
legend('Case1: \epsilon-constraint method', 'Case1: paretosearch','','Case2: \epsilon-constraint method', 'Case2: paretosearch','','Case3: paretosearch', '','Location', 'Northeast');
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

system('git status');
system('git add .');
system('git commit -m "Refined 2.0179m Pareto front"');
system('git push https://github.com/oliver-mx/GitMATLAB.git');
system('shutdown /s /t 30');




































%% detailled simulations
%
% 
% 4 26 166 75 181 138 171 197 77 102 62 113 27 1 5 109 7 124 31 158 3 157
% 73 68 184 125 72 74 177 64 8 2 69 126 128 129 6 34
%
clc
load("DATA_case3.mat")
rng default
Z1=[4 26 166 75 181 138 171 197 77 102 62 113 27 1 5 109 7 124 158 3 157 68 184 125 177 64 8 69 129 6 34];
Z2=[31 73 72 74 2 126 128]; % co-current
%
X_test=zeros(5,38*200)'; % 7600
X_test(:,1)=X3_sqp(1).*ones(1,38*200);
%
k=0;
%
for i=1:length(Z1)
    for t=1:200
        k=k+1;
        A=X3_P(:,Z1(i));
        X_test(k,2) = min(A(2)-.2+.4*rand(),70);
        X_test(k,4) = min(A(4)-.2+.4*rand(),20);
        X_test(k,3) = min(A(3)-.1+.2*rand(),X_test(k,4)-.1);
        X_test(k,5) = max(A(5)-.3+.6*rand(),1.01);
    end
end
for j=1:length(Z2)
    for t2=1:200
        k=k+1;
        A=X3_P(:,Z2(j));
        X_test(k,2) = min(A(2)-.2+.4*rand(),70);
        X_test(k,3) = min(A(3)-.2+.4*rand(),20);
        X_test(k,4) = min(A(4)-.1+.2*rand(),X_test(k,3)-.1);
        X_test(k,5) = max(A(5)-.3+.6*rand(),1.01);
    end
end

Y_test=zeros(5,k);
parfor i=1:7600
    Y_test(:,i)=fun_1(X_test(i,:),3,'sol',1e4,1e-6);
end
%% plots zum vergleich
scatter3(Y_test(1,:),Y_test(2,:),1:1:7600,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
%
grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
beep 

%system('git status');
%system('git add .');
%system('git commit -m "refined simulations"');
%system('git push https://github.com/oliver-mx/GitMATLAB.git');




% 59 159 73 102 204 197 47 237 272 372 620 178 479 520 426 187 770 546 849
% 847 1328 499 548 1393 1232 940 1675 1861 1776 1907 2071 2122 2103 2164
% 2285 2144 2255 2364 2363 2386 2398 2543 2478 2774 2740 2996 2894 2922
% 2904 2824 2818 3274 3328 3346 6322 6239 6314 6259 6328 6311 3721 6329
% 3644 3622 3630 6415 6568 4220 4205 4330 4272 4364 4706 4773 4334 4336 
% 4692 4757 4643 6672 6686 6712 6788 6649 6858 6958 6976 6978 6951 5003
% 5017 7041 7051 7141 7044 7037 7088 7059 7058 7014  7159 7083 7132 7104
% 7270 7467 7298 7576 7436 7575 7496 5789 5651 6004 6093 

% 1 142 107 121 162 93 67  154 52 146 28 43 42 35 53 30 14 13 196 29 200 82


%% 
    
length(a)


%%
X_neu=zeros(200,5);
a=[1 142 107 121 162 93 67 21 45 53 154 52 146 28 43 42 35 53 30 14 13 196 29 200 82  175 12 50 198];
X_neu(1:50,:)=X3_P(:,1:50)';
for i=1:length(a)
    X_neu(50+i,:)=X3_P(:,a(i));
end
b=[59 37 16 18 17 46 59 159 73 102 204 197 47 237 272 372 620 178 479 520 426 187 770 546 849  847 1328 499 548 1393 1232 940 1675 1861 1776 1907 2071 2122 2103 2164  2285 2144 2255 2364 2363 2386 2398 2543 2478 2774 2740 2996 2894 2922  2904 2824 2818 3274 3328 3346 6322 6239 6314 6259 6328 6311 3721 6329  3644 3622 3630 6415 6568 4220 4205 4330 4272 4364 4706 4773 4334 4336   4692 4757 4643 6672 6686 6712 6788 6649 6858 6958 6976 6978 6951 5003  5017 7041 7051 7141 7044 7037 7088 7059 7058 7014  7159 7083 7132 7104  7270 7467 7298 7576 7436 7575 7496 5789 5651 6004 6093];
for i=1:length(b)
    X_neu(50+length(a)+i,:)=X_test(b(i),:);
end
X3_P_initial=X_neu;

%% delete too low flow rates from the plot
load("DATA_case3.mat")
for i=1:200
    if Y3_P(1,i)>-0.69
        Y3_P(:,i)=Y3_P(:,1);
        X3_P(:,i)=X3_P(:,1);
    end
end
scatter3(Y3_P(1,:),Y3_P(2,:),1:1:200,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.9290 0.6940 0.1250]); hold on

%save DATA_case3 X3_sqp Y3_sqp  X3_P Y3_P time3_sqp time3_P X3_sqp_initial X3_P_initial

% 4 26 181 62 27 1 5 7 31 179
N=[4 26 181 62 27 1 5 7 31 179];
for j=1:length(N)
disp(X3_P(:,j)')
end

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




%% test locally around best point so far
%
clc
load DATA_Rev_test.mat
n=20000; %16000 = 5h at home
a1=      1.8.*ones(1,n) + 1.4.*rand(1,n);
a2=      70.*ones(1,n);
a4=      20.*ones(1,n);
a3=      12.*ones(1,n) +  6.*rand(1,n);
a5=   1.001.*ones(1,n) + 5.*rand(1,n);
%
X_Rt=[a1; a2; a3; a4; a5];
%
parfor i=1:n
Y_Rt(:,i)=fun_1(X_Rt(:,i), 3 ,'sol', 1e4 , 1e-6);
end
R=Y_Rt(3,:);
R = R(~isnan(R))';
R = R(~isinf(R))';
a=max(R);
c=find(Y_Rt(3,:)== a);
%
disp(X_Rt(:,c))
disp(Y_Rt(:,c)')
save DATA_Rev_test.mat X_Rt Y_Rt

% we try to beat:
%    2.0179
%   70.0000
%   15.9461
%   20.0000
%    1.7051

%   -2.1330    0.3588    0.4880   43.0727   92.1842

%system('git status');
%system('git add .');
%system('git commit -m "testing again 1.8 - 2.2m haha"');
%system('git push https://github.com/oliver-mx/GitMATLAB.git');
%system('shutdown /s /t 30');

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

%    2.0367
%   70.0000
%   15.6380
%   20.0000
%    1.8386
%   -2.2019    0.3633    0.4866   42.0110   84.4000









%% Simulate hybrid data - for paretosearch initial guess
n=3000;
a1=       50*ones(1,n)  +   13.*rand(1,n);
a3=      12.*ones(1,n)  +   2.*rand(1,n);
a2=  a3 - 1.*ones(1,n)  -  .6.*rand(1,n);
a4=     1.7.*ones(1,n)  +      rand(1,n);
%
X_test=[a1; a2; a3; a4];
%
parfor i=1:n
i,
Y_test(:,i)=fun_1(X_test(:,i), 3 ,'sol', 1e4 , 1e-6);
end
%
scatter3(Y_test(1,:),Y_test(2,:), 1:1:3000,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','g'); hold on
%

load("DATA_case1.mat")
f=figure(1); f.Position = [1200 500 800 500];
scatter(Y1(1,:),Y1(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
scatter(Y1_R(1),Y1_R(2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','y'); hold on
load("DATA_case2.mat")
scatter(Y2(1,:),Y2(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
scatter(Y2_R(1),Y2_R(2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','y'); hold on
load("DATA_case3.mat")
scatter(Y3_R(1),Y3_R(2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','y'); hold on
view(2);
grid on










































%% in the end comparison with "old" besd curve
load DATA_Pareto.mat
scatter(SEC_Pareto_1,FW_Pareto_1);hold on
scatter(SEC_Pareto_2,FW_Pareto_2);hold on
load DATA_Approx.mat
scatter(SEC_approx, FW_approx);hold on
Rev= 2.5.*FW_approx -0.5.*FW_approx.*SEC_approx;
[a,b]=max(Rev),













%% test locally around best point so far
%
clc
load("DATA_case3.mat")
load DATA_testy.mat
rng default
n=30000;
a1=      X3_sqp(1).*ones(1,n);
a2=      58.*ones(1,n) +  12.*rand(1,n);
a2(1:10000)=70.*ones(1,10000);
a4=      17.*ones(1,n) +  3.*rand(1,n);
a4(1+5000:10000+5000)=20.*ones(1,10000);
a3=  a4-    .9.*ones(1,n) +  4.1.*rand(1,n);
a5=   1.001.*ones(1,n) + 3.9.*rand(1,n);
%
X_testy=[a1; a2; a3; a4; a5];
%
parfor i=1:n
Y_testy(:,i)=fun_1(X_testy(:,i), 3 ,'sol', 1e4 , 1e-6);
end

save DATA_testy.mat X_testy Y_testy

%system('git status');
%system('git add .');
%system('git commit -m "Simulations for optimal Ratio"');
%system('git push https://github.com/oliver-mx/GitMATLAB.git');
%system('shutdown /s /t 30');

% optimal revenue
%    2.0179
%   70.0000
%   15.9461
%   20.0000
%    1.7051

%   -2.1330    0.3588    0.4880   43.0727   92.1842
%%
load DATA_testy.mat
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
scatter3(Y_testy(1,:),Y_testy(2,:),1:1:30000,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.9290 0.6940 0.1250]); hold on
scatter(Y3_sqp(1),Y3_sqp(2),'y', 'filled'); hold on
grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
legend('Case1: \epsilon-constraint method', 'Case1: paretosearch','','Case2: \epsilon-constraint method', 'Case2: paretosearch','','Case3: paretosearch', '','Location', 'Northeast');
ylabel('FW [m^3/h]','FontSize',16);xlabel('SEC_{net} [kWh/m^3]','FontSize',16);
%
load DATA_20.mat X_20 Y_20
scatter3(Y_20(1,:),Y_20(2,:),1:1:200,'MarkerEdgeColor','none','MarkerFaceColor','k'); hold on
grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);

%scatter(Y3_P(1,:),Y3_P(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.8 0.5 0.1]); hold on

%%
%
% optimal revenue 
% X3_sqp
%
% von schwarz:
% 2  28 73 98 190 161 25
% 19 61 89 70 36 167 186 84 152 95 160 104 121 18 132 137 62 119 64 176 44
% 138 131 189 145 133 65 155 57 60 184 78 171 158 93 105 41
%
% von orange:
% 29467 11438 28573 18320 14054 10298 14091 23681 11936 29277 13727 10984
% 10466 6693 13360 14893 24984 9352 6966 14226 12681 8128 6474 5656 8028
% 9101 8266 7385 7178 8815 5112 12289 12206 5040 8821 8139 8463 8143 9140
% 8189 8730 9198 6668
%

a1=[1]
a2=[71 191 91 186 22 94 112 82 61 153 142 95 59 120 30 75 123 45 64 21 152 20 104 133 55 63 66 157 169 154 74 180 84 174 144 192 128 193 118 189 184 165 100 58 161 68 62 67 176 135 90];
a4=[ 11873 18320 14054 14091 29736 18910 14472 28418 11303 25103 12020 27472 11500 526 3841 14261 29467 11438 28573 18320 14054 10298 14091 23681 11936 29277 13727 10984  10466 6693 13360 14893 24984 9352 6966 14226 12681 8128 6474 5656 8028  9101 8266 7385 7178 8815 5112 12289]

X3_P_initial = [X3_sqp'];
for i=1:length(a2)
    X3_P_initial = [X3_P_initial; X_20(:,a2(i))'];
end

for j=1:length(a4)
    X3_P_initial = [X3_P_initial; X_testy(:,a4(j))'];
end
X3_P_initial=[X3_P_initial;X3_P_initial];
size(X3_P_initial)
