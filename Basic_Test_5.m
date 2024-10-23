%% Simulations for Hybrid I pareto front with varying L^{PRO}
% press [ctrl]+[enter] to run code sections
addpath('Input_DATA','Scaled_model','Unscaled_model','Output_DATA')

%% Pareto front plot
clc;close all
clear all
load("Output_DATA/DATA_Case_1.mat");load("Output_DATA/DATA_Case_2.mat");
load Test_4_DATA.mat
% plot
f=figure(1);f.Position=[679 388.3333 1.6927e+03 802.6667];
scatter(Y1_pareto(:,1),Y1_pareto(:,2),'red');hold on; scatter(Y2_pareto(:,1),Y2_pareto(:,2),'blue');hold on;
xlim([-5.5 0]);ylim([0.1 1.45]);grid on;title('Pareto front','FontSize',14);xlabel('SEC_{net} [kWh/m^3]','FontSize',12);ylabel('FW [m^3/h]','FontSize',12);
%
%
lengths=[1,5,7,10,15,20,21,22,23];
for l=1:length(lengths)
    index=find(X_test(5,:)==L(lengths(l)));
    X=X_test(:,index);Y=Y_test(:,index);
    for j=1:length(index)
            if Y(5,j)<35
                Y(:,j)=zeros(28,1);
            end
    end
    X=[X; Y(1:2,:)]; % 7x8800 matrix
    Fmin=min(X(7,:));[Fmax, Fmax_index]=max(X(7,:));F=linspace(Fmin,Fmax,50);X5_init=zeros(5,50);X5_init(:,end)=X(1:5,Fmax_index);Y5_init=zeros(2,50);Y5_init(:,end)=X(6:7,Fmax_index);
    for i= 1:49
        f1=find(X(7,:)>=F(i) & X(7,:)<F(i+1));
        X5=X(:,f1);[a,b]=max(X5(6,:));
        if isempty(b)==1
            X5_init(:,i)=X5_init(:,i-1);
            Y5_init(:,i)=Y5_init(:,i-1);
        else
        X5_init(:,i)=X5(1:5,b);
        Y5_init(:,i)=X5(6:7,b);
        end
    end
    if l==1; X_1=X5_init;Y_1=Y5_init;end
    if l==2; X_2=X5_init;Y_2=Y5_init;end
    if l==3; X_3=X5_init;Y_3=Y5_init;end
    if l==4; X_4=X5_init;Y_4=Y5_init;end
    if l==5; X_5=X5_init;Y_5=Y5_init;end
    if l==6; X_6=X5_init;Y_6=Y5_init;end
    if l==7; X_7=X5_init;Y_7=Y5_init;end
    if l==8; X_8=X5_init;Y_8=Y5_init;end
    if l==9; X_9=X5_init;Y_9=Y5_init;end
end
%
Color=turbo(30);
Color=flipud(Color());
%
scatter(Y_1(1,:),Y_1(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(1,:));hold on
scatter(Y_2(1,:),Y_2(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(4,:));hold on
scatter(Y_3(1,:),Y_3(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(5,:));hold on
scatter(Y_4(1,:),Y_4(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(10,:));hold on
scatter(Y_5(1,:),Y_5(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(15,:));hold on
scatter(Y_6(1,:),Y_6(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(20,:));hold on
scatter(Y_7(1,:),Y_7(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(23,:));hold on
scatter(Y_8(1,:),Y_8(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(27,:));hold on
scatter(Y_9(1,:),Y_9(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(30,:));hold on
%
lgd=legend('SWRO','SWRO+ERD','L^{PRO} = 0.9626m \times 1','L^{PRO} = 0.9626m \times 5','L^{PRO} = 0.9626m \times 7','L^{PRO} = 0.9626m \times 10','L^{PRO} = 0.9626m \times 20','L^{PRO} = 0.9626m \times 30','L^{PRO} = 0.9626m \times 40','L^{PRO} = 0.9626m \times 50','L^{PRO} = 0.9626m \times 60','Location', 'bestoutside');
title(lgd,'Always: L^{RO} = 0.9626m \times 7')
%

%% test PRO recovery rate
clc
k=3; %bis 50 
option_mesh = 1e3; option_BVP = 1e-3; option_data = .6;
a=fun_scaled(X_9(:,k),option_data,'fig',option_mesh,option_BVP);
ev(a,1:5)

%% Berechnung der Pareto front 
%
for i=1:1
    if i==1; X_init=X_1';elseif i==2; X_init=X_2';elseif i==3; X_init=X_3';elseif i==4; X_init=X_4';elseif i==5; X_init=X_5';elseif i==6; X_init=X_6';elseif i==7; X_init=X_7';elseif i==8; X_init=X_8';else; X_init=X_9';end
    L_pro=L(lengths(i));
    disp('_____________________________________________________________')
    disp(['Starting Pareto front calculation for L_pro = ',num2str(L_pro),'.'])
    % out X has 5 elements: Pd0, PdL, pd0, pf0, Lpro
    A= [-1 1 0 0 0; 1 -1 0 0 0]; b= [0; 3.4]; Aeq=[]; beq=[]; lb = [30;30;4.99;1.00000001;L_pro]; ub = [70;70;20;5;L_pro];
    option_mesh = 1e3; option_BVP = 1e-3; option_data = .6;
    % paretosearch
    options = optimoptions('paretosearch','ParetoSetSize',50, 'InitialPoints',X_init,'Display','iter','MaxFunctionEvaluations',10000, 'ParetoSetChangeTolerance',1e-6,'UseParallel', true);
    X = paretosearch(@(x)fun_scaled(x,option_data,'Pareto',option_mesh,option_BVP),5,A,b,Aeq,beq,lb,ub,[],options);
    % evaluate optimal points
    X_pareto=X';Y_pareto=zeros(50,28);
    parfor t=1:50
        Y_pareto(t,:)=fun_scaled(X(t,:),option_data,'sol',option_mesh,option_BVP);
    end
    if i==1; Xp_1=X_pareto; Yp_1=Y_pareto;Xp_2=Xp_1;Xp_3=Xp_1;Xp_4=Xp_1;Xp_5=Xp_1;Xp_6=Xp_1;Xp_7=Xp_1;Xp_8=Xp_1;Xp_9=Xp_1;Yp_9=Yp_1;Yp_8=Yp_1;Yp_7=Yp_1;Yp_6=Yp_1;Yp_5=Yp_1;Yp_4=Yp_1;Yp_3=Yp_1;Yp_2=Yp_1;
        elseif i==2; Xp_2=X_pareto; Yp_2=Y_pareto;
            elseif i==3; Xp_3=X_pareto; Yp_3=Y_pareto;
                elseif i==4; Xp_4=X_pareto; Yp_4=Y_pareto;
                    elseif i==5; Xp_5=X_pareto; Yp_5=Y_pareto;
                        elseif i==6; Xp_6=X_pareto; Yp_6=Y_pareto;
                            elseif i==7; Xp_7=X_pareto; Yp_7=Y_pareto;
                                elseif i==8; Xp_8=X_pareto; Yp_8=Y_pareto;
    else; Xp_9=X_pareto; Yp_9=Y_pareto;
    end
save Output_DATA/Hybrid_L.mat Xp_1 Xp_2 Xp_3 Xp_4 Xp_5 Xp_6 Xp_7 Xp_8 Xp_9 Yp_1 Yp_2 Yp_3 Yp_4 Yp_5 Yp_6 Yp_7 Yp_8 Yp_9
end
system('git add .'); system('git commit -m "Hybrid-paretosearch"');system('git push https://github.com/oliver-mx/GitMATLAB.git');































%% plot for different PRO lengths
clc;close all
load("Output_DATA/DATA_Case_1.mat");load("Output_DATA/DATA_Case_2.mat");
load Test_4_DATA.mat
% plot
f=figure(1);f.Position=[1.0183e+03 349.6667 1.1887e+03 888.0000];
% farben:
%Color=summer(23);
%Color=jet(30);
Color=turbo(30);
Color=flipud(Color());
% remove some PRO lengths from the plot
%R=[];r=0;
R=[2,3,4,6,8,9,11,12,13,14,16,17,18,19];r=0;
%R=[1:22];r=0;
% scatter plots
for i=1:23
    if any(R==i)
        r=1;
    else
        index=find(X_test(5,:)==L(i));
        farbe=Color(i,:);
        %if i==20;farbe='yellow';end
        if i==23;farbe='black';end
        for j=1:length(index)
            if Y_test(5,index(j))<35
                Y_test(:,index(j))=zeros(28,1);
            end
        end
        scatter(Y_test(1,index),Y_test(2,index),'MarkerEdgeColor','none','MarkerFaceColor',farbe);hold on
    end
end
% plot SWRO and SWRO+ERD
scatter(Y1_pareto(:,1),Y1_pareto(:,2),'red','filled');hold on; 
scatter(Y2_pareto(:,1),Y2_pareto(:,2),'blue','filled');hold on;
xlim([-5.5 0]);ylim([0.1 1.45]);
grid on;title('Pareto front','FontSize',14);
xlabel('SEC_{net} [kWh/m^3]','FontSize',12);ylabel('FW [m^3/h]','FontSize',12);
% create legend
if r==0
lgd=legend('L^{PRO} = 0.9626m \times 1',' ',' ',' ',' ',' ',' ',' ',' ','L^{PRO} = 0.9626m \times 10',' ',' ',' ',' ','L^{PRO} = 0.9626m \times 20',' ',' ',' ',' ','L^{PRO} = 0.9626m \times 30','L^{PRO} = 0.9626m \times 40','L^{PRO} = 0.9626m \times 50','L^{PRO} = 0.9626m \times 60','Case1: SWRO','Case2: SWRO+ERD','Location', 'bestoutside');
title(lgd,'Legend')
end
if sum(R)==253
lgd=legend('Hybrid: L^{RO} \approx 6.7m,  L^{PRO} \approx 57.7m','SWRO: L^{RO} \approx 6.7m','SWRO+ERD L^{RO} \approx 6.7m','Location', 'bestoutside');
title(lgd,'Legend')
end
figure(1); set(gcf,'color','w'); f = gcf; exportgraphics(f,'PRO_length_test.png');



















