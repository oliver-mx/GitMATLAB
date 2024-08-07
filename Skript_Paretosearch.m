%% plot of the 3 poitnts with maximum revenue:
%
f=figure(1); f.Position = [1200 500 800 500];
scatter(-4,0.3,'MarkerEdgeColor',[.3 0 0],'MarkerFaceColor','r'); hold on
scatter(-2.5887,0.3307,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
scatter(-2.0153,0.3394,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.9290 0.6940 0.1250]); hold on
%%%%
grid on;lgd=legend('Single SWRO unit','SWRO with ERD', 'Hybrid system', 'Location', 'Northeast');ylabel('FW [m^3/h]','FontSize',16);xlabel('SEC_{net} [kWh/m^3]','FontSize',16);
xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
set(gcf,'color','w'); f = gcf;exportgraphics(f,'Figure_rev.png');

%% max revenue mit pw=2.4 pe=0.14
load('DATA_Pareto.mat'); load('DATA_paretosearch');
pe=0.14;pw=2.4;
%clc,
F=FW_Pareto_2;
S=SEC_Pareto_2;
F=-Y_paretosearch_2(:,2);
S=-Y_paretosearch_2(:,1);
R1=pw.*F+pe.*S.*F;
[a,b]=max(R1);
display(['Revenue = ', num2str(a)]);
d=F(b);
c=S(b);
display(['    SEC = ', num2str(-c)]);
display(['     FW = ', num2str(d)]);

%%
load('DATA_Pareto.mat'); load('DATA_paretosearch');
pe=0.14;pw=2.4;
%clc,
F=FW_Pareto_1;
S=SEC_Pareto_1;
F=-Y_paretosearch_1(:,2);
S=-Y_paretosearch_1(:,1);
R1=pw.*F+pe.*S.*F;
[a,b]=max(R1);
display(['Revenue = ', num2str(a)]);
d=F(b);
c=S(b);
display(['    SEC = ', num2str(-c)]);
display(['     FW = ', num2str(d)]);

%%
load('DATA_Pareto.mat'); load('DATA_paretosearch');
pe=0.14;pw=2.4;
%clc,
F=FW_Pareto_1;
S=SEC_Pareto_1;
F=-Y_paretosearch_3(:,2);
S=-Y_paretosearch_3(:,1);
R1=pw.*F+pe.*S.*F;
[a,b]=max(R1);
display(['Revenue = ', num2str(a)]);
d=F(b);
c=S(b);
display(['    SEC = ', num2str(-c)]);
display(['     FW = ', num2str(d)]);


%% Vortrag scatter plot Pareto + Parameter space
close all;clear all;  
f=figure(1);f.Position = [1200 500 800 500];
load('DATA_Pareto.mat'); load('DATA_paretosearch');
scatter(SEC_Pareto_2,FW_Pareto_2,'MarkerEdgeColor',[.3 0 0],'MarkerFaceColor','r'); hold on
scatter(SEC_Pareto_1,FW_Pareto_1,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
grid on;lgd=legend('Single SWRO unit','SWRO with ERD');ylabel('FW [m^3/h]','FontSize',16);xlabel('SEC_{net} [kWh/m^3]','FontSize',16);%title('Pareto Front','FontSize',16);
xlim([-5.501 -1]); ylim([0 0.48]);
figure(1); set(gcf,'color','w'); f = gcf; exportgraphics(f,'Figure_25.png');

f=figure(2);f.Position=[1200 500 800 500];
scatter(X_Pareto_2(1,:),X_Pareto_2(2,:),'MarkerEdgeColor',[.3 0 0],'MarkerFaceColor','r'); hold on
scatter(X_Pareto_1(1,:),X_Pareto_1(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
grid on;lgd=legend('Single SWRO unit','SWRO with ERD', 'Location','Northwest');xlabel('Inlet pressure [bar]','FontSize',16);ylabel('Outlet pressure [bar]','FontSize',16);%title('Parameter space','FontSize',16);
xlim([30 70]); ylim([30 70]);
figure(2); set(gcf,'color','w'); f = gcf; exportgraphics(f,'Figure_25b.png');

%% SWRO recovery 
load('DATA_Pareto.mat'); load('DATA_paretosearch');
load DATA_REC
for i=201:200
    i,
    x0=X_paretosearch_3(i,:);
    [output1, output2]=fun(x0,3,'sol',1e4,1e-4);
    REC(i)=output1(4);
end
save DATA_REC REC
close all;
scatter3(-Y_paretosearch_3(:,1),-Y_paretosearch_3(:,2),ones(1,200),50*ones(1,200),REC); hold on
cb=colorbar()
view(2)
REC,

%% scatter plot Pareto Front: single SWRO and SWRO+ERD
close all;clear all;  
f=figure(1);f.Position =[1.2737e+03 444.3333 800 500.0000];
load('DATA_Pareto.mat'); load('DATA_paretosearch');
scatter(SEC_Pareto_2,FW_Pareto_2,'MarkerEdgeColor',[.3 0 0],'MarkerFaceColor','r'); hold on
scatter(SEC_Pareto_1,FW_Pareto_1,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
grid on;lgd=legend('Single SWRO unit','SWRO with ERD');ylabel('FW [m^3/h]');xlabel('SEC_{net} [kWh/m^3]');title('Pareto Front');
xlim([-5.501 -1]); ylim([0 0.48]);
%figure(1); set(gcf,'color','w'); f = gcf; exportgraphics(f,'Figure_25.png');

f=figure(2);f.Position=[1.2723e+03 445.6667 1.2347e+03 797.3333];
%left original
subplot(2,2,1)
scatter(SEC_Pareto_2,FW_Pareto_2,'MarkerEdgeColor',[.3 0 0],'MarkerFaceColor','r'); hold on
scatter(SEC_Pareto_1,FW_Pareto_1,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
grid on;lgd=legend('Single SWRO unit','SWRO with ERD');ylabel('FW [m^3/h]');xlabel('SEC_{net} [kWh/m^3]');title('Pareto Front');
xlim([-5.501 -1]); ylim([0 0.48]);
subplot(2,2,3)
scatter(X_Pareto_2(1,:),X_Pareto_2(2,:),'MarkerEdgeColor',[.3 0 0],'MarkerFaceColor','r'); hold on
scatter(X_Pareto_1(1,:),X_Pareto_1(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
grid on;lgd=legend('Single SWRO unit','SWRO with ERD', 'Location','Northwest');xlabel('Inlet pressure');ylabel('Outlet pressure');title('Parameter space');
xlim([30 70]); ylim([30 70]);
%right paretosearch
subplot(2,2,2)
scatter(-Y_paretosearch_2(:,1),-Y_paretosearch_2(:,2),'MarkerEdgeColor',[.3 0 0],'MarkerFaceColor','m'); hold on
scatter(-Y_paretosearch_1(:,1),-Y_paretosearch_1(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','c'); hold on
%%%%%
%scatter(SEC_Pareto_2,FW_Pareto_2,'MarkerEdgeColor',[.3 0 0],'MarkerFaceColor','r'); hold on
%scatter(SEC_Pareto_1,FW_Pareto_1,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
%%%%%
grid on;lgd=legend('Single SWRO unit','SWRO with ERD');ylabel('FW [m^3/h]');xlabel('SEC_{net} [kWh/m^3]');title('Pareto Front');
xlim([-5.501 -1]); ylim([0 0.48]);
subplot(2,2,4)
scatter(X_paretosearch_2(:,1),X_paretosearch_2(:,2),'MarkerEdgeColor',[.3 0 0],'MarkerFaceColor','m'); hold on
scatter(X_paretosearch_1(:,1),X_paretosearch_1(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','c'); hold on
grid on;lgd=legend('Single SWRO unit','SWRO with ERD', 'Location','Northwest');xlabel('Inlet pressure');ylabel('Outlet pressure');title('Parameter space');
xlim([30 70]); ylim([30 70]);


%% Paretosearch algorithim
%
close all; clear all; clc
rng default % For reproducibility
A= []; b= []; Aeq=[]; beq=[]; lb = [30;10;10;1.05;1.5]; ub = [70;20;20;2;1.5];
% initial points:
load('DATA_Approx.mat');
X0 = X_approx';
%remove infeasible poitns of X0:
z=length(X_approx(1,:));
for i=1:z
    k=0;
    for j=1:4
    if X_approx(j,z+1-i) > ub(j); k=1;end
    if X_approx(j,z+1-i) < lb(j); k=1;end
    end
    if k==1;X0(z+1-i,:)=[]; end
end
%
load('DATA_paretosearch');
% options
% MaxFunctionEvaluations = {'3000*(numberOfVariables+numberOfObjectives)'} for paretosearch
% suggested: 3000*6=18000, i.e around 1 week computation
%
% first i test with: 5000
%
options = optimoptions('paretosearch','ParetoSetSize',200, 'InitialPoints',X0,'Display','iter', 'MaxFunctionEvaluations',5000);
% x = paretosearch(fun,nvars,A,b,Aeq,beq,lb,ub,nonlcon,options)
[X_paretosearch_3, Y_paretosearch_3] = paretosearch(@(x)fun(x,3,'Pareto',1e4,1e-4),5,A,b,Aeq,beq,lb,ub,@(x)nonlcon(x,'default'),options); 
save DATA_paretosearch X_paretosearch_1 Y_paretosearch_1 X_paretosearch_2 Y_paretosearch_2 X_paretosearch_3 Y_paretosearch_3

% how to make more accurate:
%
% MeshTolerance             Nonnegative scalar | {1e-6}
% ParetoSetChangeTolerance  Nonnegative scalar | {1e-4}
%
%

%Iter   F-count   NumSolutions  Spread       Volume 
%   0        71        55          -         3.1712e+00
%   1       164        66       8.8451e-01   3.2579e+00
%
%   .
%   .
%   .
%  13      1045       180       6.0705e-01   3.4878e+00
%  14      1275       200       5.7480e-01   3.4882e+00

%Pareto set found that satisfies the constraints. 

%Optimization completed because the relative change in the volume of the Pareto set 
%is less than 'options.ParetoSetChangeTolerance' and constraints are satisfied to within 
%'options.ConstraintTolerance'.


%% plot
close all;clear all;clc
f=figure(1);f.Position=[515 434.3333 1992 808.6667];
%left original
load('DATA_Pareto.mat');load('DATA_paretosearch.mat');load('DATA_Approx.mat');
subplot(2,3,1)
scatter(SEC_Pareto_2,FW_Pareto_2,'MarkerEdgeColor',[.3 0 0],'MarkerFaceColor','r'); hold on
scatter(SEC_Pareto_1,FW_Pareto_1,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
load('DATA_Approx.mat');
scatter(SEC_approx, FW_approx, 'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','y'); hold on
grid on;lgd=legend('Single SWRO unit','SWRO with ERD');ylabel('FW [m^3/h]');xlabel('SEC_{net} [kWh/m^3]');title('Preto front');
xlim([-5.501 -1]); ylim([0 0.48]);view(2);

subplot(2,3,4)
scatter(X_Pareto_2(1,:),X_Pareto_2(2,:),'MarkerEdgeColor',[.3 0 0],'MarkerFaceColor','r'); hold on
scatter(X_Pareto_1(1,:),X_Pareto_1(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
grid on;lgd=legend('Pareto front:  Single SWRO unit','Pareto front:  SWRO with ERD', 'Location','Northwest');xlabel('Inlet pressure');ylabel('Outlet pressure');title('Parameter space');

%right paretosearch
subplot(2,3,2)
scatter(-Y_paretosearch_2(:,1),-Y_paretosearch_2(:,2),'MarkerEdgeColor',[.3 0 0],'MarkerFaceColor','m'); hold on
scatter(-Y_paretosearch_1(:,1),-Y_paretosearch_1(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','c'); hold on
scatter(-Y_paretosearch_3(:,1),-Y_paretosearch_3(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.9290 0.6940 0.1250]); hold on
%%%%%
scatter(SEC_Pareto_2,FW_Pareto_2,'MarkerEdgeColor',[.3 0 0],'MarkerFaceColor','r'); hold on
scatter(SEC_Pareto_1,FW_Pareto_1,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
scatter(SEC_approx, FW_approx, 'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','y'); hold on
%%%%%
grid on;lgd=legend('Single SWRO unit','SWRO with ERD');ylabel('FW [m^3/h]');xlabel('SEC_{net} [kWh/m^3]');title('Pareto front');
xlim([-5.501 -1]); ylim([0 0.48]);view(2);
%right paretosearch
subplot(2,3,3)
scatter(-Y_paretosearch_2(:,1),-Y_paretosearch_2(:,2),'MarkerEdgeColor',[.3 0 0],'MarkerFaceColor','m'); hold on
scatter(-Y_paretosearch_1(:,1),-Y_paretosearch_1(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','c'); hold on
scatter(-Y_paretosearch_3(:,1),-Y_paretosearch_3(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.9290 0.6940 0.1250]); hold on
%%%%%
grid on;lgd=legend('Single SWRO unit','SWRO with ERD');ylabel('FW [m^3/h]');xlabel('SEC_{net} [kWh/m^3]');title('Pareto front');
xlim([-5.501 -1]); ylim([0 0.48]);view(2);
subplot(2,3,5)
scatter(X_paretosearch_2(:,1),X_paretosearch_2(:,2),'MarkerEdgeColor',[.3 0 0],'MarkerFaceColor','m'); hold on
scatter(X_paretosearch_1(:,1),X_paretosearch_1(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','c'); hold on
grid on;lgd=legend('Single SWRO unit','SWRO with ERD', 'Location','Northwest');xlabel('Inlet pressure');ylabel('Outlet pressure');title('Parameter space');
xlim([30 70]); ylim([30 70]);
% hybrid initial conditions
z=subplot(2,3,6);
%scatter3(X_paretosearch_3(:,1),X_paretosearch_3(:,2),X_paretosearch_3(:,3),[],X_paretosearch_3(:,4));hold on;xlabel('SWRO inlet');ylabel('PRO outlet');zlabel('PRO inlet');a=colorbar;a.Label.String = 'PRO fresh inlet';title('Parameter space');
scatter3(X_paretosearch_3(:,1),X_paretosearch_3(:,2),X_paretosearch_3(:,3),[],-Y_paretosearch_3(:,2));hold on;xlabel('SWRO inlet');ylabel('PRO outlet');zlabel('PRO inlet');a=colorbar;a.Label.String = 'FW';colormap(z,cool);title('Parameter space');



%% Paretosearch algorithim
%
close all; clear all; clc
rng default % For reproducibility
A= []; b= []; Aeq=[]; beq=[]; lb = [30;10;10;1.05;1.5]; ub = [70;20;20;2;1.5];
% initial points:
load('DATA_Approx.mat');
X0 = X_approx';
%remove infeasible poitns of X0:
z=length(X_approx(1,:));
for i=1:z
    k=0;
    for j=1:4
    if X_approx(j,z+1-i) > ub(j); k=1;end
    if X_approx(j,z+1-i) < lb(j); k=1;end
    end
    if k==1;X0(z+1-i,:)=[]; end
end
%
load('DATA_paretosearch');
load('Data_neu.mat')
% options
% MaxFunctionEvaluations = {'3000*(numberOfVariables+numberOfObjectives)'} for paretosearch
% suggested: 3000*6=18000, i.e around 1 week computation
%
% first i test with: 5000
% ParetoSetChangeTolerance  Nonnegative scalar | {1e-8}
%
options = optimoptions('paretosearch','ParetoSetSize',200, 'InitialPoints',X0,'Display','iter', 'MaxFunctionEvaluations',5000,'ParetoSetChangeTolerance',1e-8);
% x = paretosearch(fun,nvars,A,b,Aeq,beq,lb,ub,nonlcon,options)
[X_neu, Y_neu] = paretosearch(@(x)fun(x,3,'Pareto',1e4,1e-4),5,A,b,Aeq,beq,lb,ub,@(x)nonlcon(x,'default'),options); 
save Data_neu.mat X_neu Y_neu

%Iter   F-count   NumSolutions  Spread       Volume 
%   
%  31      4213       200       4.6796e-01   3.6505e+00
%  32      4393       200       4.9537e-01   3.6510e+00
%
%Pareto set found that satisfies the constraints. 
%
%Optimization completed because the relative change in the volume of the Pareto set 
%is less than 'options.ParetoSetChangeTolerance' and constraints are satisfied to within 
%'options.ConstraintTolerance'.
%

%% new plot
close all;clear all;clc
f=figure(1);f.Position=[515 434.3333 1992 808.6667];
load('Data_neu.mat')
%left original
load('DATA_Pareto.mat');load('DATA_paretosearch.mat');load('DATA_Approx.mat');
subplot(2,3,1)
scatter(SEC_Pareto_2,FW_Pareto_2,'MarkerEdgeColor',[.3 0 0],'MarkerFaceColor','r'); hold on
scatter(SEC_Pareto_1,FW_Pareto_1,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
load('DATA_Approx.mat');
scatter(SEC_approx, FW_approx, 'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','y'); hold on
grid on;lgd=legend('Single SWRO unit','SWRO with ERD', 'Hybrid','Location','SouthWest');ylabel('FW [m^3/h]');xlabel('SEC_{net} [kWh/m^3]');title('Pareto front (solving optimisation/ simulations)');
xlim([-5.501 -.74]); ylim([0 0.48]);view(2);

subplot(2,3,4)
scatter(X_Pareto_2(1,:),X_Pareto_2(2,:),'MarkerEdgeColor',[.3 0 0],'MarkerFaceColor','r'); hold on
scatter(X_Pareto_1(1,:),X_Pareto_1(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
grid on;lgd=legend('Single SWRO unit','SWRO with ERD', 'Location','Northwest');xlabel('Inlet pressure');ylabel('Outlet pressure');title('Parameter space');

%right paretosearch
subplot(2,3,2)
scatter(-Y_paretosearch_2(:,1),-Y_paretosearch_2(:,2),'MarkerEdgeColor',[.3 0 0],'MarkerFaceColor','m'); hold on
scatter(-Y_paretosearch_1(:,1),-Y_paretosearch_1(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','c'); hold on
%scatter(-Y_paretosearch_3(:,1),-Y_paretosearch_3(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.9290 0.6940 0.1250]); hold on
scatter(-Y_neu(:,1),-Y_neu(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.9290 0.6940 0.1250]); hold on
%%%%%
scatter(SEC_Pareto_2,FW_Pareto_2,'MarkerEdgeColor',[.3 0 0],'MarkerFaceColor','r'); hold on
scatter(SEC_Pareto_1,FW_Pareto_1,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
scatter(SEC_approx, FW_approx, 'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','y'); hold on
%%%%%
grid on;ylabel('FW [m^3/h]');xlabel('SEC_{net} [kWh/m^3]');title('Pareto front');%lgd=legend( 'Paretosearch: Single SWRO', 'Paretosearch: SWRO+ERD', 'Paretosearch: Hybrid','Single SWRO unit','SWRO with ERD', 'Hybrid', 'Location','SouthWest');
xlim([-5.501 -.74]); ylim([0 0.48]);view(2);
%right paretosearch
subplot(2,3,3)
scatter(-Y_paretosearch_2(:,1),-Y_paretosearch_2(:,2),'MarkerEdgeColor',[.3 0 0],'MarkerFaceColor','m'); hold on
scatter(-Y_paretosearch_1(:,1),-Y_paretosearch_1(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','c'); hold on
%scatter(-Y_paretosearch_3(:,1),-Y_paretosearch_3(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.9290 0.6940 0.1250]); hold on
scatter(-Y_neu(:,1),-Y_neu(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.9290 0.6940 0.1250]); hold on
%%%%%
grid on;lgd=legend('Single SWRO unit','SWRO with ERD', 'Hybrid','Location','SouthWest');ylabel('FW [m^3/h]');xlabel('SEC_{net} [kWh/m^3]');title('Pareto front (using paretosearch)');
xlim([-5.501 -.74]); ylim([0 0.48]);view(2);
subplot(2,3,5)
scatter(X_paretosearch_2(:,1),X_paretosearch_2(:,2),'MarkerEdgeColor',[.3 0 0],'MarkerFaceColor','m'); hold on
scatter(X_paretosearch_1(:,1),X_paretosearch_1(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','c'); hold on
grid on;lgd=legend('Single SWRO unit','SWRO with ERD', 'Location','Northwest');xlabel('Inlet pressure');ylabel('Outlet pressure');title('Parameter space');
xlim([30 70]); ylim([30 70]);
% hybrid initial conditions
z=subplot(2,3,6);
%scatter3(X_paretosearch_3(:,1),X_paretosearch_3(:,2),X_paretosearch_3(:,3),[],-Y_paretosearch_3(:,2));hold on;xlabel('SWRO inlet');ylabel('PRO outlet');zlabel('PRO inlet');a=colorbar;a.Label.String = 'FW';colormap(z,cool);title('Parameter space');
scatter3(X_neu(:,1),X_neu(:,2),X_neu(:,3),[],-Y_neu(:,2));hold on;xlabel('SWRO inlet');ylabel('PRO outlet');zlabel('PRO inlet');a=colorbar;a.Label.String = 'FW';colormap(z,cool);title('Parameter space');

%% plots für Vortrag
close all;clear all;clc; f=figure(1); f.Position = [1200 500 800 500];
% Hybrid (Simulation)
k=1;
load neu_Data Xneu Yneu;
if k==1
    % remove all points where pressures constraints aren't satisfied
    for i=1:length(Xneu)
        k2=0;
        if Xneu(1,i)>70;k2=1; end
        if Xneu(2,i)>20;k2=1; end
        if Xneu(3,i)>20;k2=1; end
        if Xneu(4,i)<1;k2=1; end
        if k2==1; Yneu(:,i)=[NaN; NaN]; end
        ZZ=[6887,17180,14668,17305,16716,7313,7625,9599];for z=1:length(ZZ);Yneu(:,ZZ(z))=[NaN; NaN];end
        if Yneu(1,i)>-2.138; if Yneu(2,i)> 0.353; Yneu(:,i)=[NaN; NaN];end; end
    end
    scatter(Yneu(1,:),Yneu(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','m'); hold on; 
end
if k==2
    % remove all points where pressures constraints aren't satisfied
    for i=1:length(Xneu)
        k2=0;
        if Xneu(1,i)>70;k2=1; end
        if Xneu(2,i)>20;k2=1; end
        if Xneu(3,i)>20;k2=1; end
        if Xneu(4,i)<1;k2=1; end
        if k2==1; Yneu(:,i)=[NaN; NaN]; end
        ZZ=[6887,17180,14668,17305,16716,7313,7625,9599];for z=1:length(ZZ);Yneu(:,ZZ(z))=[NaN; NaN];end
        if Yneu(1,i)>-2.138; if Yneu(2,i)> 0.353; Yneu(:,i)=[NaN; NaN];end; end
    end
    scatter(Yneu(1,:),Yneu(2,:), 'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','m'); hold on;
    %scatter3(Yneu(1,:),Yneu(2,:), 1:1:19200,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','m'); hold on;
end
load('DATA_Pareto.mat'); 
scatter(SEC_Pareto_2,FW_Pareto_2,'MarkerEdgeColor',[.3 0 0],'MarkerFaceColor','r'); hold on
scatter(SEC_Pareto_1,FW_Pareto_1,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
if k<2
load('DATA_Approx.mat');
scatter(SEC_approx, FW_approx, 'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','y'); hold on
end
xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
if k==2; grid on;lgd=legend('Hybrid simulations','Single SWRO unit','SWRO with ERD', 'Location', 'Northeast');ylabel('FW [m^3/h]','FontSize',16);xlabel('SEC_{net} [kWh/m^3]','FontSize',16);end
if k==1; grid on;lgd=legend('Hybrid simulations','Single SWRO unit','SWRO with ERD', 'Nondominated points', 'Location', 'Northeast');ylabel('FW [m^3/h]','FontSize',16);xlabel('SEC_{net} [kWh/m^3]','FontSize',16);end
if k==0; grid on;lgd=legend('Single SWRO unit','SWRO with ERD', 'Nondominated points', 'Location', 'Northeast');ylabel('FW [m^3/h]','FontSize',16);xlabel('SEC_{net} [kWh/m^3]','FontSize',16);end
%
if k==2;figure(1); set(gcf,'color','w'); f = gcf;exportgraphics(f,'Figure_21a.png'); end
if k==1;figure(1); set(gcf,'color','w'); f = gcf;exportgraphics(f,'Figure_21b.png'); end
if k==0;figure(1); set(gcf,'color','w'); f = gcf;exportgraphics(f,'Figure_21c.png'); end
%
load('Data_neu.mat');load('DATA_Pareto.mat');load('DATA_paretosearch.mat');load('DATA_Approx.mat');
f=figure(2); f.Position = [1200 500 800 500];
scatter(-Y_paretosearch_2(:,1),-Y_paretosearch_2(:,2),'MarkerEdgeColor',[.3 0 0],'MarkerFaceColor','m'); hold on
scatter(-Y_paretosearch_1(:,1),-Y_paretosearch_1(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','c'); hold on
%scatter(-Y_paretosearch_3(:,1),-Y_paretosearch_3(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.9290 0.6940 0.1250]); hold on
scatter(-Y_neu(:,1),-Y_neu(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.9290 0.6940 0.1250]); hold on
%%%%
grid on;lgd=legend('Paretosearch: Single SWRO unit','Paretosearch: SWRO with ERD', 'Paretosearch: Hybrid system', 'Location', 'Northeast');ylabel('FW [m^3/h]','FontSize',16);xlabel('SEC_{net} [kWh/m^3]','FontSize',16);
xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
figure(2); set(gcf,'color','w'); f = gcf;exportgraphics(f,'Figure_21d.png');
f=figure(3); f.Position = [1200 500 800 500];
scatter(-Y_paretosearch_2(:,1),-Y_paretosearch_2(:,2),'MarkerEdgeColor',[.3 0 0],'MarkerFaceColor','m'); hold on
scatter(-Y_paretosearch_1(:,1),-Y_paretosearch_1(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','c'); hold on
%scatter(-Y_paretosearch_3(:,1),-Y_paretosearch_3(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.9290 0.6940 0.1250]); hold on
scatter(-Y_neu(:,1),-Y_neu(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.9290 0.6940 0.1250]); hold on
%%%%%
scatter(SEC_Pareto_2,FW_Pareto_2,'MarkerEdgeColor',[.3 0 0],'MarkerFaceColor','r'); hold on
scatter(SEC_Pareto_1,FW_Pareto_1,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
%scatter(SEC_approx, FW_approx, 'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','y'); hold on
%%%%%
grid on;lgd=legend('Single SWRO unit','SWRO with ERD', 'Hybrid system', 'Location', 'Northeast');ylabel('FW [m^3/h]','FontSize',16);xlabel('SEC_{net} [kWh/m^3]','FontSize',16);
xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
figure(3); set(gcf,'color','w'); f = gcf;exportgraphics(f,'Figure_21e.png');
f=figure(4); f.Position = [1200 500 800 500];
scatter(-Y_paretosearch_2(:,1),-Y_paretosearch_2(:,2),'MarkerEdgeColor',[.3 0 0],'MarkerFaceColor','r'); hold on
scatter(-Y_paretosearch_1(:,1),-Y_paretosearch_1(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
%scatter(-Y_paretosearch_3(:,1),-Y_paretosearch_3(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.9290 0.6940 0.1250]); hold on
scatter(-Y_neu(:,1),-Y_neu(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.9290 0.6940 0.1250]); hold on
%%%%%
grid on;lgd=legend('Single SWRO unit','SWRO with ERD', 'Hybrid system', 'Location', 'Northeast');ylabel('FW [m^3/h]','FontSize',16);xlabel('SEC_{net} [kWh/m^3]','FontSize',16);
xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
figure(4); set(gcf,'color','w'); f = gcf;exportgraphics(f,'Figure_21f.png');
close all;clear all;clc; 
f=figure(5); f.Position = [1200 500 800 500];
load neu_Data Xneu Yneu;
load('DATA_Pareto.mat'); 
scatter(SEC_Pareto_2,FW_Pareto_2,'MarkerEdgeColor',[.3 0 0],'MarkerFaceColor','r'); hold on
scatter(SEC_Pareto_1,FW_Pareto_1,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
grid on;lgd=legend('Single SWRO unit','SWRO with ERD', 'Location', 'Northeast');ylabel('FW [m^3/h]','FontSize',16);xlabel('SEC_{net} [kWh/m^3]','FontSize',16);
%
figure(5); set(gcf,'color','w'); f = gcf;exportgraphics(f,'Figure_21.png');



%%
x0=[67.4641,16.3597,10.8084,1.0806,1.5000];
[a,b]=fun(x0,3,'sol',1e4,1e-4),


%% Test Paretosearch algorithim with bad initial data
%
% take small set of random initial data
%
% how good will the paretosearch alg perform
%
close all; clear all; clc
rng default % For reproducibility
A= []; b= []; Aeq=[]; beq=[]; lb = [30;10;10;1.05;1.5]; ub = [70;20;20;2;1.5];
% initial points:
load('DATA_Approx.mat');
X0 = [65.9160,11.4496,16.0490,1.2047,1.5000;
      61.9221,14.9150,10.3609,1.2338,1.5000;
      48.6616,19.6141,17.4707,1.2823,1.5000; 
      65.1053,15.7382,12.7348,1.2955,1.5000;
      67.4641,16.3597,10.8084,1.0806,1.5000];
% options
% MaxFunctionEvaluations = {'3000*(numberOfVariables+numberOfObjectives)'} for paretosearch
% suggested: 3000*6=18000, i.e around 1 week computation
%
% first i test with: 5000
% ParetoSetChangeTolerance  Nonnegative scalar | {1e-8}
%
load Data_neu_test 
options = optimoptions('paretosearch','ParetoSetSize',200, 'InitialPoints',X0,'Display','iter', 'MaxFunctionEvaluations',5000,'ParetoSetChangeTolerance',1e-8);
% x = paretosearch(fun,nvars,A,b,Aeq,beq,lb,ub,nonlcon,options)
[XXX, YYY] = paretosearch(@(x)fun(x,3,'Pareto',1e4,1e-4),5,A,b,Aeq,beq,lb,ub,@(x)nonlcon(x,'default'),options); 
save Data_neu_test XXX YYY

%Iter   F-count   NumSolutions  Spread       Volume 
%
%30      3835       200       5.2373e-01   3.1567e+00
%31      4005       200       9.3853e-01   3.2380e+00

%Pareto set found that satisfies the constraints. 
%
%Optimization completed because the relative change in the volume of the Pareto set 
%is less than 'options.ParetoSetChangeTolerance' and constraints are satisfied to within 
%'options.ConstraintTolerance'.

%% new plot
close all;clear all;clc
f=figure(1);f.Position=[37 631 2.4993e+03 506.6667];
load('Data_neu.mat');load('Data_neu_test.mat') 
%left original
load('DATA_Pareto.mat');load('DATA_paretosearch.mat');load('DATA_Approx.mat');
subplot(1,3,1)
scatter(SEC_Pareto_2,FW_Pareto_2,'MarkerEdgeColor',[.3 0 0],'MarkerFaceColor','r'); hold on
scatter(SEC_Pareto_1,FW_Pareto_1,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
XX=[-4.4450, -4.9748, -2.9145, -2.8798, -5.6522];YY=[0.2680,0.2151,0.1204,0.3051,0.2504];
scatter(XX,YY,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','g'); hold on
load('DATA_Approx.mat');
scatter(SEC_approx, FW_approx, 'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','y'); hold on
grid on;lgd=legend('Single SWRO unit','SWRO with ERD', 'Hybrid','Location','SouthWest');ylabel('FW [m^3/h]');xlabel('SEC_{net} [kWh/m^3]');title('Pareto front (solving optimisation/ simulations)');
xlim([-6 -.74]); ylim([0 0.48]);view(2);
%right paretosearch
subplot(1,3,2)
scatter(-Y_paretosearch_2(:,1),-Y_paretosearch_2(:,2),'MarkerEdgeColor',[.3 0 0],'MarkerFaceColor','m'); hold on
scatter(-Y_paretosearch_1(:,1),-Y_paretosearch_1(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','c'); hold on
%scatter(-Y_paretosearch_3(:,1),-Y_paretosearch_3(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.9290 0.6940 0.1250]); hold on
scatter(-Y_neu(:,1),-Y_neu(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.9290 0.6940 0.1250]); hold on
scatter(-YYY(:,1),-YYY(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','g'); hold on
%%%%%
scatter(SEC_Pareto_2,FW_Pareto_2,'MarkerEdgeColor',[.3 0 0],'MarkerFaceColor','r'); hold on
scatter(SEC_Pareto_1,FW_Pareto_1,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
%%%%%
grid on;ylabel('FW [m^3/h]');xlabel('SEC_{net} [kWh/m^3]');title('Pareto front');%lgd=legend( 'Paretosearch: Single SWRO', 'Paretosearch: SWRO+ERD', 'Paretosearch: Hybrid','Single SWRO unit','SWRO with ERD', 'Hybrid', 'Location','SouthWest');
xlim([-5.501 -.74]); ylim([0 0.48]);view(2);
%right paretosearch
subplot(1,3,3)
scatter(-Y_paretosearch_2(:,1),-Y_paretosearch_2(:,2),'MarkerEdgeColor',[.3 0 0],'MarkerFaceColor','m'); hold on
scatter(-Y_paretosearch_1(:,1),-Y_paretosearch_1(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','c'); hold on
scatter(-Y_neu(:,1),-Y_neu(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.9290 0.6940 0.1250]); hold on
scatter(-YYY(:,1),-YYY(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','g'); hold on
%%%%%
grid on;lgd=legend('Single SWRO unit','SWRO with ERD', 'Hybrid','Location','SouthWest');ylabel('FW [m^3/h]');xlabel('SEC_{net} [kWh/m^3]');title('Pareto front (using paretosearch)');
xlim([-5.501 -.74]); ylim([0 0.48]);view(2);
