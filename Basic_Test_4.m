%% Simulations for Hybrid I pareto front with varying L^{PRO}
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
load Test_4_DATA.mat
% define new input data
n2=202400;
X_test=[repmat(X2_pareto,1,n2/200); zeros(3,n2)];Y_test=zeros(28,n2);k=0;
% 4*11*23 = 1012
P=[5 10 15 20]; % 4
p=[1.1 1.05 1.01 1.005 1.001 1.0005 1.0001 1.00005 1.00001 1.000005 1.000001]; % 11
L= 0.9626.*[1 2 3 4 5 6 7 8 9 10 12 14 16 18 20 22 24 26 28 30 40 50 60]; % 23
% 200*1012
for i=1:200
    for j=1:length(P)
        for t=1:length(p)
            for z=1:length(L)
                k=k+1;
                X_test(3,k)=P(j);
                X_test(4,k)=p(t);
                X_test(5,k)=L(z);
            end
        end
    end
end
% simulate the hybrid system (low accuracy!!!!!!)
option_mesh = 1e3; option_BVP = 1e-3; option_data = .6;
parfor i=1:n2
    Y_test(:,i)=fun_scaled(X_test(:,i),option_data,'sol',option_mesh,option_BVP);
end
save Test_4_DATA.mat X_test Y_test P p L n2
system('git add .'); system('git commit -m "PRO length test n \approx 200k"');system('git push https://github.com/oliver-mx/GitMATLAB.git');


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
if sum(R)==152
lgd=legend('L^{PRO} = 0.9626m \times 1','L^{PRO} = 0.9626m \times 5','L^{PRO} = 0.9626m \times 7','L^{PRO} = 0.9626m \times 10','L^{PRO} = 0.9626m \times 20','L^{PRO} = 0.9626m \times 30','L^{PRO} = 0.9626m \times 40','L^{PRO} = 0.9626m \times 50','L^{PRO} = 0.9626m \times 60','Case1: SWRO','Case2: SWRO+ERD','Location', 'bestoutside');
title(lgd,'Legend')
end
if sum(R)==253
lgd=legend('Hybrid: L^{RO} \approx 6.7m,  L^{PRO} \approx 57.7m','SWRO: L^{RO} \approx 6.7m','SWRO+ERD L^{RO} \approx 6.7m','Location', 'bestoutside');
title(lgd,'Legend')
end
figure(1); set(gcf,'color','w'); f = gcf; exportgraphics(f,'PRO_length_test.png');





















