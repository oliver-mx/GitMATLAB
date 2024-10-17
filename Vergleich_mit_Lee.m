%% Basic simulations for all cases (ideal, ICP, ICP+ECP)
% press [ctrl]+[enter] to run code sections
addpath('Input_DATA','Scaled_model','Unscaled_model','Output_DATA')

%% test
% flow rates bei C=0%
%
% [9 10 12 13 14].*(1000/60/17.7165/3600)
% = 0.0024    0.0026    0.0031    0.0034    0.0037 [kg/ms]
%
clc,
[a,b]=fun_scaled([-0.0024, 5, 1.000003, 0.05],-1,'fig',1e4,1e-4);
ev(a,[5 6 12 17])

%% Plot figure 10
close all;f=figure(1);f.Position = [597 635 1024 604.6667];
load Output_DATA/Vergleich_mit_Lee.mat
lw=1.3;ms=9;
plot([9 10 12 13 14],Y(1,:)./100,'k-o','Linewidth',lw,'MarkerSize',ms,'MarkerFaceColor','black'); hold on % 5 bar
plot([9 10 12 13 14],Y(2,:)./100,'k-o','Linewidth',lw,'MarkerSize',ms,'MarkerFaceColor','white'); hold on % 10 bar
plot([9 10 12 13 14],Y(3,:)./100,'k-v','Linewidth',lw,'MarkerSize',ms,'MarkerFaceColor','black'); hold on % 15 bar
plot([9 10 12 13 14],Y(4,:)./100,'k-^','Linewidth',lw,'MarkerSize',ms,'MarkerFaceColor','white'); hold on % 20 bar
xlim([8 15]);ylim([0 1]);
yticks([0 0.2 0.4 0.6 0.8 1]);yticklabels({'0.0','0.2','0.4','0.6','0.8','1.0'})
xlabel('PRO feed flow rate (m^3/h)','Fontsize',14);ylabel('PRO recovery','Fontsize',14)
lgd = legend('5 bar','10 bar','15 bar','20 bar', 'Location','eastoutside');
set(lgd,'FontSize',12);
%figure(1); set(gcf,'color','w'); f = gcf; exportgraphics(f,'Vergleich_Lee.png');

%% Calculation for Fig 10
P = [5 10 15 20];
total_flow = [-0.0024 -0.0026 -0.0031 -0.0034 -0.0037]; % [9 10 12 13 14] m^3/h --> kg/s/m
C=0.06;
Y=zeros(4,5); YY=zeros(20,21);
X=1.000003.*ones(4,5); % 1.000003 works well
k=0;
for i=1:5
    for j=1:4
        k=k+1;
        out=fun_scaled([total_flow(i), P(j), X(j,i), C],-1,'sol',1e4,1e-4);
        YY(k,:)=out;
        Y(j,i)=out(5);
    end
end
disp(Y)
save Output_DATA/Vergleich_mit_Lee.mat X Y YY

%% PROBLEM with PRO fresh side pressure
P = [5 10 15 20]; C=0.06; total_flow = [-0.0024 -0.0026 -0.0031 -0.0034 -0.0037];
Y=zeros(4,5); YY=zeros(20,21);close all;clc
f=figure(1);f.Position = [692.3333 273 1.6879e+03 1.0574e+03]; tiledlayout(2,2)
load Output_DATA/Vergleich_mit_Lee.mat
nexttile; lw=1.1;ms=7;
plot([9 10 12 13 14],Y(1,:)./100,'k-o','Linewidth',lw,'MarkerSize',ms,'MarkerFaceColor','black'); hold on % 5 bar
plot([9 10 12 13 14],Y(2,:)./100,'k-o','Linewidth',lw,'MarkerSize',ms,'MarkerFaceColor','white'); hold on % 10 bar
plot([9 10 12 13 14],Y(3,:)./100,'k-v','Linewidth',lw,'MarkerSize',ms,'MarkerFaceColor','black'); hold on % 15 bar
plot([9 10 12 13 14],Y(4,:)./100,'k-^','Linewidth',lw,'MarkerSize',ms,'MarkerFaceColor','white'); hold on % 20 bar
xlim([8 15]);ylim([0 1]); yticks([0 0.2 0.4 0.6 0.8 1]);yticklabels({'0.0','0.2','0.4','0.6','0.8','1.0'})
xlabel('PRO feed flow rate (m^3/h)','Fontsize',14);ylabel('PRO recovery','Fontsize',14);lgd = legend('5 bar','10 bar','15 bar','20 bar', 'Location','eastoutside');set(lgd,'FontSize',12);
lgd.Title.String = 'Pd_f(0)=100000.3 Pa'; lgd.Title.FontSize = 12;
% _______________________________________________________________
nexttile
X=1.000001.*ones(4,5);
for i=1:5
    for j=1:4
        out=fun_scaled([total_flow(i), P(j), X(j,i), C],-1,'sol',1e4,1e-4);Y(j,i)=out(5);
    end
end
plot([9 10 12 13 14],Y(1,:)./100,'k-o','Linewidth',lw,'MarkerSize',ms,'MarkerFaceColor','black'); hold on % 5 bar
plot([9 10 12 13 14],Y(2,:)./100,'k-o','Linewidth',lw,'MarkerSize',ms,'MarkerFaceColor','white'); hold on % 10 bar
plot([9 10 12 13 14],Y(3,:)./100,'k-v','Linewidth',lw,'MarkerSize',ms,'MarkerFaceColor','black'); hold on % 15 bar
plot([9 10 12 13 14],Y(4,:)./100,'k-^','Linewidth',lw,'MarkerSize',ms,'MarkerFaceColor','white'); hold on % 20 bar
xlim([8 15]);ylim([0 1]); yticks([0 0.2 0.4 0.6 0.8 1]);yticklabels({'0.0','0.2','0.4','0.6','0.8','1.0'})
xlabel('PRO feed flow rate (m^3/h)','Fontsize',14);ylabel('PRO recovery','Fontsize',14);lgd = legend('5 bar','10 bar','15 bar','20 bar', 'Location','eastoutside');set(lgd,'FontSize',12);
lgd.Title.String = 'Pd_f(0)=100000.1 Pa'; lgd.Title.FontSize = 12;
% _______________________________________________________________
nexttile
X=1.00001.*ones(4,5);
for i=1:5
    for j=1:4
        out=fun_scaled([total_flow(i), P(j), X(j,i), C],-1,'sol',1e4,1e-4);Y(j,i)=out(5);
    end
end
plot([9 10 12 13 14],Y(1,:)./100,'k-o','Linewidth',lw,'MarkerSize',ms,'MarkerFaceColor','black'); hold on % 5 bar
plot([9 10 12 13 14],Y(2,:)./100,'k-o','Linewidth',lw,'MarkerSize',ms,'MarkerFaceColor','white'); hold on % 10 bar
plot([9 10 12 13 14],Y(3,:)./100,'k-v','Linewidth',lw,'MarkerSize',ms,'MarkerFaceColor','black'); hold on % 15 bar
plot([9 10 12 13 14],Y(4,:)./100,'k-^','Linewidth',lw,'MarkerSize',ms,'MarkerFaceColor','white'); hold on % 20 bar
xlim([8 15]);ylim([0 1]); yticks([0 0.2 0.4 0.6 0.8 1]);yticklabels({'0.0','0.2','0.4','0.6','0.8','1.0'})
xlabel('PRO feed flow rate (m^3/h)','Fontsize',14);ylabel('PRO recovery','Fontsize',14);lgd = legend('5 bar','10 bar','15 bar','20 bar', 'Location','eastoutside');set(lgd,'FontSize',12);
lgd.Title.String = 'Pd_f(0)=100001 Pa'; lgd.Title.FontSize = 12;
% _______________________________________________________________
nexttile
X=1.00005.*ones(4,5);
for i=1:5
    for j=1:4
        out=fun_scaled([total_flow(i), P(j), X(j,i), C],-1,'sol',1e4,1e-4);Y(j,i)=out(5);
    end
end
plot([9 10 12 13 14],Y(1,:)./100,'k-o','Linewidth',lw,'MarkerSize',ms,'MarkerFaceColor','black'); hold on % 5 bar
plot([9 10 12 13 14],Y(2,:)./100,'k-o','Linewidth',lw,'MarkerSize',ms,'MarkerFaceColor','white'); hold on % 10 bar
plot([9 10 12 13 14],Y(3,:)./100,'k-v','Linewidth',lw,'MarkerSize',ms,'MarkerFaceColor','black'); hold on % 15 bar
plot([9 10 12 13 14],Y(4,:)./100,'k-^','Linewidth',lw,'MarkerSize',ms,'MarkerFaceColor','white'); hold on % 20 bar
xlim([8 15]);ylim([0 1]); yticks([0 0.2 0.4 0.6 0.8 1]);yticklabels({'0.0','0.2','0.4','0.6','0.8','1.0'})
xlabel('PRO feed flow rate (m^3/h)','Fontsize',14);ylabel('PRO recovery','Fontsize',14);lgd = legend('5 bar','10 bar','15 bar','20 bar', 'Location','eastoutside');set(lgd,'FontSize',12);
lgd.Title.String = 'Pd_f(0)=100005 Pa'; lgd.Title.FontSize = 12;
%
%figure(1); set(gcf,'color','w'); f = gcf; exportgraphics(f,'Vergleich_Lee2.png');















