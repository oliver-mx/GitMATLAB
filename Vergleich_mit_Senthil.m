%% Basic simulations for all cases (ideal, ICP, ICP+ECP)
% press [ctrl]+[enter] to run code sections
addpath('Input_DATA','Scaled_model','Unscaled_model','Output_DATA')

%% Vergleich mit Senthil Table 2
% default simulation:
clc
output=fun_unscaled([0; 50.47; 48.23],0,'sol',1e4,1e-3);
%
disp('Operating cond:  P_d(0) = 50.47 bar      P_d(L) = 48.23 bar')
disp(' ')
disp('Experiment:      Recovery = 33.1 %       C_Permeate = 44 ppm')
disp(['Our model:       Recovery = ',num2str(output(4)),' %    C_Permeate = ',num2str(output(9)),' ppm   ==>   J_d(0) = ',num2str(131.7146216211166*output(6)),' kg/sm'])
disp('_______________')
%
J_d_in=131.7146216211166*output(6);
output=fun_unscaled([J_d_in; 55.81; 0],0,'sol',1e4,1e-3);
disp(['Operating cond:  P_d(0) = 55.81 bar      J_d(0) = ',num2str(131.7146216211166*output(6)),' kg/sm'] )
disp(' ')
disp('Experiment:      Recovery = 39.6 %       C_Permeate = 49 ppm         P_d(L) = 54.72')
disp(['Our model:       Recovery = ',num2str(output(4)),' %     C_Permeate = ',num2str(output(9)),' ppm    P_d(L) = ',num2str(output(11)),' bar'])
%
disp('_______________')
output=fun_unscaled([J_d_in; 60.28; 0],0,'sol',1e4,1e-3);
disp(['Operating cond:  P_d(0) = 60.28 bar      J_d(0) = ',num2str(131.7146216211166*output(6)),' kg/sm'] )
disp(' ')
disp('Experiment:      Recovery = 44.5 %       C_Permeate = 52 ppm         P_d(L) = 59.22')
disp(['Our model:       Recovery = ',num2str(output(4)),' %    C_Permeate = ',num2str(output(9)),' ppm    P_d(L) = ',num2str(output(11)),' bar'])
disp('_______________')

%% Vergleich mit Senthil Fig 2
X=30:1:70;XX=linspace(.44*0.080316925377597,1.1*0.080316925377597,length(X));
for x=1:length(X)
    out=fun_unscaled([J_d_in; X(x); 0],0,'sol',1e4,1e-3);
    Y(x)= out(4); Z(x)=out(10);
    YY(x)= -out(1); ZZ(x)=out(9);
    out=fun_unscaled([XX(x); 50.47; 0],0,'sol',1e4,1e-3);
    YYY(x)= out(9); ZZZ(x)=-out(1);
end
%%
close all;f=figure(1);f.Position = [1.2977e+03 635.6667 1.0987e+03 639.3333];
tiledlayout(2,2); nexttile; 
nexttile % Feed pressure vs Recovery & Brine salinity
lw=1.5;lc='k';rc='#0b0fe3';
plot(X, Y, 'Color', lc, 'LineWidth',lw);xlabel('Feed pressure (bar)','Fontsize',10);  ylabel('Recovery (%)','Fontsize',10);ay=gca;hold on
yyaxis right
plot(X, Z, 'Color', rc, 'LineWidth',lw); ylabel('Reject Concentration (%)','Fontsize',10); xlim([28 72]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
nexttile % Feed pressure vs SEC & Permeate salinity
plot(X, YY, 'Color', lc, 'LineWidth',lw);xlabel('Feed pressure (bar)','Fontsize',10);  ylabel('SEC (kWh/m^3)','Fontsize',10);ay=gca;hold on
yyaxis right
plot(X, ZZ, 'Color', rc, 'LineWidth',lw); ylabel('Product Concentration (ppm)','Fontsize',10); xlim([28 72]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
nexttile % Feed flow rate vs Permeate salinity & SEC
plot(XX, YYY, 'Color', lc, 'LineWidth',lw);xlabel('Feed mass flow rate (kg/sm)','Fontsize',10);  ylabel('Product Concentration (ppm)','Fontsize',10);ay=gca;hold on
yyaxis right
plot(XX, ZZZ, 'Color', rc, 'LineWidth',lw); ylabel('SEC (kWh/m^3)','Fontsize',10); xlim([.4*0.080316925377597 1.15*0.080316925377597]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;

%% save figure
figure(1); set(gcf,'color','w'); f = gcf; exportgraphics(f,'Vergleich_Senthil.png');

%% load/save data for plot
%load("Output_DATA/Senthil_verleich.mat")
save Output_DATA/Senthil_verleich.mat X Y Z XX YY ZZ YYY ZZZ


