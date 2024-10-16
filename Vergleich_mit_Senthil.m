%% Basic simulations for all cases (ideal, ICP, ICP+ECP)
% press [ctrl]+[enter] to run code sections
addpath('Input_DATA','Scaled_model','Unscaled_model','Output_DATA')

%% Vergleich mit Senthil Table 2
clc
output=fun_unscaled([0; 50.47; 48.23],0,'fig',1e4,1e-3);
disp('Operating cond:  P_d(0) = 50.47 bar      P_d(L) = 48.23 bar')
disp(' ')
disp('Experiment:      Recovery = 33.1 %       C_Permeate = 44 ppm')
disp(['Our model:       Recovery = ',num2str(output(4)),' %    C_Permeate = ',num2str(output(9)),' ppm   ==>   J_d(0) = ',num2str(131.7146216211166*output(6)),' kg/sm'])
str1=100*(output(4)-33.1)/33.1;str2=100*(output(9)-44)/44;
if str1 > 0; str1=['+',num2str(str1)]; else str1=num2str(str1); end
if str2 > 0; str2=['+',num2str(str2)]; else str2=num2str(str2); end
warning(['                   ',str1,' %                ',str2,' %'])
disp('_______________')
J_d_in=131.7146216211166*output(6);
output=fun_unscaled([J_d_in; 55.81; 0],0,'sol',1e4,1e-3);
disp(['Operating cond:  P_d(0) = 55.81 bar      J_d(0) = ',num2str(131.7146216211166*output(6)),' kg/sm'] )
disp(' ')
disp('Experiment:      Recovery = 39.6 %       C_Permeate = 49 ppm         P_d(L) = 54.72')
disp(['Our model:       Recovery = ',num2str(output(4)),' %     C_Permeate = ',num2str(output(9)),' ppm    P_d(L) = ',num2str(output(11)),' bar'])
str1=100*(output(4)-39.6)/39.6;str2=100*(output(9)-49)/49;str3=100*(output(11)-54.72)/54.72;
if str1 > 0; str1=['+',num2str(str1)]; else str1=num2str(str1); end
if str2 > 0; str2=['+',num2str(str2)]; else str2=num2str(str2); end
if str3 > 0; str3=['+',num2str(str3)]; else str3=num2str(str3); end
warning(['                   ',str1,' %                 ',str2,' %               ',str3,' %'])
disp('_______________')
output=fun_unscaled([J_d_in; 60.28; 0],0,'sol',1e4,1e-3);
disp(['Operating cond:  P_d(0) = 60.28 bar      J_d(0) = ',num2str(131.7146216211166*output(6)),' kg/sm'] )
disp(' ')
disp('Experiment:      Recovery = 44.5 %       C_Permeate = 52 ppm         P_d(L) = 59.22')
disp(['Our model:       Recovery = ',num2str(output(4)),' %    C_Permeate = ',num2str(output(9)),' ppm    P_d(L) = ',num2str(output(11)),' bar'])
str1=100*(output(4)-44.5)/44.5;str2=100*(output(9)-52)/52;str3=100*(output(11)-59.22)/59.22;
if str1 > 0; str1=['+',num2str(str1)]; else str1=num2str(str1); end
if str2 > 0; str2=['+',num2str(str2)]; else str2=num2str(str2); end
if str3 > 0; str3=['+',num2str(str3)]; else str3=num2str(str3); end
warning(['                   ',str1,' %                 ',str2,' %               ',str3,' %'])
disp('_______________')

%% plot figure
close all;f=figure(1);f.Position = [1.2977e+03 635.6667 1.0987e+03 639.3333];
load("Output_DATA/Vergleich_mit_Senthil.mat")
Z=Z*0.01;ZZ=ZZ/1e6;YYY=YYY/1e6;
for i=1:length(Z)
    Z(i)=433000*Z(i)/(200*Z(i)+433); % kg/m^3
    Z(i)=1000*Z(i); % g/m^3
    ZZ(i)=433000*ZZ(i)/(200*ZZ(i)+433); % kg/m^3
    ZZ(i)=1000*ZZ(i); % g/m^3
    YYY(i)=433000*YYY(i)/(200*YYY(i)+433); % kg/m^3
    YYY(i)=1000*YYY(i); % g/m^3
end
tiledlayout(2,2); nexttile; 
nexttile % Feed pressure vs Recovery & Brine salinity
lw=1.5;lc='k';rc='#0b0fe3';
plot(X, Y, 'Color', lc, 'LineWidth',lw);xlabel('Feed pressure (bar)','Fontsize',10);  ylabel('Recovery (%)','Fontsize',10);ylim([-6 60]);ay=gca;hold on
yyaxis right
plot(X, Z, 'Color', rc, 'LineWidth',lw); ylabel('Reject Concentration (g/m^3)','Fontsize',10); xlim([28 72]);ay.YAxis(2).Exponent = 0;ylim([30000 80000]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
nexttile % Feed pressure vs SEC & Permeate salinity
%plot(X, YY, 'Color', lc, 'LineWidth',lw);xlabel('Feed pressure (bar)','Fontsize',10);  ylabel('SEC (kWh/m^3)','Fontsize',10);ay=gca;hold on
plot(X, abs(QQ)./7/0.9626/7.7210, 'Color', lc, 'LineWidth',lw);xlabel('Feed pressure (bar)','Fontsize',10);ylim([7 56]);ylabel('SEC (W/m^2)','Fontsize',10);ay=gca;hold on
yyaxis right
plot(X, ZZ, 'Color', rc, 'LineWidth',lw); ylabel('Product Concentration (g/m^3)','Fontsize',10); xlim([28 72]);ylim([40 90]);ay.YAxis(2).Exponent = 0;ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
nexttile % Feed flow rate vs Permeate salinity & SEC
plot(XX, YYY, 'Color', lc, 'LineWidth',lw);xlabel('Feed mass flow rate (kg/sm)','Fontsize',10);  ylabel('Product Concentration (g/m^3)','Fontsize',10);ay.YAxis(2).Exponent = 0;ylim([46 85]);ay=gca;hold on
yyaxis right
plot(XX, ZZZ, 'Color', rc, 'LineWidth',lw); ylabel('SEC (kWh/m^3)','Fontsize',10); xlim([.4*0.080316925377597 1.15*0.080316925377597]);ylim([1.6 1.85]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
%plot(XX, abs(QQQ)./7/0.9626/7.7210, 'Color', rc, 'LineWidth',lw); ylabel('SEC (W/m^2)','Fontsize',10); xlim([.4*0.080316925377597 1.15*0.080316925377597]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;

%% save figure
figure(1); set(gcf,'color','w'); f = gcf; exportgraphics(f,'Vergleich_Senthil.png');








%% Calculation for Fig 2
X=30:1:70;XX=linspace(.44*J_d_in,1.1*J_d_in,length(X));
for x=1:length(X)
    out=fun_unscaled([J_d_in; X(x); 0],0,'sol',1e4,1e-3);
    Y(x)= out(4); Z(x)=out(10);Q(x)=out(13);
    YY(x)= -out(1); ZZ(x)=out(9);QQ(x)=out(13);
    out=fun_unscaled([XX(x); 50.47; 0],0,'sol',1e4,1e-3);
    YYY(x)= out(9); ZZZ(x)=-out(1);QQQ(x)=out(13);
end
save Output_DATA/Vergleich_mit_Senthil.mat X Y Z XX YY ZZ YYY ZZZ Q QQ QQQ



