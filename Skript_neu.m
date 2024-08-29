%% Neues Skript nur für optimale PRO Längenverhältnis
%
% Ansatz:
% für jedes L^{PRO} eine Pareto front berechnen
%
% Punkt mit maximalem Revenue entlang jeder Kurve festhalten
%
% --> vergleiche die maximalen Rev
%
% in der nähe vom besten Punkt simulieren
%
% dann fmincon --> finale Länge L^{PRO}
%

%

%% Teste of initial data for paretosearch
i2=i2+1; close all;clc;
L= .5:.1:5;
option_mesh = 1e4; option_BVP = 1e-6; option_data = .3;    
% initial guess
% rng default
I1 = L(i2).*ones(1,70);
I1(1)
I2 = 52.*ones(1,70)+18.*rand(1,70);
I4 = 18.*ones(1,70)+2.*rand(1,70);
I3 = I4-.9.*rand(1,70)- .4.*ones(1,70);
I5 = 3.*ones(1,70)+ 2*rand(1,70);
%
X_init=[I1;I2;I3;I4;I5];
Y_init=zeros(5,70);
% test simulation
parfor i=1:70
    Y_init(:,i)=fun_1(X_init(:,i),option_data,'sol',option_mesh,option_BVP);
end
%% plots zum vergleich
load("DATA_case1.mat")
load("DATA_case2.mat")
f=figure(1); f.Position = [1282.3 0746.3 1277.3 0614.7];
scatter(Y1(1,:),Y1(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
scatter(Y2(1,:),Y2(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
scatter3(Y_init(1,:),Y_init(2,:),1:1:70,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','g'); hold on
%
grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
legend('Case1: \epsilon-constraint method', 'Case2: paretosearch',['Now using L^{PRO} = '  ,num2str(L(i2))],'Location', 'Northeast');
ylabel('FW [m^3/h]','FontSize',16);xlabel('SEC_{net} [kWh/m^3]','FontSize',16);
beep 

%% extract some data
load('DATA_Initial.mat')
I=[59 25 24 22 61 62 13 63 40 9];
Z=[];
for i=1:10
Z=[Z X_init(:,i)];
end
XL_25=Z';
save DATA_Initial.mat XL_05 XL_06 XL_07 XL_08 XL_09 XL_10 XL_11 XL_12 XL_13 XL_14 ...
                      XL_15 XL_16 XL_17 XL_18 XL_19 XL_20 XL_21 XL_22 XL_23 XL_24 XL_25
clc,

%% calculate all Pareto fronts:
%
close all;clc;
load DATA_Initial.mat
L=linspace(.5,2.5,21);
for j=16:21
    A= [0 0 1 -1 0]; b= -.01; Aeq=[]; beq=[]; lb = [L(j);30;1.01;1.01;1]; ub = [L(j);70;20;20;5];
    option_mesh = 1e4; option_BVP = 1e-6; option_data = .3;
    rng default
        if j==1; X_init=XL_05; end
        if j==2; X_init=XL_06; end
        if j==3; X_init=XL_07; end
        if j==4; X_init=XL_08; end
        if j==5; X_init=XL_09; end
        if j==6; X_init=XL_10; end
        if j==7; X_init=XL_11; end
        if j==8; X_init=XL_12; end
        if j==9; X_init=XL_13; end
        if j==10; X_init=XL_14; end
        if j==11; X_init=XL_15; end
        if j==12; X_init=XL_16; end
        if j==13; X_init=XL_17; end
        if j==14; X_init=XL_18; end
        if j==15; X_init=XL_19; end
        if j==16; X_init=XL_20; end
        if j==17; X_init=XL_21; end
        if j==18; X_init=XL_22; end
        if j==19; X_init=XL_23; end
        if j==20; X_init=XL_24; end
        if j==21; X_init=XL_25; end
    X_init=repmat(X_init, 20, 1); %200 x 5 double
    options = optimoptions('paretosearch','ParetoSetSize',200, 'InitialPoints',X_init,'Display','off', 'MaxFunctionEvaluations',10000, 'ParetoSetChangeTolerance',1e-8,'UseParallel', true);
    X = paretosearch(@(x)fun_1(x,option_data,'Pareto',option_mesh,option_BVP),5,A,b,Aeq,beq,lb,ub,@(x)nonlcon(x,'default'),options); 
    %
        if j==1
            parfor i=1:200
                    X_05(:,i)=X(i,:);
                    Y_05(:,i)=fun_1(X_05(:,i),option_data,'sol',option_mesh,option_BVP);
            end
            save DATA_05.mat X_05 Y_05
        end
        if j==2
            parfor i=1:200
                    X_06(:,i)=X(i,:);
                    Y_06(:,i)=fun_1(X_06(:,i),option_data,'sol',option_mesh,option_BVP);
            end
            save DATA_06.mat X_06 Y_06
        end
        if j==3
            parfor i=1:200
                    X_07(:,i)=X(i,:);
                    Y_07(:,i)=fun_1(X_07(:,i),option_data,'sol',option_mesh,option_BVP);
            end
            save DATA_07.mat X_07 Y_07
        end
        if j==4
            parfor i=1:200
                    X_08(:,i)=X(i,:);
                    Y_08(:,i)=fun_1(X_08(:,i),option_data,'sol',option_mesh,option_BVP);
            end
            save DATA_08.mat X_08 Y_08
        end
        if j==5
            parfor i=1:200
                    X_09(:,i)=X(i,:);
                    Y_09(:,i)=fun_1(X_09(:,i),option_data,'sol',option_mesh,option_BVP);
            end
            save DATA_09.mat X_09 Y_09
        end
        if j==6
            parfor i=1:200
                    X_10(:,i)=X(i,:);
                    Y_10(:,i)=fun_1(X_10(:,i),option_data,'sol',option_mesh,option_BVP);
            end
            save DATA_10.mat X_10 Y_10
        end
        if j==7
            parfor i=1:200
                    X_11(:,i)=X(i,:);
                    Y_11(:,i)=fun_1(X_11(:,i),option_data,'sol',option_mesh,option_BVP);
            end
            save DATA_11.mat X_11 Y_11
        end
        if j==8
            parfor i=1:200
                    X_12(:,i)=X(i,:);
                    Y_12(:,i)=fun_1(X_12(:,i),option_data,'sol',option_mesh,option_BVP);
            end
            save DATA_12.mat X_12 Y_12
        end
        if j==9
            parfor i=1:200
                    X_13(:,i)=X(i,:);
                    Y_13(:,i)=fun_1(X_13(:,i),option_data,'sol',option_mesh,option_BVP);
            end
            save DATA_13.mat X_13 Y_13
        end
        if j==10
            parfor i=1:200
                    X_14(:,i)=X(i,:);
                    Y_14(:,i)=fun_1(X_14(:,i),option_data,'sol',option_mesh,option_BVP);
            end
            save DATA_14.mat X_14 Y_14
        end
        if j==11
            parfor i=1:200
                    X_15(:,i)=X(i,:);
                    Y_15(:,i)=fun_1(X_15(:,i),option_data,'sol',option_mesh,option_BVP);
            end
            save DATA_15.mat X_15 Y_15
        end
        if j==12
            parfor i=1:200
                    X_16(:,i)=X(i,:);
                    Y_16(:,i)=fun_1(X_16(:,i),option_data,'sol',option_mesh,option_BVP);
            end
            save DATA_16.mat X_16 Y_16
        end
        if j==13
            parfor i=1:200
                    X_17(:,i)=X(i,:);
                    Y_17(:,i)=fun_1(X_17(:,i),option_data,'sol',option_mesh,option_BVP);
            end
            save DATA_17.mat X_17 Y_17
        end
        if j==14
            parfor i=1:200
                    X_18(:,i)=X(i,:);
                    Y_18(:,i)=fun_1(X_18(:,i),option_data,'sol',option_mesh,option_BVP);
            end
            save DATA_18.mat X_18 Y_18
        end
        if j==15
            parfor i=1:200
                    X_19(:,i)=X(i,:);
                    Y_19(:,i)=fun_1(X_19(:,i),option_data,'sol',option_mesh,option_BVP);
            end
            save DATA_19.mat X_19 Y_19
        end
        if j==16
            parfor i=1:200
                    X_20(:,i)=X(i,:);
                    Y_20(:,i)=fun_1(X_20(:,i),option_data,'sol',option_mesh,option_BVP);
            end
            save DATA_20.mat X_20 Y_20
        end
        if j==17
            parfor i=1:200
                    X_21(:,i)=X(i,:);
                    Y_21(:,i)=fun_1(X_21(:,i),option_data,'sol',option_mesh,option_BVP);
            end
            save DATA_21.mat X_21 Y_21
        end
        if j==18
            parfor i=1:200
                    X_22(:,i)=X(i,:);
                    Y_22(:,i)=fun_1(X_22(:,i),option_data,'sol',option_mesh,option_BVP);
            end
            save DATA_22.mat X_22 Y_22
        end
        if j==19
            parfor i=1:200
                    X_23(:,i)=X(i,:);
                    Y_23(:,i)=fun_1(X_23(:,i),option_data,'sol',option_mesh,option_BVP);
            end
            save DATA_23.mat X_23 Y_23
        end
        if j==20
            parfor i=1:200
                    X_24(:,i)=X(i,:);
                    Y_24(:,i)=fun_1(X_24(:,i),option_data,'sol',option_mesh,option_BVP);
            end
            save DATA_24.mat X_24 Y_24
        end
        if j==21
            parfor i=1:200
                    X_25(:,i)=X(i,:);
                    Y_25(:,i)=fun_1(X_25(:,i),option_data,'sol',option_mesh,option_BVP);
            end
            save DATA_25.mat X_25 Y_25
        end
        clc,
        j,
end

%% interactive scatter plot
close all;
scatter_switcher()

%% old plot
% from .. > Matlab > Hybrid Pareto front > Paretosearch Skipt
close all
fig = figure(1);
fig.Position=[1200 500 800 500];
load('DATA_paretosearch.mat') % old blue and red
load('Data_neu.mat') % old orange
scatter(-Y_paretosearch_2(:,1),-Y_paretosearch_2(:,2),'MarkerEdgeColor',[.3 0 0],'MarkerFaceColor','m'); hold on
scatter(-Y_paretosearch_1(:,1),-Y_paretosearch_1(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','c'); hold on
%scatter(-Y_paretosearch_3(:,1),-Y_paretosearch_3(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.9290 0.6940 0.1250]); hold on
scatter(-Y_neu(:,1),-Y_neu(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.9290 0.6940 0.1250]); hold on
%%%%%
grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
lgd=legend('Single SWRO unit','SWRO with ERD', 'Hybrid','Location','SouthWest');ylabel('FW [m^3/h]');xlabel('SEC_{net} [kWh/m^3]');title('Pareto front (using paretosearch)');
xlim([-5.501 -.74]); ylim([0 0.48]);view(2);

%% calculate f_fun(X:_neu)
X_neu=[X_neu(:,5) X_neu(:,1:4)]; % need to change positions
Y_Pareto=zeros(5,200);
parfor i=1:20
   X_Pareto(:,i)=X_neu(i,:);
   Y_Pareto(:,i)=fun_1(X_Pareto(:,i),.3,'sol',1e4,1e-4); %ACHTUNG not 1e-6 !!!!
   i,
end
%scatter(-Y_neu(:,1),-Y_neu(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.9290 0.6940 0.1250]); hold on

scatter(Y_Pareto(1,:),Y_Pareto(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','y'); hold on

% With exact same setup:
% 1e4; 1e-4; and KK=71.3; KD=71.3; KF=71.3;
% --> i can reproduce everything
%
% now try with
% 1e4; 1e-6; and KK=71.3; KD=71.3; KF=71.3;
% --> 103 NaN's in Y_Pareto (top and bottom are still ok, mostly missing in middle of the curve)
%

% now try with
% 1e4; 1e-4; and KK=7.13e2; KD=71.3; KF=71.3;
% --> 68 NaN (still pretty much the same)
%

% now try with (only 26 points)
% 1e4; 1e-4; and KK=7.13e2; KD=1/KK; KF=1/KK;
% --> now we see shift to the left (9 NaN's)
%

%sum(isnan(Y_Pareto(1,:)))







%% plot of best Pareto front

load DATA_case1.mat Y1
load DATA_case2.mat Y2
scatter(Y1(1,:),Y1(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
scatter(Y2(1,:),Y2(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on


%
% change the fun_1 such that i can reproduce the old pareto front
%
% find data X and Y data from old front
%
% calculate Y_new using fun_1
%
% check that Y_new == Y
%
% then run skript for one areto front again (L=1.5m)
%
% then run skript for all the others again
%
%
%



%% optimales revenue
%
% bester Revenue von den paretofronten ist ....
% ... mit L=
%
% annahme stetigkeit: anschließend optimierung innerhalb I=[ . , . ]
%
% hopefully optimal L
%
% Dann Pareto front von finaler Länge bestimmen
%
% nochmal bestes revenue mit neuer Pareto front vergleichen
%
% Dann in anderes Skript Überarbeiten
%
% __> PAPER FERTIG MACHEN !!!!!!!!!!!!!!!!!!!!!!!!!!

%%
%%

















close all;
scatter_switcher()
 function scatter_switcher()   
    % Create a figure and axis
    fig = figure('Name', 'Scatter Plot Switcher', 'NumberTitle', 'off');
    fig.Position=[1200 500 800 500];
    grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
    ax = axes('Parent', fig);
    % Create a slider
    slider = uicontrol('Style', 'slider', ...
                       'Min', 1, 'Max', 15, 'Value', 1, ...
                       'Units', 'normalized', ...
                       'Position', [0.15, 0.01, 0.7, 0.05]);
    % Add a listener to the slider
    addlistener(slider, 'Value', 'PreSet', @(src, event) updatePlot(ax, round(get(slider, 'Value'))));
    % Initial plot
    updatePlot(ax, round(get(slider,'Value')));
end

function updatePlot(ax, plotType)
    %Color=colorcube(21);
    Color=summer(30);
    Color=flipud(Color());
    cla(ax);
    % Plot based on slider value
    if plotType == 1
        load DATA_case1.mat Y1
        load DATA_case2.mat Y2
        scatter(Y1(1,:),Y1(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
        scatter(Y2(1,:),Y2(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
        load DATA_05.mat X_05 Y_05
        scatter(Y_05(1,:),Y_05(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(1,:)); hold on
        grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
        legend('Case1: \epsilon-constraint method', 'Case2: \epsilon-constraint method', ['L^{PRO} = ',num2str(X_05(1,1)), 'm'],'Location', 'Northeast');
        title(ax, ['SWRO vs PRO length ratio:    1 : ',num2str(X_05(1,1)/4)]);
    end
    if plotType == 2
        load DATA_case1.mat Y1
        load DATA_case2.mat Y2
        scatter(Y1(1,:),Y1(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
        scatter(Y2(1,:),Y2(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
        load DATA_06.mat X_06 Y_06
        scatter(Y_06(1,:),Y_06(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(2,:)); hold on
        grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
        legend('Case1: \epsilon-constraint method', 'Case2: \epsilon-constraint method', ['L^{PRO} = ',num2str(X_06(1,1)), 'm'],'Location', 'Northeast');
        title(ax, ['SWRO vs PRO length ratio:    1 : ',num2str(X_06(1,1)/4)]);
    end
    if plotType == 3
        load DATA_case1.mat Y1
        load DATA_case2.mat Y2
        scatter(Y1(1,:),Y1(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
        scatter(Y2(1,:),Y2(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
        load DATA_07.mat X_07 Y_07
        scatter(Y_07(1,:),Y_07(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(3,:)); hold on
        grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
        legend('Case1: \epsilon-constraint method', 'Case2: \epsilon-constraint method', ['L^{PRO} = ',num2str(X_07(1,1)), 'm'],'Location', 'Northeast');
        title(ax, ['SWRO vs PRO length ratio:    1 : ',num2str(X_07(1,1)/4)]);
    end
    if plotType == 4
        load DATA_case1.mat Y1
        load DATA_case2.mat Y2
        scatter(Y1(1,:),Y1(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
        scatter(Y2(1,:),Y2(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
        load DATA_08.mat X_08 Y_08
        scatter(Y_08(1,:),Y_08(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(4,:)); hold on
        grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
        legend('Case1: \epsilon-constraint method', 'Case2: \epsilon-constraint method', ['L^{PRO} = ',num2str(X_08(1,1)), 'm'],'Location', 'Northeast');
        title(ax, ['SWRO vs PRO length ratio:    1 : ',num2str(X_08(1,1)/4)]);
    end
    if plotType == 5
        load DATA_case1.mat Y1
        load DATA_case2.mat Y2
        scatter(Y1(1,:),Y1(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
        scatter(Y2(1,:),Y2(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
        load DATA_09.mat X_09 Y_09
        scatter(Y_09(1,:),Y_09(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(5,:)); hold on
        grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
        legend('Case1: \epsilon-constraint method', 'Case2: \epsilon-constraint method', ['L^{PRO} = ',num2str(X_09(1,1)), 'm'],'Location', 'Northeast');
        title(ax, ['SWRO vs PRO length ratio:    1 : ',num2str(X_09(1,1)/4)]);
    end
    if plotType == 6
        load DATA_case1.mat Y1
        load DATA_case2.mat Y2
        scatter(Y1(1,:),Y1(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
        scatter(Y2(1,:),Y2(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
        load DATA_10.mat X_10 Y_10
        scatter(Y_10(1,:),Y_10(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(6,:)); hold on
        grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
        legend('Case1: \epsilon-constraint method', 'Case2: \epsilon-constraint method', ['L^{PRO} = ',num2str(X_10(1,1)), 'm'],'Location', 'Northeast');
        title(ax, ['SWRO vs PRO length ratio:    1 : ',num2str(X_10(1,1)/4)]);
    end
    if plotType == 7
        load DATA_case1.mat Y1
        load DATA_case2.mat Y2
        scatter(Y1(1,:),Y1(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
        scatter(Y2(1,:),Y2(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
        load DATA_11.mat X_11 Y_11
        scatter(Y_11(1,:),Y_11(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(7,:)); hold on
        grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
        legend('Case1: \epsilon-constraint method', 'Case2: \epsilon-constraint method', ['L^{PRO} = ',num2str(X_11(1,1)), 'm'],'Location', 'Northeast');
        title(ax, ['SWRO vs PRO length ratio:    1 : ',num2str(X_11(1,1)/4)]);
    end
    if plotType == 8
        load DATA_case1.mat Y1
        load DATA_case2.mat Y2
        scatter(Y1(1,:),Y1(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
        scatter(Y2(1,:),Y2(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
        load DATA_12.mat X_12 Y_12
        scatter(Y_12(1,:),Y_12(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(8,:)); hold on
        grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
        legend('Case1: \epsilon-constraint method', 'Case2: \epsilon-constraint method', ['L^{PRO} = ',num2str(X_12(1,1)), 'm'],'Location', 'Northeast');
        title(ax, ['SWRO vs PRO length ratio:    1 : ',num2str(X_12(1,1)/4)]);
    end
    if plotType == 9
        load DATA_case1.mat Y1
        load DATA_case2.mat Y2
        scatter(Y1(1,:),Y1(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
        scatter(Y2(1,:),Y2(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
        load DATA_13.mat X_13 Y_13
        scatter(Y_13(1,:),Y_13(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(9,:)); hold on
        grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
        legend('Case1: \epsilon-constraint method', 'Case2: \epsilon-constraint method', ['L^{PRO} = ',num2str(X_13(1,1)), 'm'],'Location', 'Northeast');
        title(ax, ['SWRO vs PRO length ratio:    1 : ',num2str(X_13(1,1)/4)]);
    end
    if plotType == 10
        load DATA_case1.mat Y1
        load DATA_case2.mat Y2
        scatter(Y1(1,:),Y1(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
        scatter(Y2(1,:),Y2(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
        load DATA_14.mat X_14 Y_14
        scatter(Y_14(1,:),Y_14(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(10,:)); hold on
        grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
        legend('Case1: \epsilon-constraint method', 'Case2: \epsilon-constraint method', ['L^{PRO} = ',num2str(X_14(1,1)), 'm'],'Location', 'Northeast');
        title(ax, ['SWRO vs PRO length ratio:    1 : ',num2str(X_14(1,1)/4)]);
    end
    if plotType == 11
        load DATA_case1.mat Y1
        load DATA_case2.mat Y2
        scatter(Y1(1,:),Y1(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
        scatter(Y2(1,:),Y2(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
        load DATA_15.mat X_15 Y_15
        scatter(Y_15(1,:),Y_15(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(11,:)); hold on
        grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
        legend('Case1: \epsilon-constraint method', 'Case2: \epsilon-constraint method', ['L^{PRO} = ',num2str(X_15(1,1)), 'm'],'Location', 'Northeast');
        title(ax, ['SWRO vs PRO length ratio:    1 : ',num2str(X_15(1,1)/4)]);
    end
    if plotType == 12
        load DATA_case1.mat Y1
        load DATA_case2.mat Y2
        scatter(Y1(1,:),Y1(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
        scatter(Y2(1,:),Y2(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
        load DATA_16.mat X_16 Y_16
        scatter(Y_16(1,:),Y_16(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(12,:)); hold on
        grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
        legend('Case1: \epsilon-constraint method', 'Case2: \epsilon-constraint method', ['L^{PRO} = ',num2str(X_16(1,1)), 'm'],'Location', 'Northeast');
        title(ax, ['SWRO vs PRO length ratio:    1 : ',num2str(X_16(1,1)/4)]);
    end
    if plotType == 13
        load DATA_case1.mat Y1
        load DATA_case2.mat Y2
        scatter(Y1(1,:),Y1(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
        scatter(Y2(1,:),Y2(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
        load DATA_17.mat X_17 Y_17
        scatter(Y_17(1,:),Y_17(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(13,:)); hold on
        grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
        legend('Case1: \epsilon-constraint method', 'Case2: \epsilon-constraint method', ['L^{PRO} = ',num2str(X_17(1,1)), 'm'],'Location', 'Northeast');
        title(ax, ['SWRO vs PRO length ratio:    1 : ',num2str(X_17(1,1)/4)]);
    end
    if plotType == 14
        load DATA_case1.mat Y1
        load DATA_case2.mat Y2
        scatter(Y1(1,:),Y1(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
        scatter(Y2(1,:),Y2(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
        load DATA_18.mat X_18 Y_18
        scatter(Y_18(1,:),Y_18(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(14,:)); hold on
        grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
        legend('Case1: \epsilon-constraint method', 'Case2: \epsilon-constraint method', ['L^{PRO} = ',num2str(X_18(1,1)), 'm'],'Location', 'Northeast');
        title(ax, ['SWRO vs PRO length ratio:    1 : ',num2str(X_18(1,1)/4)]);
    end
    if plotType == 15
        load DATA_case1.mat Y1
        load DATA_case2.mat Y2
        scatter(Y1(1,:),Y1(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
        scatter(Y2(1,:),Y2(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
        load DATA_19.mat X_19 Y_19
        scatter(Y_19(1,:),Y_19(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(15,:)); hold on
        grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
        legend('Case1: \epsilon-constraint method', 'Case2: \epsilon-constraint method', ['L^{PRO} = ',num2str(X_19(1,1)), 'm'],'Location', 'Northeast');
        title(ax, ['SWRO vs PRO length ratio:    1 : ',num2str(X_19(1,1)/4)]);
    end
    if 1==0
    if plotType == 16
        load DATA_case1.mat Y1
        load DATA_case2.mat Y2
        scatter(Y1(1,:),Y1(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
        scatter(Y2(1,:),Y2(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
        load DATA_20.mat X_20 Y_20
        scatter(Y_20(1,:),Y_20(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(16,:)); hold on
        grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
        legend('Case1: \epsilon-constraint method', 'Case2: \epsilon-constraint method', ['L^{PRO} = ',num2str(X_20(1,1)), 'm'],'Location', 'Northeast');
        title(ax, ['SWRO vs PRO length ratio:    1 : ',num2str(X_20(1,1)/4)]);
    end
    if plotType == 17
        load DATA_case1.mat Y1
        load DATA_case2.mat Y2
        scatter(Y1(1,:),Y1(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
        scatter(Y2(1,:),Y2(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
        load DATA_21.mat X_21 Y_21
        scatter(Y_21(1,:),Y_21(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(17,:)); hold on
        grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
        legend('Case1: \epsilon-constraint method', 'Case2: \epsilon-constraint method', ['L^{PRO} = ',num2str(X_21(1,1)), 'm'],'Location', 'Northeast');
        title(ax, ['SWRO vs PRO length ratio:    1 : ',num2str(X_21(1,1)/4)]);
    end
    if plotType == 18
        load DATA_case1.mat Y1
        load DATA_case2.mat Y2
        scatter(Y1(1,:),Y1(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
        scatter(Y2(1,:),Y2(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
        load DATA_22.mat X_22 Y_22
        scatter(Y_22(1,:),Y_22(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(18,:)); hold on
        grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
        legend('Case1: \epsilon-constraint method', 'Case2: \epsilon-constraint method', ['L^{PRO} = ',num2str(X_22(1,1)), 'm'],'Location', 'Northeast');
        title(ax, ['SWRO vs PRO length ratio:    1 : ',num2str(X_22(1,1)/4)]);
    end
    if plotType == 19
        load DATA_case1.mat Y1
        load DATA_case2.mat Y2
        scatter(Y1(1,:),Y1(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
        scatter(Y2(1,:),Y2(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
        load DATA_23.mat X_23 Y_23
        scatter(Y_23(1,:),Y_23(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(19,:)); hold on
        grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
        legend('Case1: \epsilon-constraint method', 'Case2: \epsilon-constraint method', ['L^{PRO} = ',num2str(X_23(1,1)), 'm'],'Location', 'Northeast');
        title(ax, ['SWRO vs PRO length ratio:    1 : ',num2str(X_23(1,1)/4)]);
    end
    if plotType == 20
        load DATA_case1.mat Y1
        load DATA_case2.mat Y2
        scatter(Y1(1,:),Y1(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
        scatter(Y2(1,:),Y2(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
        load DATA_24.mat X_24 Y_24
        scatter(Y_24(1,:),Y_24(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(20,:)); hold on
        grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
        legend('Case1: \epsilon-constraint method', 'Case2: \epsilon-constraint method', ['L^{PRO} = ',num2str(X_24(1,1)), 'm'],'Location', 'Northeast');
        title(ax, ['SWRO vs PRO length ratio:    1 : ',num2str(X_24(1,1)/4)]);
    end
    if plotType == 21
        load DATA_case1.mat Y1
        load DATA_case2.mat Y2
        scatter(Y1(1,:),Y1(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
        scatter(Y2(1,:),Y2(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
        load DATA_25.mat X_25 Y_25
        scatter(Y_25(1,:),Y_25(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(21,:)); hold on
        grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
        legend('Case1: \epsilon-constraint method', 'Case2: \epsilon-constraint method', ['L^{PRO} = ',num2str(X_25(1,1)), 'm'],'Location', 'Northeast');
        title(ax, ['SWRO vs PRO length ratio:    1 : ',num2str(X_25(1,1)/4)]);
    end    
    end
end
