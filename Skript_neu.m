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
close all;clc;
%i2=0; rng default % bei Initialisierung
i2=i2+1; 
L= .5:.1:5;
option_mesh = 1e3; option_BVP = 1e-4; option_data = .3;  
k=14*8;
% initial guess
I1 = L(i2).*ones(1,k);
I1(1)
I2 = 60.*ones(1,k)+10.*rand(1,k);
I4 = 16.*ones(1,k)+4.*rand(1,k);
I3 = I4-.9.*rand(1,k)- .2.*ones(1,k);
I5 = 1.2.*ones(1,k)+ 2*rand(1,k);
%
X_init=[I1;I2;I3;I4;I5];
Y_init=zeros(5,k);
% test simulation
parfor i=1:k
    Y_init(:,i)=fun_1(X_init(:,i),option_data,'sol',option_mesh,option_BVP);
end
% plots zum vergleich
load("DATA_case1.mat")
load("DATA_case2.mat")
f=figure(1); f.Position = [1282.3 0746.3 1277.3 0614.7];
scatter(Y1(1,:),Y1(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
scatter(Y2(1,:),Y2(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
scatter3(Y_init(1,:),Y_init(2,:),1:1:k,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','g'); hold on
%
grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
legend('Case1: \epsilon-constraint method', 'Case2: paretosearch',['Now using L^{PRO} = '  ,num2str(L(i2))],'Location', 'Northeast');
ylabel('FW [m^3/h]','FontSize',16);xlabel('SEC_{net} [kWh/m^3]','FontSize',16);
beep 

%% extract some data
load('DATA_Initial.mat')
I=[107 9 4 66];
a=length(I);
while a<10
    I=[I I(end)];
    a=length(I);
end
Z=[];% keep it empty :)
for i=1:10
Z=[Z X_init(:,i)];
end
XL_25=Z';
save DATA_Initial.mat XL_05 XL_06 XL_07 XL_08 XL_09 XL_10 XL_11 XL_12 XL_13 XL_14 ...
                      XL_15 XL_16 XL_17 XL_18 XL_19 XL_20 XL_21 XL_22 XL_23 XL_24 XL_25
clc,


%% remove some ausreißer
load('DATA_Initial.mat')
load('DATA_15.mat')
I=[42 6 25 160 85 28 30 12];
I=sort(I);
Z=X_15;
a=1;
while a<length(I)+1
    index=I(a);
    Z(:,index)=XL_15(:,1);
    a=a+1;
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
for j=22:22   
    %
    A= []; b=[]; Aeq=[]; beq=[]; lb = [L(j);30;1.01;1.01;1]; ub = [L(j);70;20;20;5];
    option_mesh = 1e3; option_BVP = 1e-4; option_data = .3;
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
    options = optimoptions('paretosearch','ParetoSetSize',200, 'InitialPoints',X_init','Display','iter', 'MaxFunctionEvaluations',10000, 'ParetoSetChangeTolerance',1e-7,'UseParallel', true);
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
        display('finished with:')
        j,
end
%
clc

%
%system('git status');
%system('git add .');
%system('git commit -m "fmincon sqp"');
%system('git push https://github.com/oliver-mx/GitMATLAB.git');
%system('shutdown /s /t 30');

%% interactive Pareto front plot
close all;
scatter_pareto()

%% interactive revenue plot
close all;
scatter_revenue()

%% try 3d surf plot ???
close all;
surf_pareto() % evtl mark optimal revenue point

%% plto only selected curves
close all;
scatter_select([.6 1.4 2.1 2.5])

%%
clc,
Z=[159 108 86 92];
for i=1:length(Z)
disp(X_25(:,Z(i))')
end
%%
k=1200;
option_mesh = 1e3; option_BVP = 1e-4; option_data = .3;
I1=2.5.*ones(1,k);
I2=70.*ones(1,k); 
I3=10.*ones(1,k)+5.*rand(1,k);
I4=20.*ones(1,k); 
I5=4..*ones(1,k)+rand(1,k); 
X_init=[I1;I2;I3;I4;I5];
parfor i=1:k
Y_sol(:,i)=fun_1(X_init(:,i),option_data,'sol',option_mesh,option_BVP);
end
%%
scatter3(Y_sol(1,:),Y_sol(2,:),1:1:k,'r')
%

%
%%
clc,
X_25(:,1:8);
equal_to_first=all(X_25 == X_25(:,1), 1),
find(equal_to_first)

%% replace 



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

% calculate f_fun(X:_neu)
X_neu=[X_neu(:,5) X_neu(:,1:4)]; % need to change positions
Y_Pareto=zeros(5,40);
parfor i=1:40
   X_Pareto(:,i)=X_neu(i,:);
   Y_Pareto(:,i)=fun_1(X_Pareto(:,i),.3,'sol',1e3,1e-4); %ACHTUNG not 1e4 1e-6 !!!!
   i,
end
scatter(Y_Pareto(1,1:40),Y_Pareto(2,1:40),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','y'); hold on

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

% now try
% 1e4; 1e-4; and KK=7.13e2; KD=7.13e-2; KF=7.13e-2;
% --> 
%

% next simulation
% 1e4; 1e-4; and KK=7.13e2; KD=7.13e-2; KF=7.13e-2;
%
%

%sum(isnan(Y_Pareto(1,:)))

%% plot of best Pareto front
load DATA_case1.mat Y1
load DATA_case2.mat Y2
scatter(Y1(1,:),Y1(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
scatter(Y2(1,:),Y2(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on

%%
load Data_neu.mat
scatter3(-Y_neu(:,1), -Y_neu(:,2),1:1:200, 'c')

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






































%% all functions of the script:
close all
surf_pareto()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function scatter_revenue()   
    % Create a figure and axis
    fig = figure('Name', 'Interactive revenue plot', 'NumberTitle', 'off');
    fig.Position=[1200 500 800 500];
    ax = axes('Parent', fig);
    set(ax, 'Position', [0.1, 0.25, 0.8, 0.6]);
    grid on; xlim([0 3]); ylim([-.1 1.5]);view(2);
    xlabel('please move the slider')
    hLabel = uicontrol('Style', 'text', ...
                   'Position', [10, 45, 110, 20], ... % [left, bottom, width, height]
                   'String', 'Water price:', ...
                   'FontSize', 10, ...
                   'FontWeight', 'bold', ...
                   'BackgroundColor', get(0, 'defaultuicontrolbackgroundcolor'));
    h2Label = uicontrol('Style', 'text', ...
                   'Position', [10, 8, 110, 20], ... % [left, bottom, width, height]
                   'String', 'Electricity price:', ...
                   'FontSize', 10, ...
                   'FontWeight', 'bold', ...
                   'BackgroundColor', get(0, 'defaultuicontrolbackgroundcolor'));
    % Create a slider
    slider1 = uicontrol('Style', 'slider', ...
                       'Min', 0, 'Max', 4, 'Value', 2, ...
                       'Units', 'normalized', ...
                       'Position', [0.15, 0.08, 0.7, 0.05]);
    slider2 = uicontrol('Style', 'slider', ...
                       'Min', 0, 'Max', 1, 'Value', .3, ...
                       'Units', 'normalized', ...
                       'Position', [0.15, 0.01, 0.7, 0.05]);
    % Add a listener to the slider
    addlistener(slider1, 'Value', 'PreSet', @(src, event) updatePlot2(ax, round((get(slider1, 'Value')-.01)*10)/10, round(get(slider2, 'Value')*10)/10));
    addlistener(slider2, 'Value', 'PreSet', @(src, event) updatePlot2(ax, round((get(slider1, 'Value')-.01)*10)/10, round(get(slider2, 'Value')*10)/10));  
    % Initial plot
    updatePlot2(ax, round((get(slider1, 'Value')-.01)*10)/10, round(get(slider2, 'Value')*10)/10);
end
function updatePlot2(ax, s1, s2)
    % Plot based on slider value
    Color=summer(30);
    Color=flipud(Color());
    cla(ax);
    a=zeros(1,21);
    %
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
    %
    %load Data_neu.mat Y_neu 
    %Y_15=-Y_neu';
    %
    a(1)=max(s1.*Y_05(2,:) + s2.*Y_05(2,:).*Y_05(1,:));
    scatter(X_05(1,1),a(1),'MarkerEdgeColor','none','MarkerFaceColor',Color(5,:));hold on
    %
    a(2)=max(s1.*Y_06(2,:) + s2.*Y_06(2,:).*Y_06(1,:));
    scatter(X_06(1,1),a(2),'MarkerEdgeColor','none','MarkerFaceColor',Color(6,:));hold on
    %    
    a(3)=max(s1.*Y_07(2,:) + s2.*Y_07(2,:).*Y_07(1,:));
    scatter(X_07(1,1),a(3),'MarkerEdgeColor','none','MarkerFaceColor',Color(7,:));hold on
    %   
    a(4)=max(s1.*Y_08(2,:) + s2.*Y_08(2,:).*Y_08(1,:));
    scatter(X_08(1,1),a(4),'MarkerEdgeColor','none','MarkerFaceColor',Color(8,:));hold on
    %
    a(5)=max(s1.*Y_09(2,:) + s2.*Y_09(2,:).*Y_09(1,:));
    scatter(X_09(1,1),a(5),'MarkerEdgeColor','none','MarkerFaceColor',Color(9,:));hold on
    % 
    a(6)=max(s1.*Y_10(2,:) + s2.*Y_10(2,:).*Y_10(1,:));
    scatter(X_10(1,1),a(6),'MarkerEdgeColor','none','MarkerFaceColor',Color(10,:));hold on
    %
    a(7)=max(s1.*Y_11(2,:) + s2.*Y_11(2,:).*Y_11(1,:));
    scatter(X_11(1,1),a(7),'MarkerEdgeColor','none','MarkerFaceColor',Color(11,:));hold on
    %
    a(8)=max(s1.*Y_12(2,:) + s2.*Y_12(2,:).*Y_12(1,:));
    scatter(X_12(1,1),a(8),'MarkerEdgeColor','none','MarkerFaceColor',Color(12,:));hold on
    %
    a(9)=max(s1.*Y_13(2,:) + s2.*Y_13(2,:).*Y_13(1,:));
    scatter(X_13(1,1),a(9),'MarkerEdgeColor','none','MarkerFaceColor',Color(13,:));hold on
    %
    a(10)=max(s1.*Y_14(2,:) + s2.*Y_14(2,:).*Y_14(1,:));
    scatter(X_14(1,1),a(10),'MarkerEdgeColor','none','MarkerFaceColor',Color(14,:));hold on
    %
    a(11)=max(s1.*Y_15(2,:) + s2.*Y_15(2,:).*Y_15(1,:));
    scatter(X_15(1,1),a(11),'MarkerEdgeColor','none','MarkerFaceColor',Color(15,:));hold on
    %
    a(12)=max(s1.*Y_16(2,:) + s2.*Y_16(2,:).*Y_16(1,:));
    scatter(X_16(1,1),a(12),'MarkerEdgeColor','none','MarkerFaceColor',Color(16,:));hold on
    %
    a(13)=max(s1.*Y_17(2,:) + s2.*Y_17(2,:).*Y_17(1,:));
    scatter(X_17(1,1),a(13),'MarkerEdgeColor','none','MarkerFaceColor',Color(17,:));hold on
    %
    a(14)=max(s1.*Y_18(2,:) + s2.*Y_18(2,:).*Y_18(1,:));
    scatter(X_18(1,1),a(14),'MarkerEdgeColor','none','MarkerFaceColor',Color(18,:));hold on
    %
    a(15)=max(s1.*Y_19(2,:) + s2.*Y_19(2,:).*Y_19(1,:));
    scatter(X_19(1,1),a(15),'MarkerEdgeColor','none','MarkerFaceColor',Color(19,:));hold on
    %
    a(16)=max(s1.*Y_20(2,:) + s2.*Y_20(2,:).*Y_20(1,:));
    scatter(X_20(1,1),a(16),'MarkerEdgeColor','none','MarkerFaceColor',Color(20,:));hold on
    %
    a(17)=max(s1.*Y_21(2,:) + s2.*Y_21(2,:).*Y_21(1,:));
    scatter(X_21(1,1),a(17),'MarkerEdgeColor','none','MarkerFaceColor',Color(21,:));hold on
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
    [b,c]=max(a);
    c=0.4+c*.1;
    % highlight max revenue: 
    plot([c,c], [0,b*.96], 'r');
    scatter(c, b, 'r');
    %
    [~,bb]=max(s1.*Y_22(2,:) + s2.*Y_22(2,:).*Y_22(1,:));
    bb,
    disp(X_22(:,bb)),
    %
    grid on; xlim([0 3]); ylim([-.1 1.5]);view(2);
    ylabel('Revenue [$/h]')
    xlabel('L^{PRO} [m]')
    legend(['Max revenue = ',num2str(b),' $/h'], ['Highest revenue with L^{PRO} = ',num2str(c),' m'])
    % title(['water price = ',num2str(s1),' $/m^3,              electricity price = ', num2str(s2), ' $/kwh']);
    title({'Revenue comparison (L^{RO} = 4m):', ['with  p_w = ',num2str(s1),' $/m^3,   p_e = ', num2str(s2), ' $/kwh']});
    h3Label = uicontrol('Style', 'text', ...
                   'Position', [680, 45, 50, 20], ... % [left, bottom, width, height]
                   'String', num2str(s1), ...
                   'FontSize', 10, ...
                   'FontWeight', 'bold', ...
                   'BackgroundColor', get(0, 'defaultuicontrolbackgroundcolor'));
   h4Label = uicontrol('Style', 'text', ...
                   'Position', [680, 8, 50, 20], ... % [left, bottom, width, height]
                   'String', num2str(s2), ...
                   'FontSize', 10, ...
                   'FontWeight', 'bold', ...
                   'BackgroundColor', get(0, 'defaultuicontrolbackgroundcolor'));
end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function scatter_pareto()   
    % Create a figure and axis
    fig = figure('Name', 'Interactive Pareto fronts plot', 'NumberTitle', 'off');
    fig.Position=[1200 500 800 500];
    grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
    ax = axes('Parent', fig);
    % Create a slider
    slider = uicontrol('Style', 'slider', ...
                       'Min', 1, 'Max', 21, 'Value', 1, ...
                       'Units', 'normalized', ...
                       'Position', [0.15, 0.01, 0.7, 0.05]);
    % Add a listener to the slider
    addlistener(slider, 'Value', 'PreSet', @(src, event) updatePlot(ax, round(get(slider, 'Value'))));
    % Initial plot
    updatePlot(ax, round(get(slider,'Value')));
end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        scatter3(Y_05(1,:),Y_05(2,:),1:1:200,'MarkerEdgeColor','none','MarkerFaceColor',Color(4,:)); hold on
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
        scatter3(Y_06(1,:),Y_06(2,:),1:1:200,'MarkerEdgeColor','none','MarkerFaceColor',Color(5,:)); hold on
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
        scatter3(Y_07(1,:),Y_07(2,:),1:1:200,'MarkerEdgeColor','none','MarkerFaceColor',Color(6,:)); hold on
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
        scatter3(Y_08(1,:),Y_08(2,:),1:1:200,'MarkerEdgeColor','none','MarkerFaceColor',Color(7,:)); hold on
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
        scatter3(Y_09(1,:),Y_09(2,:),1:1:200,'MarkerEdgeColor','none','MarkerFaceColor',Color(8,:)); hold on
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
        scatter3(Y_10(1,:),Y_10(2,:),1:1:200,'MarkerEdgeColor','none','MarkerFaceColor',Color(9,:)); hold on
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
        scatter3(Y_11(1,:),Y_11(2,:),1:1:200,'MarkerEdgeColor','none','MarkerFaceColor',Color(10,:)); hold on
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
        scatter3(Y_12(1,:),Y_12(2,:),1:1:200,'MarkerEdgeColor','none','MarkerFaceColor',Color(11,:)); hold on
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
        scatter3(Y_13(1,:),Y_13(2,:),1:1:200,'MarkerEdgeColor','none','MarkerFaceColor',Color(12,:)); hold on
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
        scatter3(Y_14(1,:),Y_14(2,:),1:1:200,'MarkerEdgeColor','none','MarkerFaceColor',Color(13,:)); hold on
        grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
        legend('Case1: \epsilon-constraint method', 'Case2: \epsilon-constraint method', ['L^{PRO} = ',num2str(X_14(1,1)), 'm'],'Location', 'Northeast');
        title(ax, ['SWRO vs PRO length ratio:    1 : ',num2str(X_14(1,1)/4)]);
    end
    if plotType == 11
        load DATA_case1.mat Y1
        load DATA_case2.mat Y2
        scatter(Y1(1,:),Y1(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
        scatter(Y2(1,:),Y2(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
        %load Data_neu.mat Y_neu 
        load DATA_15.mat X_15 Y_15
        %Y_15=-Y_neu'; %scatter3(-Y_neu(:,1), -Y_neu(:,2),1:1:200, 'c')
        scatter3(Y_15(1,:),Y_15(2,:),1:1:200,'MarkerEdgeColor','none','MarkerFaceColor',Color(14,:)); hold on
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
        scatter3(Y_16(1,:),Y_16(2,:),1:1:200,'MarkerEdgeColor','none','MarkerFaceColor',Color(15,:)); hold on
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
        scatter3(Y_17(1,:),Y_17(2,:),1:1:200,'MarkerEdgeColor','none','MarkerFaceColor',Color(16,:)); hold on
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
        scatter3(Y_18(1,:),Y_18(2,:),1:1:200,'MarkerEdgeColor','none','MarkerFaceColor',Color(17,:)); hold on
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
        scatter3(Y_19(1,:),Y_19(2,:),1:1:200,'MarkerEdgeColor','none','MarkerFaceColor',Color(18,:)); hold on
        grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
        legend('Case1: \epsilon-constraint method', 'Case2: \epsilon-constraint method', ['L^{PRO} = ',num2str(X_19(1,1)), 'm'],'Location', 'Northeast');
        title(ax, ['SWRO vs PRO length ratio:    1 : ',num2str(X_19(1,1)/4)]);
    end
    if plotType == 16
        load DATA_case1.mat Y1
        load DATA_case2.mat Y2
        scatter(Y1(1,:),Y1(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
        scatter(Y2(1,:),Y2(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
        load DATA_20.mat X_20 Y_20
        scatter3(Y_20(1,:),Y_20(2,:),1:1:200,'MarkerEdgeColor','none','MarkerFaceColor',Color(19,:)); hold on
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
        scatter3(Y_21(1,:),Y_21(2,:),1:1:200,'MarkerEdgeColor','none','MarkerFaceColor',Color(20,:)); hold on
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
        scatter3(Y_22(1,:),Y_22(2,:),1:1:200,'MarkerEdgeColor','none','MarkerFaceColor',Color(21,:)); hold on
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
        scatter3(Y_23(1,:),Y_23(2,:),1:1:200,'MarkerEdgeColor','none','MarkerFaceColor',Color(22,:)); hold on
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
        scatter3(Y_24(1,:),Y_24(2,:),1:1:200,'MarkerEdgeColor','none','MarkerFaceColor',Color(23,:)); hold on
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
        scatter3(Y_25(1,:),Y_25(2,:),1:1:200,'MarkerEdgeColor','none','MarkerFaceColor',Color(24,:)); hold on
        grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
        legend('Case1: \epsilon-constraint method', 'Case2: \epsilon-constraint method', ['L^{PRO} = ',num2str(X_25(1,1)), 'm'],'Location', 'Northeast');
        title(ax, ['SWRO vs PRO length ratio:    1 : ',num2str(X_25(1,1)/4)]);
    end 
    hLabel = uicontrol('Style', 'text', ...
                   'Position', [10, 8, 110, 20], ... % [left, bottom, width, height]
                   'String', '       L^{PRO}:', ...
                   'FontSize', 10, ...
                   'FontWeight', 'bold', ...
                   'BackgroundColor', get(0, 'defaultuicontrolbackgroundcolor'));
    h2Label = uicontrol('Style', 'text', ...
                   'Position', [680, 8, 50, 20], ... % [left, bottom, width, height]
                   'String', num2str(0.4+plotType*.1), ...
                   'FontSize', 10, ...
                   'FontWeight', 'bold', ...
                   'BackgroundColor', get(0, 'defaultuicontrolbackgroundcolor'));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function surf_pareto()
    % Create a figure and axis
    fig = figure('Name', 'Interactive revenue plot', 'NumberTitle', 'off');
    fig.Position=[1200 500 800 500];
    L = 0.5:0.1:2.5;
    A=repmat(L,200,1);
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
    B = [Y_05(1,:); Y_06(1,:); Y_07(1,:); Y_08(1,:); Y_09(1,:); Y_10(1,:); Y_11(1,:); Y_12(1,:); Y_13(1,:); Y_14(1,:); Y_15(1,:); Y_16(1,:); Y_17(1,:); Y_18(1,:); Y_19(1,:); Y_20(1,:); Y_21(1,:); Y_22(1,:); Y_23(1,:); Y_24(1,:); Y_25(1,:)];
    C = [Y_05(2,:); Y_06(2,:); Y_07(2,:); Y_08(2,:); Y_09(2,:); Y_10(2,:); Y_11(2,:); Y_12(2,:); Y_13(2,:); Y_14(2,:); Y_15(2,:); Y_16(2,:); Y_17(2,:); Y_18(2,:); Y_19(2,:); Y_20(2,:); Y_21(2,:); Y_22(2,:); Y_23(2,:); Y_24(2,:); Y_25(2,:)];    
    h1 = surf(A,-B',C'); hold on;
    C1 = C'; 
    set(h1, 'CData', C1, 'FaceColor', 'interp');
    colormap(parula); 
    %
    load DATA_case1.mat Y1
    load DATA_case2.mat Y2
    %
    surf(A,-repmat(Y1(1,:),21,1)',repmat(Y1(2,:),21,1)','FaceColor', 'red'); hold on;
    surf(A,-repmat(Y2(1,:),21,1)',repmat(Y2(2,:),21,1)','FaceColor', 'blue'); hold on;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function scatter_select(I)
    % Create a figure and axis
    fig = figure('Name', 'Interactive revenue plot', 'NumberTitle', 'off');
    fig.Position=[1200 500 800 500];
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
    load DATA_case1.mat Y1
    load DATA_case2.mat Y2
    scatter(Y1(1,:),Y1(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); hold on
    scatter(Y2(1,:),Y2(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b'); hold on
    if any(I==.5);scatter(Y_05(1,:),Y_05(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(1,:)); hold on; end
    if any(I==.6);scatter(Y_06(1,:),Y_06(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(2,:)); hold on; end
    if any(I==.7);scatter(Y_07(1,:),Y_07(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(3,:)); hold on; end
    if any(I==.8);scatter(Y_08(1,:),Y_08(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(4,:)); hold on; end
    if any(I==.9);scatter(Y_09(1,:),Y_09(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(5,:)); hold on; end
    if any(I==1);scatter(Y_10(1,:),Y_10(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(6,:)); hold on; end
    if any(I==1.1);scatter(Y_11(1,:),Y_11(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(7,:)); hold on; end
    if any(I==1.2);scatter(Y_12(1,:),Y_12(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(8,:)); hold on; end
    if any(I==1.3);scatter(Y_13(1,:),Y_13(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(9,:)); hold on; end
    if any(I==1.4);scatter(Y_14(1,:),Y_14(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(10,:)); hold on; end
    if any(I==1.5);scatter(Y_15(1,:),Y_15(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(11,:)); hold on; end
    if any(I==1.6);scatter(Y_16(1,:),Y_16(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(12,:)); hold on; end
    if any(I==1.7);scatter(Y_17(1,:),Y_17(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(13,:)); hold on; end
    if any(I==1.8);scatter(Y_18(1,:),Y_18(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(14,:)); hold on; end
    if any(I==1.9);scatter(Y_19(1,:),Y_19(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(15,:)); hold on; end
    if any(I==2);scatter(Y_20(1,:),Y_20(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(16,:)); hold on; end
    if any(I==2.1);scatter(Y_21(1,:),Y_21(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(17,:)); hold on; end
    if any(I==2.2);scatter(Y_22(1,:),Y_22(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(18,:)); hold on; end
    if any(I==2.3);scatter(Y_23(1,:),Y_23(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(19,:)); hold on; end
    if any(I==2.4);scatter(Y_24(1,:),Y_24(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(20,:)); hold on; end
    if any(I==2.5);scatter(Y_25(1,:),Y_25(2,:),'MarkerEdgeColor','none','MarkerFaceColor',Color(21,:)); hold on; end 
    grid on; xlim([-5.501 -.5]); ylim([0 0.48]);view(2);
    if length(I)==1;legend('Case1: L^{RO}= 4m', 'Case2: L^{RO}= 4m', ['Hybrid: L^{RO}= 4m, L^{PRO}= ',num2str(I(1)),'m'],'Location', 'Northeast');end
    if length(I)==2;legend('Case1: L^{RO}= 4m', 'Case2: L^{RO}= 4m', ['Hybrid: L^{RO}= 4m, L^{PRO}= ',num2str(I(1)),'m'],['Hybrid: L^{RO}= 4m, L^{PRO}= ',num2str(I(2)),'m'],'Location', 'SouthWest');end
    if length(I)==3;legend('Case1: L^{RO}= 4m', 'Case2: L^{RO}= 4m', ['Hybrid: L^{RO}= 4m, L^{PRO}= ',num2str(I(1)),'m'],['Hybrid: L^{RO}= 4m, L^{PRO}= ',num2str(I(2)),'m'],['Hybrid: L^{RO}= 4m, L^{PRO}= ',num2str(I(3)),'m'],'Location', 'SouthWest');end
    if length(I)==4;legend('Case1: L^{RO}= 4m', 'Case2: L^{RO}= 4m', ['Hybrid: L^{RO}= 4m, L^{PRO}= ',num2str(I(1)),'m'],['Hybrid: L^{RO}= 4m, L^{PRO}= ',num2str(I(2)),'m'],['Hybrid: L^{RO}= 4m, L^{PRO}= ',num2str(I(3)),'m'],['Hybrid: L^{RO}= 4m, L^{PRO}= ',num2str(I(4)),'m'],'Location', 'SouthWest');end
    if length(I)>4;legend('Case1: L^{RO}= 4m', 'Case2: L^{RO}= 4m', ['Hybrid: L^{RO}= 4m, L^{PRO}= ',num2str(I(1)),'m'],['Hybrid: L^{RO}= 4m, L^{PRO}= ',num2str(I(2)),'m'],['Hybrid: L^{RO}= 4m, L^{PRO}= ',num2str(I(3)),'m'],['Hybrid: L^{RO}= 4m, L^{PRO}= ',num2str(I(4)),'m'],['Hybrid: L^{RO}= 4m, L^{PRO}= ',num2str(I(5)),'m'],'Location', 'SouthWest');end
    title('Pareto fronts');
end
