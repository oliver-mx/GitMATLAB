%% Basic simulations for singleSWRO/SWRO+ERD
% press [ctrl]+[enter] to run code sections
addpath('Input_DATA','Scaled_model','Unscaled_model','Output_DATA')

%% test random data
rng default
clear all; close all; clc,
n=10000;
X0=zeros(2,n); % <-- only SWRO
for i=1:n
    X0(1,i)= 31+39*rand();
    X0(2,i)= X0(1,i) - .1 - 3.3*rand();
end

%% Simulations
Y1=zeros(n,28);Y2=zeros(n,28);
startTime=datetime("now");
parfor i=1:n
        Y1(i,:)=fun_unscaled([0; X0(:,i)],.1,'sol',1e4,1e-4);
        Y2(i,:)=fun_unscaled([0; X0(:,i)],.2,'sol',1e4,1e-4);
end
endTime=datetime("now");
time = endTime - startTime; %
save Output_DATA/Test2_Output X0 Y1 Y2 n time

system('git add .'); system('git commit -m "Neue Simulation"');system('git push https://github.com/oliver-mx/GitMATLAB.git');

%% Plot Simulation results for singleSWRO/SWRO+ERD/HybridI/HybridII
load Output_DATA/Basic_Test_Output
close all; clc;
ev(Y1(1,:)); userNumber = input('Please enter a number of the set {1,2,3,...,18}: ');
f=figure(1); f.Position= [619 609.6667 1.5473e+03 648.0000];tiledlayout(1,2);
% same colorbar
cmap = jet(256); % Use the 'jet' colormap with 256 colors
sz=2.5; % size of points
%
nexttile
scatter3(X0(1,:), 10*(X0(1,:)-X0(2,:)),Y1(:,userNumber),sz,Y1(:,userNumber),'filled');hold on; title('simple SWRO'); xlabel('P_{d;0}^{RO} [bar]');axis equal;view(2); grid on;colormap(cmap);colorbar
xlim([30.7 70]);ylim([1.1 33]);yticks(10:10:30);yticklabels({'\DeltaP^{RO} = 1','\DeltaP^{RO} = 2','\DeltaP^{RO} = 3'})
nexttile
scatter3(X0(1,:), 10*(X0(1,:)-X0(2,:)),Y2(:,userNumber),sz,Y2(:,userNumber),'filled');hold on; title('SWRO with ERD'); xlabel('P_{d;0}^{RO} [bar]'); axis equal;view(2); grid on;colormap(cmap);colorbar
xlim([30.7 70]);ylim([1.1 33]);yticks(10:10:30);yticklabels({'\DeltaP^{RO} = 1','\DeltaP^{RO} = 2','\DeltaP^{RO} = 3'})

%% find good PRO operating conditions



%% Expected Pareto Fronts
load Output_DATA/Basic_Test_Output
close all; clc;
f=figure(1); f.Position= [12.3333 685.6667 2546 648];tiledlayout(1,4);
%for i=1:10000 %evtl kick out wrong data
%    if X0(1,i)-X0(2,i)>4 || X0(1,i)-X0(2,i)<.1
%        Y1(i,1)=-1;Y1(i,2)=-1;Y2(i,1)=-1;Y2(i,2)=-1;Y3(i,1)=-1;Y3(i,2)=-1;Y4(i,1)=-1;Y4(i,2)=-1;Y5(i,1)=-1;Y5(i,2)=-1;Y6(i,1)=-1;Y6(i,2)=-1;
%    end
%end
nexttile
scatter(Y1(:,1),Y1(:,2),'r','filled');hold on; title('simple SWRO'); xlabel('SEC_{net} [kWh/m^3]'); ylabel('FW [m^3/h]');
xlim([-10 0]);ylim([0 2.2]);
nexttile
scatter(Y2(:,1),Y2(:,2),'b','filled');hold on; title('SWRO with ERD'); xlabel('SEC_{net} [kWh/m^3]'); ylabel('FW [m^3/h]');
xlim([-10 0]);ylim([0 2.2]);
nexttile
scatter(Y3(:,1),Y3(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.8500 0.3250 0.0980]);hold on; title('Hybrid I'); xlabel('SEC_{net} [kWh/m^3]'); ylabel('FW [m^3/h]');
xlim([-10 0]);ylim([0 2.2]);
nexttile
scatter(Y4(:,1),Y4(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.9290 0.6940 0.1250]);hold on; title('Hybrid II'); xlabel('SEC_{net} [kWh/m^3]'); ylabel('FW [m^3/h]');
xlim([-10 0]);ylim([0 2.2]);

%evtl one figure that combines them
f=figure(2);f.Position= [823.6667 167.6667 989.3333 483.3333];
scatter(Y1(:,1),Y1(:,2),'r','filled');hold on;xlabel('SEC_{net} [kWh/m^3]'); ylabel('FW [m^3/h]');
scatter(Y2(:,1),Y2(:,2),'b','filled');hold on; title('SWRO with ERD'); xlabel('SEC_{net} [kWh/m^3]'); ylabel('FW [m^3/h]');
scatter(Y3(:,1),Y3(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.8500 0.3250 0.0980]);hold on; title('Hybrid I'); xlabel('SEC_{net} [kWh/m^3]'); ylabel('FW [m^3/h]');
scatter(Y4(:,1),Y4(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.9290 0.6940 0.1250]);hold on; title('Hybrid II'); xlabel('SEC_{net} [kWh/m^3]'); ylabel('FW [m^3/h]');
xlim([-10 0]);ylim([0 2.2]);
%scatter(-1.7763,0.56987,'yellow','filled')




















