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
load Output_DATA/Test2_Output
close all;
ev;
userNumber = input('Please enter a number of the set {1,2,3,...,18}: ');
f=figure(1); f.Position= [619 609.6667 1.5473e+03 648.0000];tiledlayout(1,2);
% same colorbar
cmap = jet(256); % Use the 'jet' colormap with 256 colors
sz=8.5; % size of points
%
nexttile
scatter3(X0(1,:), 10*(X0(1,:)-X0(2,:)),Y1(:,userNumber),sz,Y1(:,userNumber),'filled');hold on; title('simple SWRO'); xlabel('P_{d;0}^{RO} [bar]');axis equal;view(2); grid on;colormap(cmap);colorbar
xlim([30.7 70]);ylim([1.1 33]);yticks(10:10:30);yticklabels({'\DeltaP^{RO} = 1','\DeltaP^{RO} = 2','\DeltaP^{RO} = 3'})
nexttile
scatter3(X0(1,:), 10*(X0(1,:)-X0(2,:)),Y2(:,userNumber),sz,Y2(:,userNumber),'filled');hold on; title('SWRO with ERD'); xlabel('P_{d;0}^{RO} [bar]'); axis equal;view(2); grid on;colormap(cmap);colorbar
xlim([30.7 70]);ylim([1.1 33]);yticks(10:10:30);yticklabels({'\DeltaP^{RO} = 1','\DeltaP^{RO} = 2','\DeltaP^{RO} = 3'})

%% overview rejected brine (flow rate+concentration)
% .99*Y2(17) und Y2(18)
load Output_DATA/Test2_Output
J_brine=.99*Y2(:,17);close all;
C_brine=Y2(:,18); A=[J_brine C_brine];
scatter(C_brine, J_brine, 'b');hold on % overview of incomming brine with concentration:
f=figure(1); % :
title('possible brine outlet');xlabel('C_d(L) [%]');ylabel('0.99 \times J_d(L) [kg/sm]')
% select sample points
sample=[3836 4305 1668 7676 7872 9311 987 1257 7393 7338 770 8350 7964 3538 2559];
scatter(C_brine(sample), J_brine(sample), 'r', 'filled');hold on 

%% calculate best PRO operating conditions for smaple points
% for Pressures 5, 10, 15, 20 bar
test_pressures=[1.1 1.07 1.05 1.03 1.01 1.007 1.005 1.003 1.001 1.0007 1.0005 1.0003 1.0001 1.00007 1.00005 1.00003 1.00001 1.000007 1.000005 1.000003 1.000001 1.000007 1.000005 1.000003 1.000001];
P=zeros(4,length(sample)); PP=[5 10 15 20];
for p=1:0 %4
    clc; disp(p)
    for i=1:length(sample)
        output=zeros(28,length(test_pressures));
        for j=1:length(test_pressures)
        output(:,j)=fun_scaled([-J_brine(i); 0 ; PP(p) ; test_pressures(j); C_brine(i)],.3,'sol',1e3,1e-3);
        end
        [a,b]=max(output(9,:)); % we want max W_net
        P(p,i)=test_pressures(b); %save the best pressure 
    end
end
%save PRO_TEST.mat P

%% hybrid simulations
load Output_DATA/Test2_Output
fun_p(i)
startTime=datetime("now");
parfor i=1:n
        Y1(i,:)=fun_unscaled([0; X0(:,i)],.1,'sol',1e4,1e-4);
        Y2(i,:)=fun_unscaled([0; X0(:,i)],.2,'sol',1e4,1e-4);
end
endTime=datetime("now");
time = endTime - startTime; %
save Output_DATA/Test3_Output X0 Y1 Y2 n time

%%
clc,PP=[5 10 15 20];
k=1;
i=sample(k);p=1; 
out1=fun_unscaled([0; X0(:,i)],.2,'sol',1e4,1e-4);
ev(out1,[1 2 4 5 7])
[X0(:,i);PP(p);P(p,k)]

out2=fun_scaled([X0(:,i);PP(p);P(p,k)],.4,'sol',1e4,1e-3);
ev(out2,[1 2 4 5 7])















%% Expected Pareto Fronts
load Output_DATA/Test2_Output
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


















%% linear interpolation in 3d (for different pressures {5, 10, 15, 20})
function z = fun_p(i,p) % i \in 1:n; p \in 1:4
load Output_DATA/Test2_Output Y2;load PRO_TEST.mat P ;
sample=[3836 4305 1668 7676 7872 9311 987 1257 7393 7338 770 8350 7964 3538 2559];
sample_grid= [.99*Y2(sample,17) Y2(sample,18)];
x0_test=[.99*Y2(i,17) Y2(i,18)];
D1=zeros(1,length(sample));D2=D1;
    for j=1:length(sample)
        D1(j)=norm(x0_test(1)-sample_grid(j,1));
        D2(j)=norm(x0_test(2)-sample_grid(j,2));
    end
D= D1/(max(Y2(:,17))-min(Y2(:,17))) + D2/(max(Y2(:,18))-min(Y2(:,18)));
[a,b]=mink(D,3);c=1/sum(a); 
z=sum(c*a.*P(p,b));
end