%% Basic simulations for all cases (ideal, ICP, ICP+ECP)
% press [ctrl]+[enter] to run code sections
addpath('Input_DATA','Scaled_model','Unscaled_model','Output_DATA')

%% test random data
rng default
X0=zeros(2,10000);
for i=1:10000
    X0(1,i)= 31+39*rand();
    X0(2,i)= X0(1,i) - .1 - 3.3*rand();
end
%% Simulate the systems
parfor i=1:10000
    Y1(i,:)=fun_unscaled(X0(:,i),-.1,'sol',1e4,1e-3);Y2(i,:)=fun_unscaled(X0(:,i),-.2,'sol',1e4,1e-3);Y3(i,:)=fun_unscaled(X0(:,i),-.3,'sol',1e4,1e-3);
    Y4(i,:)=fun_unscaled([0; X0(:,i)],-.4,'sol',1e4,1e-3);Y5(i,:)=fun_unscaled([0; X0(:,i)],-.5,'sol',1e4,1e-3);Y6(i,:)=fun_unscaled([0; X0(:,i)],-.6,'sol',1e4,1e-3);
    %Y7(i,:)=fun_scaled(X0_hybrid,-.7,'sol',1e4,1e-3);Y8(i,:)=fun_scaled(X0_hybrid,-.8,'sol',1e4,1e-3);Y9(i,:)=fun_scaled(X0_hybrid,-.9,'sol',1e4,1e-3);
end
save Output_DATA/Vergleich_mit_ideal X0 X0_hybrid Y1 Y2 Y3 Y4 Y5 Y6 Y7 Y8 Y9

system('git add .'); system('git commit -m "Vergleich_mit_ideal"');system('git push https://github.com/oliver-mx/GitMATLAB.git');













%% Test evaluation 
load("Output_DATA/Vergleich_mit_ideal.mat")
close all; clc;
ev(Y1(1,:)); userNumber = input('Please enter a number of the set {1,2,3,...,18}: ');
f=figure(1);
f.Position= [117.6667 334.3333 1.4200e+03 999.3334];tiledlayout(2,3);
% introduce minimum pressure drop
delta_p_min=1; % 0.1 default
for i=1:10000
    if X0(1,i)-X0(2,i)<delta_p_min
            Y1(i,:)=NaN(1,18);Y2(i,:)=NaN(1,18);Y3(i,:)=NaN(1,18);Y4(i,:)=NaN(1,18);Y5(i,:)=NaN(1,18);Y6(i,:)=NaN(1,18);
    end
end
% draw level lines for 4e-4, 5e-4, 6e-4, 7e-4, 8e-4.
A1=zeros(2,1);A2=A1;A3=A1;A4=A1;A5=A1;tol=5e-6;
B1=zeros(2,1);B2=B1;B3=B1;B4=B1;B5=B1;
if 1==0
    for j=1:10000
        if Y3(j,6) > 4e-4-tol && Y3(j,6) < 4e-4+tol; A1(:,j)=X0(:,j); end
        if Y3(j,6) > 5e-4-tol && Y3(j,6) < 5e-4+tol; A2(:,j)=X0(:,j); end
        if Y3(j,6) > 6e-4-tol && Y3(j,6) < 6e-4+tol; A3(:,j)=X0(:,j); end
        if Y3(j,6) > 7e-4-tol && Y3(j,6) < 7e-4+tol; A4(:,j)=X0(:,j); end
        if Y3(j,6) > 8e-4-tol && Y3(j,6) < 8e-4+tol; A5(:,j)=X0(:,j); end
        if Y6(j,6) > 4e-4-tol && Y6(j,6) < 4e-4+tol; B1(:,j)=X0(:,j); end
        if Y6(j,6) > 5e-4-tol && Y6(j,6) < 5e-4+tol; B2(:,j)=X0(:,j); end
        if Y6(j,6) > 6e-4-tol && Y6(j,6) < 6e-4+tol; B3(:,j)=X0(:,j); end
        if Y6(j,6) > 7e-4-tol && Y6(j,6) < 7e-4+tol; B4(:,j)=X0(:,j); end
        if Y6(j,6) > 8e-4-tol && Y6(j,6) < 8e-4+tol; B5(:,j)=X0(:,j); end
    end
end
% same colorbar
all_vectors = [Y1(:,userNumber); Y2(:,userNumber); Y3(:,userNumber); Y4(:,userNumber); Y5(:,userNumber); Y6(:,userNumber)];
cmap = jet(256); % Use the 'jet' colormap with 256 colors
cmin = min(all_vectors(:));
cmax = max(all_vectors(:));
%
nexttile
scatter3(X0(1,:), 10*(X0(1,:)-X0(2,:)),Y1(:,userNumber),1.5,Y1(:,userNumber),'filled');hold on; title('simple SWRO (ideal)'); xlabel('Seawater inlet pressure [bar]'); ylabel('Pressure drop [bar]'); axis equal;view(2); grid on;colormap(cmap);colorbar
xlim([30.5 70]);ylim([10*delta_p_min 33]);yticks(10:10:30);yticklabels({'\Delta P = 1','\Delta P = 2','\Delta P = 3'})
nexttile
scatter3(X0(1,:), 10*(X0(1,:)-X0(2,:)),Y2(:,userNumber),1.5,Y2(:,userNumber),'filled');hold on; title('simple SWRO (with ICP)'); xlabel('Seawater inlet pressure [bar]'); ylabel('Pressure drop [bar]'); axis equal;view(2); grid on;colormap(cmap);colorbar
xlim([30.5 70]);ylim([10*delta_p_min 33]);yticks(10:10:30);yticklabels({'\Delta P = 1','\Delta P = 2','\Delta P = 3'})
nexttile
scatter3(X0(1,:), 10*(X0(1,:)-X0(2,:)),Y3(:,userNumber),1.5,Y3(:,userNumber),'filled');hold on; title('simple SWRO (ICP+ECP)'); xlabel('Seawater inlet pressure [bar]'); ylabel('Pressure drop [bar]'); axis equal;view(2); grid on;colormap(cmap);colorbar
xlim([30.5 70]);ylim([10*delta_p_min 33]);yticks(10:10:30);yticklabels({'\Delta P = 1','\Delta P = 2','\Delta P = 3'})
scatter([A1(1,:) A2(1,:) A3(1,:) A4(1,:) A5(1,:)], 10.*[(A1(1,:)-A1(2,:)) (A2(1,:)-A2(2,:)) (A3(1,:)-A3(2,:)) (A4(1,:)-A4(2,:)) (A5(1,:)-A5(2,:)) ],5,'k');hold on
nexttile
scatter3(X0(1,:), 10*(X0(1,:)-X0(2,:)),Y4(:,userNumber),1,Y4(:,userNumber),'filled');hold on; title('SWRO+ERD (ideal)'); xlabel('Seawater inlet pressure [bar]'); ylabel('Pressure drop [bar]'); axis equal;view(2); grid on;colormap(cmap);colorbar
xlim([30.5 70]);ylim([10*delta_p_min 33]);yticks(10:10:30);yticklabels({'\Delta P = 1','\Delta P = 2','\Delta P = 3'})
nexttile
scatter3(X0(1,:), 10*(X0(1,:)-X0(2,:)),Y5(:,userNumber),1,Y5(:,userNumber),'filled');hold on; title('SWRO+ERD (with ICP)'); xlabel('Seawater inlet pressure [bar]'); ylabel('Pressure drop [bar]'); axis equal;view(2); grid on;colormap(cmap);colorbar
xlim([30.5 70]);ylim([10*delta_p_min 33]);yticks(10:10:30);yticklabels({'\Delta P = 1','\Delta P = 2','\Delta P = 3'})
nexttile
scatter3(X0(1,:), 10*(X0(1,:)-X0(2,:)),Y6(:,userNumber),1,Y6(:,userNumber),'filled');hold on; title('SWRO+ERD (ICP+ECP)'); xlabel('Seawater inlet pressure [bar]'); ylabel('Pressure drop [bar]'); axis equal;view(2); grid on;colormap(cmap);colorbar
xlim([30.5 70]);ylim([10*delta_p_min 33]);yticks(10:10:30);yticklabels({'\Delta P = 1','\Delta P = 2','\Delta P = 3'})
scatter([B1(1,:) B2(1,:) B3(1,:) B4(1,:) B5(1,:)], 10.*[(B1(1,:)-B1(2,:)) (B2(1,:)-B2(2,:)) (B3(1,:)-B3(2,:)) (B4(1,:)-B4(2,:)) (B5(1,:)-B5(2,:)) ],5,'k');hold on
%nexttile
%scatter3(X0(1,:), 10*(X0(1,:)-X0(2,:)),Y7(:,userNumber),1,Y7(:,userNumber),'filled');hold on; title('Hybrid (ideal)'); xlabel('Seawater inlet pressure [bar]'); ylabel('Pressure drop [bar]'); axis equal;view(2); grid on;colormap(cmap);colorbar
%xlim([30.5 70]);ylim([10*delta_p_min 33]);yticks(10:10:30);yticklabels({'\Delta P = 1','\Delta P = 2','\Delta P = 3'})
%nexttile
%scatter3(X0(1,:), 10*(X0(1,:)-X0(2,:)),Y8(:,userNumber),1,Y8(:,userNumber),'filled');hold on; title('Hybrid (with ICP)'); xlabel('Seawater inlet pressure [bar]'); ylabel('Pressure drop [bar]'); axis equal;view(2); grid on;colormap(cmap);colorbar
%xlim([30.5 70]);ylim([10*delta_p_min 33]);yticks(10:10:30);yticklabels({'\Delta P = 1','\Delta P = 2','\Delta P = 3'})
%nexttile
%scatter3(X0(1,:), 10*(X0(1,:)-X0(2,:)),Y9(:,userNumber),1,Y9(:,userNumber),'filled');hold on; title('Hybrid (ICP+ECP)'); xlabel('Seawater inlet pressure [bar]'); ylabel('Pressure drop [bar]'); axis equal;view(2); grid on;colormap(cmap);colorbar
%xlim([30.5 70]);ylim([10*delta_p_min 33]);yticks(10:10:30);yticklabels({'\Delta P = 1','\Delta P = 2','\Delta P = 3'})

%% evalue pareto plot
load("Output_DATA/Simulation_output.mat")
close all; clc
for i=1:10000
    if X0(1,i)-X0(2,i)>4 || X0(1,i)-X0(2,i)<.1
        Y1(i,1)=-1;Y1(i,2)=-1;Y2(i,1)=-1;Y2(i,2)=-1;Y3(i,1)=-1;Y3(i,2)=-1;Y4(i,1)=-1;Y4(i,2)=-1;Y5(i,1)=-1;Y5(i,2)=-1;Y6(i,1)=-1;Y6(i,2)=-1;
    end
end
f=figure(1); f.Position= [578.3333 778.3333 1.5567e+03 487.3333];tiledlayout(1,3);
nexttile
scatter(Y1(:,1),Y1(:,2),'r');hold on; title('simple SWRO (ideal)'); xlabel('SEC_{net} [kWh/m^3]'); ylabel('FW [m^3/h]');
scatter(Y4(:,1),Y4(:,2),'b');hold on;grid on;
xlim([-10 0]);ylim([0 2.2]);
nexttile
scatter(Y2(:,1),Y2(:,2),'r');hold on; title('simple SWRO (ICP)'); xlabel('SEC_{net} [kWh/m^3]'); ylabel('FW [m^3/h]');
scatter(Y5(:,1),Y5(:,2),'b');hold on;grid on;
xlim([-10 0]);ylim([0 2.2]);
nexttile
scatter(Y3(:,1),Y3(:,2),'r');hold on; title('simple SWRO (ICP+ECP)'); xlabel('SEC_{net} [kWh/m^3]'); ylabel('FW [m^3/h]');
scatter(Y6(:,1),Y6(:,2),'b');hold on;grid on;
xlim([-10 0]);ylim([0 2.2]);