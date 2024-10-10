%% press [ctrl]+[enter] to run code sections
addpath('Input_DATA','Scaled_model','Unscaled_model','Output_DATA')

%% get brine concentration and flow rate
clc
%init=[0; 70; 69.2]; % 8% brine
%init=[0; 70; 68.6]; % 7% brine
%init=[0; 70; 68.2]; % 6% brine
init=[0; 70; 67.2]; % 5% brine

output=fun_unscaled(init,0,'fig',1e4,1e-3);
ev(output)

%%
clc
% for ideal
%init=[0; 4.99; 5; 1.002; 0.0823]; 
%

% for ICP+ECP
%init=[0; 4.4; 5; 1.001; 0.0823]; 


%init=[-0.0176; 0; 10; 1.00001; 0.0823]; 
%init=[-0.0301; 0; 5; 1.0001; 0.0708]; 
%init=[-0.0385; 0; 5; 1.0001; 0.0655]; 
%init=[-0.0595; 0; 5; 1.0001; 0.0567]; 
%
output=fun_scaled(init,.3,'fig',1e4,1e-3);
ev(output,[5 8 11])

%init(1:3)=[-init(1); init(3); init(2)];
%output=fun_scaled(init,.3,'fig',1e4,1e-3);
%ev(output,[5 8 11])

%% PRO simulation
clc
rng default
n=100000;
for i=1:n
    Pd_L=5+15*rand();
    Pd_0=Pd_L-3*rand();
    Pf_0=1+.7*rand();
    init(:,i)=[0; Pd_0; Pd_L; Pf_0; 0.0823];
end
%%
results=zeros(18,n);
parfor i=1:n
    results(:,i)=fun_scaled(init,.3,'fig',1e4,1e-3);
end
save PRO_TEST results init
system('git add .');
system(['git commit -m "finished PRO test"']);
system('git push https://github.com/oliver-mx/GitMATLAB.git');


%% Test evaluation 
load PRO_TEST
close all; clc;
ev(results(:,1)); userNumber = input('Please enter a number of the set {1,2,3,...,18}: ');
f=figure(1);
f.Position= [117.6667 334.3333 1.4200e+03 999.3334];
% same colorbar
all_vectors = [Y1(:,userNumber); Y2(:,userNumber); Y3(:,userNumber); Y4(:,userNumber); Y5(:,userNumber); Y6(:,userNumber)];
cmap = jet(256); % Use the 'jet' colormap with 256 colors
cmin = min(all_vectors(:));
cmax = max(all_vectors(:));
% size
m_size=5;
%
scatter3(init(3,:), 10*(init(3,:)-init(2,:)),results(:,userNumber),m_size,results(:,userNumber),'filled');hold on; title('PRO test'); xlabel('Inlet pressure [bar]'); ylabel('Pressure drop [bar]'); axis equal;view(2); grid on;colormap(cmap);colorbar
%xlim([31 70]);ylim([10*delta_p_min 33]);yticks(10:10:30);yticklabels({'\Delta P = 1','\Delta P = 2','\Delta P = 3'})

































%% unscaeld test:
clc,
%
%[a1,b1]=fun_scaled([-0.0544004*.99; 0; 20; 1.1],.3,'fig',1e4,1e-3);
%ev(a1,[5 8 11])
%
n=1000;
TTT=linspace(1.000001,1.00005,n);
Testt=zeros(n,18);
parfor i = 1:n
  Testt(i,:)=fun_scaled([-0.0544004*.99; 0; 15; TTT(i)],.3,'sol',1e4,1e-3);
end

%% best Recovery
clc,
[a,b]=max(Testt(:,5));
disp('best opeating conditions:')
format long
disp([-0.0544004*.99; 0; 15; TTT(b)]')
format short
disp('Function result:')
[c,d]=fun_scaled([-0.0544004*.99; 0; 15; TTT(b)],.3,'fig',1e4,1e-3);
ev(c,[5 8 11])

%% 

%  -0.053856396000000                   0  20.000000000000000   1.000028908908909

% PD_net     = -1.7127 [kWh/m^2]
% REC_PRO =  36.3336 [%]
% Wastewater_in = 1.9255e-05 [m^3/s]
% C_dilluted    = 0.04644 [%]




































%%
system('git pull https://github.com/oliver-mx/GitMATLAB.git');

%%
system('git status');
system('git add .');
system('git commit -m "Further work on ERD..."');
system('git push https://github.com/oliver-mx/GitMATLAB.git');
