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

%% PRO test:
clc
%
init=[-0.0595; 0; 5; 1.0003; 0.0567]; 
%
output=fun_scaled(init,.3,'fig',1e4,1e-6);
ev(output,[5 6 7 8 10 11])

%init(1:3)=[-init(1); init(3); init(2)];
%output=fun_scaled(init,.3,'fig',1e4,1e-3);
%ev(output,[5 8 11])


%% -------- hybrid I test ----------
clc
%
init=[70; 67.2; 5; 1.0003]; 
%
output=fun_scaled(init,.4,'fig',1e4,1e-3);
ev(output)





























%% fmincon PRO pressure
A= []; b= []; Aeq=[]; beq=[]; lb = [-0.0595; 0; 5; 1; 0.0567]; ub = [-0.0595; 0; 20; 5; 0.0567];
foptions = optimoptions('fmincon','Display','iter','Algorithm','interior-point','StepTolerance', 1e-16, 'OptimalityTolerance',1e-6, 'MaxFunEvals',2000);
rng default
x0=[-0.0595; 0; 5.001; 1.0001; 0.0567]; 
[x1, x2] = fmincon(@(x)fun_scaled(x,.3,'REC',1e4,1e-6),x0,A,b,Aeq,beq,lb,ub,[],foptions);
%[x1, x2] = fmincon(@(x)fun_scaled(x,.3,'PD_net',1e4,1e-6),x0,A,b,Aeq,beq,lb,ub,[],foptions);
format long
disp(x1(3:4)')
format short
disp(-x2)

%%
output=fun_scaled(x1,.3,'fig',1e4,1e-6);
%ev(output,[5 6 7 8 10 11])
ev(output,[5 6])


% L^{PRO} = 7 * 0.9626
%
% max PD_net
% REC_PRO =  14.2838 [%]
% PD_net  =  0.078084 [W^2/m^2]
%





















%% PRO simulation
clc
% [-0.0595; 0; 5; 1.0001; 0.0567]
rng default
n=100;init=zeros(5,n);
for i=1:n
    Pd_L=4.9+1*rand();
    Pf_0=1+.00001+ 0.001*rand();
    init(:,i)=[-0.0595; 0; Pd_L; Pf_0; 0.0823];
end
%%
results=zeros(18,n);
parfor i=1:n
    results(:,i)=fun_scaled(init(:,i),.3,'fig',1e4,1e-6);
end
save PRO_TEST results init p_drop

%% Test evaluation 
load PRO_TEST
close all; clc;
ev(results(:,1)); userNumber = input('Please enter a number of the set {1,2,3,...,18}: ');
f=figure(1);
f.Position= [117.6667 334.3333 1.4200e+03 999.3334];
% same colorbar
cmap = jet(256); % Use the 'jet' colormap with 256 colors
% size
m_size=5;
%
scatter3(init(3,:), 10*p_drop^(-1).*(init(3,:)-init(2,:)),results(userNumber,:),m_size,results(userNumber,:),'filled');hold on; title('PRO test'); xlabel('Inlet pressure [bar]'); ylabel('Pressure drop [bar]'); axis equal;view(2); grid on;colormap(cmap);colorbar
xlim([5 20]);

%%
clc
[a,b1]=max(results(5,:));
disp(['highest Recovery rate = ',num2str(a),' [%]'])
[a,b2]=max(results(13,:));
disp(['highest PD_net        = ',num2str(a/37.1612),' [kWh]'])
disp(' ')
disp(init(2,b1))
disp(init(3,b1)-init(2,b1))
disp(init(2,b2))
disp(init(3,b2)-init(2,b2))









































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
