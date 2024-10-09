%% press [ctrl]+[enter] to run code sections
addpath('Input_DATA','Scaled_model','Unscaled_model','Output_DATA')

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
