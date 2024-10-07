%% press [ctrl]+[enter] to run code sections
addpath('Input_DATA','Scaled_model','Unscaled_model','Output_DATA')

%% SWRO with ERD:
clc;[a1,b1]=fun_unscaled([55.81;54.72],.2,'fig',1e4,1e-3);
ev(a1,[1 2 3 4 6 7 9 10 12])
%
% brine flow rate:    0.0266213 [kg/sm]
% brine concentation: 0.0558582 [1]
%

%% test only (counter-current) PRO:
%clc;[a2,b2]=fun_scaled([-0.0266213; 0 ; 5.4 ; 1.00002; 14.99],.3,'fig',1e4,1e-3);
% --> problem with wrong pressure in draw
clc;[a2,b2]=fun_scaled([0; 0 ; 5.4 ; 1.00003; 5.39942],.3,'fig',1e4,1e-3);
ev(a2,5)
%
% draw_in_flow       = -0.0266213 [kg/sm]
% draw concentration =  0.0509729 [1]
%

%% hybrid system mit same data:
clc,
output=fun_scaled([55.81,54.72,5.485,1.00002],.4,'fig',1e4,1e-3);
ev(output)





















%% test hybrid system
clc,
output=fun_scaled([55.81,54.72,5.485,1.000029760646511],.4,'fig',1e4,1e-3);
ev(output)

%% Simulation of the hybrid system I
clc;rng default;n=500;a2=zeros(n,18);
PD_in2 = 5.*ones(1,n) + 15.*rand(1,n);
PF_in2 = 1.00001.*ones(1,n) + 0.001.*rand(1,n);
parfor i=1:n
a2(i,:)=fun_scaled([55.81,54.72,PD_in2(i),PF_in2(i)],.4,'sol',1e4,1e-3);
end
%%
beep; clc,[t4, t5]=max(a2(:,5));
ev(a2(t5))


























