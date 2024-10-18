%% press [ctrl]+[enter] to run code sections
addpath('Input_DATA','Scaled_model','Unscaled_model','Output_DATA')

%% single SWRO:
clc;
[a1,b1]=fun_unscaled([0;55.6;54.72],.1,'fig',1e4,1e-4);
%[a1,b1]=fun_scaled([0;55.6;54.72],.1,'fig',1e4,1e-4);
ev(a1,[1 2 4])
%
% SEC_net = -3.0825 [kWh/m^3]
% FW      =  0.68044 [m^3/h]
% REC_RO  =  53.5756 [%]
%


%% SWRO with ERD:
clc;
[a2,b2]=fun_unscaled([0;55.6;54.72],.2,'fig',1e4,1e-4);
%[a2,b2]=fun_scaled([0;55.6;54.72],.2,'fig',1e4,1e-4);
ev(a2,[1 2 4 7])
%
% SEC_net = -1.7798 [kWh/m^3]
% FW      =  0.66103 [m^3/h]
% REC_RO  =  52.5707 [%]
% mix_M1  =  0.47003 
%

%% PRO 
clc;
[a3,b3]=fun_scaled([-0.0024; 0 ; 5 ; 1.00003; 0.05],.3,'fig',1e4,1e-3);
ev(a3,[5 6])
%
% REC_PRO =  94.4012 [%]
% PD_net  =  0.25498 [W/m^2]
%

%% PRO 
clc;
[a3,b3]=fun_scaled([-0.0217466; 0 ; 5.4 ; 1.00003; 0.0573256],.3,'fig',1e4,1e-6);
ev(a3,[5 6])
%
% REC_PRO =  96.853 [%]
% PD_net  =  0.1579 [W^2/m^2]
%

%% hybrid system I:
clc;
[a4,b4]=fun_scaled([55.6,54.72,10,1.001],.4,'fig',1e4,1e-6);
ev(a4,[1 2 4 5 18 19])
%
% SEC_net = -1.7746 [kWh/m^3]
% FW      =  0.58614 [m^3/h]
% REC_RO  =  48.8954 [%]
% REC_PRO =  9.0944 [%]
% mix_M1  = 0.50721 
% F_1-F_2 = 0.0082248 [kg m^2/s^3]
%

%% hybrid system II:
clc;
[a5,b5]=fun_scaled([55.6,54.72,10,1.001],.5,'fig',1e4,1e-6);
ev(a5,[1 2 4 5 18 19 20 21])
% 
% SEC_net = -1.8681 [kWh/m^3]
% FW      =  0.57091 [m^3/h]
% REC_RO  =  48.0191 [%]
% REC_PRO =  9.1264 [%]
% mix_M1  = 0.51605 
% F_1-F_2 = 0.52522 [kg m^2/s^3]
% mix_M2  = 0 
% f_1-f_2 = 3.0622 [kg m^2/s^3]
%



%%
clc;close all;
I=[1.001 1.0003, 1.0001];
for i=1:length(I)
a=fun_scaled([55.6,54.72,5,I(i)],.5,'sol',1e4,1e-6);
format long
disp(I(i))
format short
disp(a(1))
disp('______________')
end















