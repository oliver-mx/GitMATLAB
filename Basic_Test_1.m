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
[a3,b3]=fun_scaled([-.99*0.021369; 0 ; 5 ; 1.00003; 0.064655],.3,'fig',1e4,1e-4);
ev(a3,[5 6 9])
%
% REC_PRO =  97.6886 [%]
% PD_net  =  0.16009 [W/m^2]
% W_net = 5.949 [kWh]
%

%% hybrid system I:
clc;
[a4,b4]=fun_scaled([55.6, 54.72, 5, 1.00003],.4,'fig',1e4,1e-4);
ev(a4,[1 2 4 5 7])
%
% SEC_net = -1.7718 [kWh/m^3]
% FW      =  0.66334 [m^3/h]
% REC_RO  =  52.6892 [%]
% REC_PRO =  92.9808 [%]
% mix_M1  =  0.46881 
%

%% hybrid system II:
clc;
[a5,b5]=fun_scaled([55.6, 54.72, 5, 1.00003],.5,'fig',1e4,1e-3);
ev(a5,[1 2 4 5 7 8])
% 
% SEC_net = -1.8172 [kWh/m^3]
% FW      =  0.64235 [m^3/h]
% REC_RO  =  51.582 [%]
% REC_PRO =  92.9951 [%]
% mix_M1  =  0.48007 
% mix_M2  =  0 
%







