%% press [ctrl]+[enter] to run code sections
addpath('Input_DATA','Scaled_model','Unscaled_model','Output_DATA')

%% single SWRO:
clc;
[a1,b1]=fun_unscaled([0;55.6;54.72],.1,'fig',1e4,1e-4);
ev(a1,[1 2 4 5 12])
%
% SEC_net = -3.0825 [kWh/m^3]
% FW      =  0.68044 [m^3/h]
% REC_RO  =  53.5756 [%]
% REC_PRO =  NaN [%]
% mix_M1  =  NaN 
%

%% SWRO with ERD:
clc;
[a2,b2]=fun_unscaled([0;55.6;54.72],.2,'fig',1e4,1e-4);
ev(a2,[1 2 4 5 12])
%
% SEC_net = -1.7798 [kWh/m^3]
% FW      =  0.66086 [m^3/h]
% REC_RO  =  52.5618 [%]
% REC_PRO =  NaN [%]
% mix_M1  =  0.46533 
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
[a4,b4]=fun_scaled([55.6,54.72,5.4,1.00003],.4,'fig',1e4,1e-6);
ev(a4,[1 2 4 5 12])
% no mixing:                          % with mixing:
% SEC_net = -1.7832 [kWh/m^3]         % SEC_net = -1.7867 [kWh/m^3]
% FW      =  0.6003 [m^3/h]           % FW      =  0.58599 [m^3/h]
% REC_RO  =  49.6975 [%]              % REC_RO  =  48.8871 [%]
% REC_PRO =  92.1988 [%]              % REC_PRO =  92.0681 [%]
% mix_M1  =  0.49412                  % mix_M1  =  0.50214 
%

%% hybrid system II:
clc;
[a5,b5]=fun_scaled([55.6,54.72,5.4,1.00003],.5,'fig',1e4,1e-6);
ev(a5,[1 2 4 5 12])
% no mixing:                          % with mixing:
% SEC_net = -1.7687 [kWh/m^3]         % SEC_net = -1.7763 [kWh/m^3]
% FW      =  0.6003 [m^3/h]           % FW      =  0.56987 [m^3/h]
% REC_RO  =  49.6975 [%]              % REC_RO  =  47.9588 [%]
% REC_PRO =  92.46 [%]                % REC_PRO =  92.3774 [%]
% mix_M1        = 0.49412             % mix_M1        = 0.51141 
%

