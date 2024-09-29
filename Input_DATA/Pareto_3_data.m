function [H,Z,swro_Z,ro_water,ro_salt,Mw,Ms,Rw,T0,eta,sigma,p_r,rho_r,C_r,swro_L,swro_alpha,swro_R,swro_KK,swro_x_r,swro_b1,swro_b2,J_r,swro_gamma,swro_gamma2,swro_W_r,L,alpha,R,KK,x_r,b1,b2,Q_r,gamma,gamma2,W_r,cE,pE,rho_E,J_sf_0,J_wf_0,Pd_0,Pd_L,Pf_L,Q_sf_0,pd_0,pf_0,pd_L,pf_L,HP_eff,LP_eff,T_eff,V_m,ERD_eff,ERD_fric,A_ERD,eta_ERD,mix_density,pw,pe,swro_beta_fix,beta_fix,mixer_ERD,version,fig,swro_KF,swro_KD,KF,KD]= Pareto_3_data(input1, input2)
%%  Pareto_3_data(input)
%
%   Data for Berechung der dritten Pareto front
%
%   Hybrid setup 2
%
%   option_data = 3
%
%   Input:
%       input1        -   [L, Pd_0, pd_0, pd_L, pf_0]
%       input2        -   x

%% model versions
version=zeros(1,10);
% version(1)=0 if co-current, 1 otherwise
version(2)=1;  % 0 = SWRO beta fixed
version(3)=0;  % 0 = PRO beta fixed
version(4)=1;  % 0 = ideal SWRO
version(5)=1;  % 0 = ideal PRO
% ERD
version(6)=4;
% version(6) = 0 --> only SWRO (no ERD)
% version(6) = 1 --> only SWRO (with ERD)
% version(6) = 2 --> only PRO
% version(6) = 3 --> SWRO-PRO hybrid system (with one ERD)
% version(6) = 4 --> SWRO-PRO hybrid system (with two ERDs) 
version(7)=1; % 1 = ICP and ECP for SWRO
version(8)=1; % 1 = ICP and ECP for PRO

%% Membrane unit properties
H = 1e-3;           % height of the membrane [m]
ro_water = 1000;    % mass density of water  [kg/m^3] 
ro_salt = 2165;     % mass density of salt [kg/m^3] 
Mw = 18;            % molecular weight of water [kg/kmol]
Ms = 58.44;         % molecular weight of salt  [kg/kmol]
Rw = 462;           % gas constant of water  [J/(kg K)] 
T0 = 297;           % temperature [K] 
eta = 1.3e-3;       % seawater viscosity [kg/(m s)]
p_r = 1e5;          % pressure [Pa]=[kg/ms^2]
rho_r = 1e3;        % density [kg/m^3]
C_r = 1;            % salt concentration [%]

%% SWRO
swro_Z=1;               % width of SWRO membrane [m]
swro_L=4;               % length of SWRO membrane [m]
if version(6)== 2; swro_L=1; end
swro_alpha = 5.0815e-9; % SWRO water permeablity co-efficient [s/m]
swro_R = 0.96;          % SWRO salt rejection rate 
swro_KK = 1e-2;         % SWRO ICP mass transfer coefficient
swro_KD = 1/swro_KK;    % SWRO ECP draw side mass transfer coefficient
swro_KF = 1/swro_KK;    % SWRO ECP fresh side mass transfer coefficient
swro_x_r= swro_L^2;     % x_r=swro_L^2 since x = linspace(0,1,n) (if x = linspace(0,swro_L,n) then x_r=swro_L;)     
swro_b1 = H/swro_x_r;                            % H/swro_L ratio
swro_b2 = swro_Z/swro_x_r;                       % Z/swro_L ratio
J_r = sqrt(H^3/swro_x_r*p_r*rho_r);              % flux [kg/s^2]
swro_gamma = swro_x_r * p_r *  swro_alpha /J_r;  % SWRO scaling factor - mass balance
swro_gamma2 = J_r^2./(swro_x_r^2 * p_r * rho_r); % SWRO scaling factor - momentum balance
swro_W_r = J_r*p_r/rho_r;                        % net work [W/m^2]
sigma = 0.999 ;                                  % Rejection coefficient
swro_beta_fix = 4.43e-4/J_r*swro_x_r;            % value for fixed SWRO beta [kg/sm^2]

%% PRO
Z = 1;              % width of the PRO membrane [m]
L = input1(1);      % length of the PRO membrane [m]      
alpha = 5.47e-9;    % water permeablity co-efficient [s/m]
R = 0.94;           % salt rejection rate [1]
KK = 7.13e2;        % mass transfer coefficient [sm^2/kg]
KD = 7.13e-2;       % PRO ECP draw side mass transfer coefficient [sm^2/kg]
KF = 7.13e-2;       % PRO ECP fresh side mass transfer coefficient [sm^2/kg]
x_r = L^2;          % x_r=L^2 since x = linspace(0,1,n) (if x = linspace(0,L,n) then x_r=L;)              
b1 = H/x_r;         % H/L ratio
b2 = Z/x_r;                            % Z/L ratio
Q_r = sqrt(H^3/x_r*p_r*rho_r);         % flux [kg/s^2]
gamma = x_r * p_r * alpha /Q_r ;       % PRO scaling factor - mass balance
gamma2 = Q_r^2./(x_r^2 * p_r * rho_r); % PRO scaling factor - momentum balance
W_r = Q_r*p_r/rho_r;                   % net work [W/m^2]
beta_fix = 1.71e-4/Q_r*x_r;            % value for fixed PRO beta [kg/sm^2]

%% Sea Water
cE= 35/983/C_r;                            % salt concentration in seawater
pE= 1e5/p_r;                               % external pressure
rho_E=(cE + 1)./(ro_water*cE/ro_salt + 1); % density of incomming seawater

%% SWRO operating conditions
J_sf_0 = 0;              % salt flux in fresh side at 0 
J_wf_0 = 0;              % water flux in fresh side at L               
Pf_L = pE;               % pressure fresh side at L

Pd_0 =  input1(2);       % pressure draw side at 0
Pd_L =  0;               % not needed in the hybrid system !

%% PRO operation conditions 
Q_sf_0 = 0;              % salt flux in fresh side at 0
pf_L = pE;               % pressure of fresh side at L

pd_0 = input1(3);         % pressure draw side at 0
pd_L = input1(4);       % pressure of fresh side at 0
pf_0 = input1(5);        % pressure draw side at L 

%% ERD/Turbine/Pump parameters
T_eff  = .95;               % turbine efficiency
HP_eff = .9;            	% high pressure pump efficiency
LP_eff = .95;               % low pressure pump efficiency
V_m = 0.052;                % Volumetric mixing
ERD_eff = .96;              % ERD unit pressure efficiency
ERD_fric = 5e-04;           % ERD friction coefficient
A_ERD = H*swro_Z;           % cross sectional area of ERD inflows/outflows
eta_ERD = 0.01;             % leak of high pressure brine
mix_density = 997/rho_r;  	% density of mixture in ERD
pw = 2.0;                   % water price [$/m^3]
pe = 0.3;                   % electricity price [$/kWh]
mixer_ERD = 1;              % PRO Draw outlet mixer adjustment (only if 2nd ERDs) (mixer_ERD=1 --> all flow to ERD2 no turbine needed)    

%% display figures
fig=[1,1,1,1]; % f(i)=1 --> figure i will be displayed

%% model specific changes:
% co-current
if pd_0 > pd_L; version(1)=0; else; version(1)=1; end            
% no SWRO needed (--> use trivial data)
if version(6) == 2; version(4)=0; Pd_0=40; Pd_L=35; end

end