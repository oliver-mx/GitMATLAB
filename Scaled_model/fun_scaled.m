function [output1, output2] = fun_scaled(input,option_data,obj,option_mesh,option_BVP)
%%  fun   Solves BVP for SWRO-PRO hybrid system using bvp5c
%
%       fun(input,option_data,obj,option_mesh,option_BVP) will solve the
%       nonlinear BVP of the selected SWRO-PRO system. In "option_data" 
%       all Informations of the system are given with exception of the 
%       additional "input". 
%
%   Input:
%       input         -   Input for the Data function
%       option_data   -   Selects Data function used: data(input)
%       obj           -   'SEC', 'FW', 'Rev', 'Pareto','sol', 'fig';
%       option_mesh   -   NMax of bvp5c
%       option_BVP    -   RelTol of bvp5c
%
%   Output (using 'sol', or 'fig'):
%       output1 = [SEC_net, FW, Rev, SWRO_Rec, PRO_Rec, ...
%                  RO_inflow, Permeate_outflow, Wastewater_inflow, C_permeate, C_brine, C_dilluted, mix_M1...
%                  W_net, -W_p1, -W_p2,  -W_p3, -W_p4, W_t];
%       output2 = sol.stats;
%
%   Example:
%       x0=[55.81e5;54.72e5];       
%       data=.11;      
%       obj='sol';  
%       mesh=1e4;
%       tol=1e-4;  
%       [a,b]=fun_unscaled(x0,data,obj,mesh,tol);
%

if nargout>2
    tic;
end

%% Read data
if option_data == .1; DATA = @(x)Test_01_data(input); end
if option_data == .2; DATA = @(x)Test_02_data(input); end
if option_data == .3; DATA = @(x)Test_03_data(input); end
if option_data == .4; DATA = @(x)Test_04_data(input); end
if option_data == .5; DATA = @(x)Test_05_data(input); end
%
if option_data == -.1; DATA = @(x)data1(input); end
if option_data == -.2; DATA = @(x)data2(input); end
if option_data == -.3; DATA = @(x)data3(input); end
if option_data == -.4; DATA = @(x)data4(input); end
if option_data == -.5; DATA = @(x)data5(input); end
if option_data == -.6; DATA = @(x)data6(input); end
if option_data == -.7; DATA = @(x)data7(input); end
if option_data == -.8; DATA = @(x)data8(input); end
if option_data == -.9; DATA = @(x)data9(input); end
%
if option_data == -1; DATA = @(x)Lee_data(input);end
if option_data == 0; DATA = @(x)Senthil_data(input);end
if option_data == 1; DATA = @(x)Case_1_data(input); end
if option_data == 2; DATA = @(x)Case_2_data(input); end
if option_data == 3; DATA = @(x)Case_3_data(input); end
%

[H, Z, swro_Z, ro_water, ro_salt, Mw, Ms, Rw, T0, eta, sigma, p_r, rho_r, C_r, swro_L, swro_alpha, swro_KK, swro_x_r, swro_b1, swro_b2, J_r, swro_gamma, swro_gamma2, swro_W_r, L, alpha, KK, x_r, b1, b2, Q_r, gamma, gamma2, W_r, cE, pE, rho_E, J_sf_0, J_wf_0, Pd_0, Pd_L, Pf_L, Q_sf_0, pd_0, pf_0, pd_L, pf_L, HP_eff, LP_eff, T_eff, V_m, ERD_eff, ERD_fric, A_ERD, eta_ERD, mix_density, pw, pe, swro_beta_fix, beta_fix, mix_M1, mix_M3, version, fig, swro_KF, swro_KD, KF, KD]...
    = DATA(1);

%% Set options
ode_options = bvpset('Stats','off','NMax', option_mesh, 'RelTol', option_BVP);      

n=25;
x = linspace(0,1,n);    
  
%% Create initial guess for the BVP
J_wd_0 = (983/1018)/J_r;          % guess for J_wd(0)
J_wf_0 = 0.01353/J_r;             % should be close to 0
Pf_0   = 1.1*Pf_L;                % guess for P_f(0)
Q_wd_0 = 0.01353*(983/1018)/Q_r;  % guess for Q_wd(0) 
Q_wf_0 = 0.01353/Q_r;             % guess for Q_wf(0)

if version(1)==0; y_init = [cE; J_wd_0; J_sf_0; J_wf_0; Pd_0; Pf_0; cE;  Q_wd_0; Q_sf_0; Q_wf_0; pd_0; pf_0]; end
if version(1)==1; y_init = [cE; J_wd_0; J_sf_0; J_wf_0; Pd_0; Pf_0; cE; -Q_wd_0; Q_sf_0; Q_wf_0; pd_0; pf_0]; end
solinit = bvpinit(x,y_init);

try
%% Solve the BVP
if version(6) == 0; sol = bvp5c(@(x,J_p)ODEsystem(x, J_p, DATA), @(ya,yb)Boundary1(ya,yb, DATA),solinit,ode_options); end
if version(6) == 1; sol = bvp5c(@(x,J_p)ODEsystem(x, J_p, DATA), @(ya,yb)Boundary2(ya,yb, DATA),solinit,ode_options); end
if version(6) == 2; sol = bvp5c(@(x,J_p)ODEsystem(x, J_p, DATA), @(ya,yb)Boundary1(ya,yb, DATA),solinit,ode_options); end
if version(6) == 3; sol = bvp5c(@(x,J_p)ODEsystem(x, J_p, DATA), @(ya,yb)Boundary3(ya,yb, DATA),solinit,ode_options); end       
if version(6) == 4; sol = bvp5c(@(x,J_p)ODEsystem(x, J_p, DATA), @(ya,yb)Boundary4(ya,yb, DATA),solinit,ode_options); end   

if version(6)>0 && version(6)~=2
    y = deval(sol,x); Y = y'; Y = real(Y);
    mix_M1 = ((Y(end,1).*Y(end,2)+Y(end,2)))/(Y(1,1).*Y(1,2)+Y(1,2));  
    % data with new mix_M1
    if option_data == .1; DATA = @(x)Test_01_data(input,mix_M1); end
    if option_data == .2; DATA = @(x)Test_02_data(input,mix_M1); end
    if option_data == .3; DATA = @(x)Test_03_data(input,mix_M1); end
    if option_data == .4; DATA = @(x)Test_04_data(input,mix_M1); end
    if option_data == .5; DATA = @(x)Test_05_data(input,mix_M1); end
    %
    if option_data == -.1; DATA = @(x)data1(input,mix_M1); end
    if option_data == -.2; DATA = @(x)data2(input,mix_M1); end
    if option_data == -.3; DATA = @(x)data3(input,mix_M1); end
    if option_data == -.4; DATA = @(x)data4(input,mix_M1); end
    if option_data == -.5; DATA = @(x)data5(input,mix_M1); end
    if option_data == -.6; DATA = @(x)data6(input,mix_M1); end
    if option_data == -.7; DATA = @(x)data7(input,mix_M1); end
    if option_data == -.8; DATA = @(x)data8(input,mix_M1); end
    if option_data == -.9; DATA = @(x)data9(input,mix_M1); end
    %
    %if option_data == -1; DATA = @(x)Lee_data(input,mix_M1);end
    if option_data == 0; DATA = @(x)Senthil_data(input,mix_M1);end
    if option_data == 1; DATA = @(x)Case_1_data(input,mix_M1); end
    if option_data == 2; DATA = @(x)Case_2_data(input,mix_M1); end
    if option_data == 3; DATA = @(x)Case_3_data(input,mix_M1); end
    % solve with bvp
    if version(6) == 0; sol = bvp5c(@(x,J_p)ODEsystem(x, J_p, DATA), @(ya,yb)Boundary1(ya,yb, DATA),solinit,ode_options); end
    if version(6) == 1; sol = bvp5c(@(x,J_p)ODEsystem(x, J_p, DATA), @(ya,yb)Boundary2(ya,yb, DATA),solinit,ode_options); end
    if version(6) == 2; sol = bvp5c(@(x,J_p)ODEsystem(x, J_p, DATA), @(ya,yb)Boundary1(ya,yb, DATA),solinit,ode_options); end
    if version(6) == 3; sol = bvp5c(@(x,J_p)ODEsystem(x, J_p, DATA), @(ya,yb)Boundary3(ya,yb, DATA),solinit,ode_options); end       
    if version(6) == 4; sol = bvp5c(@(x,J_p)ODEsystem(x, J_p, DATA), @(ya,yb)Boundary4(ya,yb, DATA),solinit,ode_options); end
end

% Evaluate the solution
if any([0,0,0,0] ~= fig)
    n=100; % higher resolution plots
end
x=linspace(0,1,n);
y = deval(sol,x); Y = y'; Y = real(Y);

%% SWRO evaluation
C_d = Y(:,1);
C_f = Y(:,3)./Y(:,4);
J_sd = Y(:,1).*Y(:,2);
J_wd = Y(:,2); 
J_d = Y(:,1).*Y(:,2)+Y(:,2);
J_sf = Y(:,3);
J_wf = Y(:,4);
J_f = Y(:,3)+Y(:,4);
P_d = Y(:,5); 
P_f = Y(:,6); 
% same quantities as in "ODEsystem.m" 
swro_p_osm_d = ro_water*Rw*T0*log(1 + 2*Mw*Y(:,1)/Ms)/p_r;
swro_p_osm_f = ro_water*Rw*T0*log(1 + 2*Mw*Y(:,3)/Ms./Y(:,4))/p_r;
% local density
swro_local_ro_d = (Y(:,1) + 1)./(Y(:,1)./ro_salt + ones(n,1)./ro_water)/rho_r;
swro_local_ro_f = (Y(:,3) + Y(:,4))./(Y(:,3)./ro_salt + Y(:,4)./ro_water)/rho_r;
    % Salt permeability
    if version(4) == 0
        swro_beta = zeros(n,1);
    else
        swro_beta = swro_beta_fix/p_r/swro_alpha.*ones(n,1);
    end
    % Water Permeate flux J_win(x)
    if version(4) == 0 % ideal
        J_cross = ((Y(:,5) - Y(:,6)) - sigma*(swro_p_osm_d - swro_p_osm_f));
    elseif version(7) == 0 % ICP
        J_cross = ((Y(:,5) - Y(:,6)) - sigma.*(swro_p_osm_d - swro_p_osm_f))./(1 + p_r*swro_alpha*swro_KK*sigma.*(swro_p_osm_d - swro_p_osm_f));
    else % ICP+ECP
        J_cross = ((Y(:,5) - Y(:,6)) .* (ones(n,1) + swro_beta_fix.*ones(n,1).* (-1 ./ swro_KD + 1.*swro_KK + 1 ./ swro_KF)) - sigma .* (swro_p_osm_d - swro_p_osm_f)) ./ (ones(n,1) + swro_beta_fix.*ones(n,1) .* (-1 / swro_KD + swro_KK + 1 / swro_KF) - p_r .* swro_alpha .*sigma .* (-swro_p_osm_d ./ swro_KD + swro_p_osm_f .* swro_KK + swro_p_osm_f ./ swro_KF));
    end
    % Salt Permeate J_sin(x)
    if version(7) == 1 % ICP+ECP
        J_sin = swro_beta .* ((C_d - C_f) + C_d .* J_cross.*p_r.*swro_alpha ./ swro_KD - (C_f) .* J_cross.*p_r.*swro_alpha .* (swro_KK + ones(n,1) ./ swro_KF)) ./ (ones(n,1) + swro_beta_fix.*ones(n,1) .* (-1 ./ swro_KD + swro_KK + 1 ./ swro_KF));
    else
        J_sin=swro_beta_fix.*(C_d-C_f);
    end

%% PRO evaluation
c_d = Y(:,1+6);                    % concentration of draw side
c_f = Y(:,3+6)./Y(:,4+6);          % concentration of fresh side
Q_sd = Y(:,1+6).*Y(:,2+6);          % salt flow rate in draw side
Q_wd = Y(:,2+6);                   % water flow rate in draw side
Q_d = Y(:,1+6).*Y(:,2+6)+Y(:,2+6); % total flow rate draw side
Q_sf = Y(:,3+6);                   % salt flow rate in fresh side
Q_wf = Y(:,4+6);                   % water flow rate in fresh side
Q_f = Y(:,3)+Y(:,4+6);            % total flow rate fresh side
p_d = Y(:,5+6);                     % Pressure draw side
p_f = Y(:,6+6);                     % Pressure fresh side
% same quantities as in "ODEsystem.m" 
p_osm_d = ro_water*Rw*T0*log(1 + 2*Mw*Y(:,1+6)/Ms)/p_r;
p_osm_f = ro_water*Rw*T0*log(1 + 2*Mw*Y(:,3+6)/Ms./Y(:,4+6))/p_r;
% local density
local_ro_d = (Y(:,1+6) + 1)./(Y(:,1+6)./ro_salt + ones(n,1)./ro_water)/rho_r;
local_ro_f = (Y(:,3+6) + Y(:,4+6))./(Y(:,3+6)./ro_salt + Y(:,4+6)./ro_water)/rho_r;
    % Salt permeability
    if version(5) == 0
        beta = zeros(n,1);
    else
        beta = beta_fix/p_r/alpha.*ones(n,1);
    end
    % Water Permeate flux Q_win(x)
    if version(5) == 0 % ideal
        Q_cross = ((p_osm_d - p_osm_f) - (Y(:,5) - Y(:,6)));
    elseif version(8) == 0 % ICP
        Q_cross = ( -(Y(:,5+6) - Y(:,6+6)) +(p_osm_d - p_osm_f))./(1 + p_r*alpha*KK*sigma.*(p_osm_d - p_osm_f));
    else % ICP+ECP
        Q_cross = ((p_osm_d - p_osm_f) - (Y(:,5+6) - Y(:,6+6)) .* (ones(n,1) + beta_fix.*ones(n,1).* (1 ./ KD + 1.*KK + 1 ./ KF))) ./ (ones(n,1) + beta_fix.*ones(n,1) .* (1 / KD + KK + 1 / KF) + p_r .* alpha .* (p_osm_d ./ KD + p_osm_f .* KK + p_osm_f ./ KF));
    end
    % Salt Permeate Q_sin(x)
    if version(8) == 1 % ICP+ECP
        Q_sin = beta .* ((c_d - c_f) - c_d .* Q_cross.*p_r.*alpha ./ KD - (c_f) .* Q_cross.*p_r.*alpha .* (KK + ones(n,1) ./ KF)) ./ (ones(n,1) + beta_fix.*ones(n,1) .* (1 ./ KD + KK + 1 ./ KF));
    else
        Q_sin = beta_fix.*(c_d-c_f);
    end

%% Output of System
% seawater enters the system
J_E = J_d(1);
J_wE = J_E/(cE+1);
J_sE = J_E-J_wE;

%% Version(6)=0 -->  only SWRO (no ERD)
if version(6)==0
    W_p1 = 1/HP_eff * (pE - P_d(1))*(J_d(1) *swro_Z)/rho_E; W_p1 = W_p1*swro_W_r;
    W_p2 = 0; W_p3=0; W_p4=0; W_t=0;
end

%% Version(6)=1 -->  only SWRO (with ERD)
if version(6)==1
    % flow from M1 to ERD1
    J_M1 = J_E * (1*mix_M1);
    J_w_M1 = J_wE* (1*mix_M1);
    J_s_M1 = J_sE* (1*mix_M1);
    % flow from ERD1 to RO
    rho_ERD= V_m*(swro_local_ro_d(end) - rho_E/rho_r)+rho_E/rho_r;
    C_ERD= -ro_salt/rho_r*(rho_ERD-ro_water/rho_r)/(ro_water/rho_r*(rho_ERD-ro_salt/rho_r));
    J_ERD = J_M1;
    J_wERD= J_ERD/(C_ERD+1);
    J_sERD= J_ERD - J_wERD;
    % flow from ERD1 to EXIT
    J_exit = J_d(end)*(1-eta_ERD);
    c_exit= (J_s_M1 + J_sd(end)*(1-eta_ERD) - J_sERD)/(J_w_M1 + J_wd(end)*(1-eta_ERD) - J_wERD);
    J_w_exit= J_exit/(c_exit+1);
    rho_exit = (J_exit)/(((c_exit*J_w_exit)/ro_salt*rho_r)+(J_w_exit/ro_water*rho_r));
    % energy balance
    f_1 = ERD_fric * mix_density * (J_exit*swro_b2/rho_exit)*((J_M1*swro_b2)/(rho_E/rho_r*swro_b1*swro_b2))^2; 
    f_2 = ERD_fric * mix_density * (J_ERD*swro_b2/rho_ERD)*((J_d(end)*swro_b2)/(swro_local_ro_d(end)*swro_b1*swro_b2))^2;
    f_r= (J_r*p_r*swro_x_r/rho_r) / (J_r^3/rho_r^2/swro_x_r); % scaling factor between f and P*J*Z/rho
    pERD = rho_ERD*(ERD_eff*(P_d(end)*J_d(end)*swro_b2*(1-eta_ERD)/swro_local_ro_d(end)-pE*J_exit*swro_b2/rho_exit - f_2/f_r) + pE*J_M1*swro_b2/rho_E*rho_r + f_1/f_r)/J_ERD/swro_b2;
    norm_f1=min(abs(f_1-f_2),abs(f_2-f_1))/J_r/p_r/swro_x_r*rho_r; % dimensions: J_r*p_r*swro_x_r*/rho_r [W];
    % pumps
    W_p1 = 1/HP_eff * (pE - P_d(1))*((J_E-J_M1) *swro_Z)/rho_E; W_p1 = W_p1*swro_W_r;
    W_p3 = 1/HP_eff * (pERD - P_d(1))*(J_ERD*swro_Z)/rho_ERD/rho_r; W_p3 = W_p3*swro_W_r;
    W_p2 = 0; W_p4=0; W_t=0;
end

%% Version(6)=2 --> only PRO
if version(6)==2
    if version(1) ==0
        W_p2= 1/LP_eff * (pE - p_d(1))*(Q_d(1)*Z)./local_ro_d(1)/rho_r; W_p2 = W_p2*W_r;
        W_t = T_eff * (p_d(end)-pE)*(Q_d(end)*Z)/local_ro_d(end)/rho_r; W_t = W_t*W_r;
    else
        W_p2= 1/LP_eff * (pE - p_d(end))*(abs(Q_d(end))*Z)./local_ro_d(end)/rho_r; W_p2 = W_p2*W_r;
        W_t = T_eff * (p_d(end)-pE)*(abs(Q_d(1))*Z)/local_ro_d(1)/rho_r; W_t = W_t*W_r;
    end
    W_p4 = 1/LP_eff * (pE - p_f(1))*(Q_f(1) *Z)./rho_E; W_p4 = W_p4*W_r;
    W_p1=0; W_p3=0;
end

%% Version(6)=3 --> SWRO-PRO hybrid system (with one ERD)
if version(6)==3
    % flow from M1 to ERD1
    J_M1 = J_E * (1*mix_M1);
    J_w_M1 = J_wE* (1*mix_M1);
    J_s_M1 = J_sE* (1*mix_M1);
    % flow from ERD1 to RO
    rho_ERD= V_m*(swro_local_ro_d(end) - rho_E/rho_r)+rho_E/rho_r;
    C_ERD= -ro_salt/rho_r*(rho_ERD-ro_water/rho_r)/(ro_water/rho_r*(rho_ERD-ro_salt/rho_r));
    J_ERD = J_M1;
    J_wERD= J_ERD/(C_ERD+1);
    J_sERD= J_ERD - J_wERD;
    qq=J_r/Q_r;
    % not affected by co / counter
    f_2 = ERD_fric * mix_density * (J_ERD*swro_b2/rho_ERD)*((J_d(end)*swro_b2)/(swro_local_ro_d(end)*swro_b1*swro_b2))^2;
    f_r= (J_r*p_r*swro_x_r/rho_r) / (J_r^3/rho_r^2/swro_x_r); % scaling factor between f and P*J*Z/rho  
    W_p1 = 1/HP_eff * (pE - P_d(1))*((J_E-J_M1) *swro_Z)/rho_E; W_p1 = W_p1*swro_W_r;
    W_p2 = 0; 
    W_p4= 1/LP_eff * (pE - p_f(1))*(Q_f(1)*Z)./local_ro_f(1)/rho_r; W_p4 = W_p4*W_r;  
    if version(6) == 0 % co-current
        f_1 = ERD_fric * mix_density * (qq*Q_d(1)*b2/local_ro_d(1))*((J_M1*swro_b2)/(rho_E/rho_r*swro_b1*swro_b2))^2; 
        pERD = rho_ERD*(ERD_eff*(P_d(end)*J_d(end)*swro_b2*(1-eta_ERD)/swro_local_ro_d(end)-p_d(1)*qq*Q_d(1)*b2/local_ro_d(1) - f_2/f_r) + pE*J_M1*swro_b2/rho_E*rho_r + f_1/f_r)/J_ERD/swro_b2;
        W_t = T_eff * (p_d(end)-pE)*(Q_d(end)*Z)/local_ro_d(end)/rho_r; W_t = W_t*W_r;
    else % counter-current
        f_1 = ERD_fric * mix_density * (qq*abs(Q_d(end))*b2/local_ro_d(end))*((J_M1*swro_b2)/(rho_E/rho_r*swro_b1*swro_b2))^2; 
        pERD = rho_ERD*(ERD_eff*(P_d(end)*J_d(end)*swro_b2*(1-eta_ERD)/swro_local_ro_d(end)-p_d(end)*qq*abs(Q_d(end))*b2/local_ro_d(end) - f_2/f_r) + pE*J_M1*swro_b2/rho_E*rho_r + f_1/f_r)/J_ERD/swro_b2;
        W_t = T_eff * (p_d(1)-pE)*(abs(Q_d(1))*Z)/local_ro_d(1)/rho_r; W_t = W_t*W_r;
    end
    norm_f1=min(abs(f_1-f_2),abs(f_2-f_1))/J_r/p_r/swro_x_r*rho_r; % dimensions: J_r*p_r*swro_x_r*/rho_r [W];
    W_p3 = 1/HP_eff * (pERD - P_d(1))*(J_ERD*swro_Z)/rho_ERD/rho_r; W_p3 = W_p3*swro_W_r;
end

%% Version(6)=4 --> SWRO-PRO hybrid system (with two ERDs)
if version(6)==4 
    % calculate concentration atand total flow at 2nd ERD
    if version(1)==0
        rho_d2= real((c_d(end)+1)./(c_d(end)/ro_salt + 1/ro_water));
    else
        rho_d2= real((c_d(1)+1)./(c_d(1)/ro_salt + 1/ro_water));
    end
    rho_ERD2 = V_m*(rho_d2 - rho_E)+rho_E;
    C_ERD2= -ro_salt*(rho_ERD2-ro_water)/(ro_water*(rho_ERD2-ro_salt));
    % at first ERD
    rho_d1= real((C_d(1)+1)./(C_d(1)/ro_salt + 1/ro_water));
    rho_ERD1= V_m*(rho_d1 - rho_ERD2) + rho_ERD2;
    C_ERD1= -ro_salt*(rho_ERD1-ro_water)/(ro_water*(rho_ERD1-ro_salt));
    % rescale density
    rho_d1=rho_d1/rho_r;
    rho_d2=rho_d2/rho_r;
    rho_ERD1=rho_ERD1/rho_r;
    rho_ERD2=rho_ERD2/rho_r;
    % seawater after passing 2nd ERD 
    J_E1 = J_d(1);
    J_ERD1 = J_E1*(1-mix_M1);
    J_wE1 = J_E1/(C_ERD2+1);
    J_wERD1= J_ERD1/(C_ERD1+1);
    % scaling factor
    qq=Q_r/J_r;
    % seawater before 2nd ERD
    J_E2 = J_d(1);
    J_wE2 = J_E1/(cE+1);
    if version(1)==0
        Q_exit = mix_M3*Q_d(end)*(1-eta_ERD);
        c_exit = (cE*J_wE2 + qq*Q_sd(end)*(1-eta_ERD) - C_ERD2*J_wE2)/(J_wE2 + qq*Q_wd(end)*(1-eta_ERD) - J_wE2);
    else
        Q_exit = (1-mix_M3)*abs(Q_d(1))*(1-eta_ERD);
        c_exit= (cE*J_wE2 + qq*abs(Q_sd(1))*(1-eta_ERD) - C_ERD2*J_wE2)/(J_wE2 + qq*abs(Q_wd(1))*(1-eta_ERD) - J_wE2);
    end
    rho_exit= (ro_salt*(c_exit+1))/(c_exit*ro_water+ro_salt); p_exit=pE;
    % Calculation for 2nd ERD ---------------------------------------------
    f_r= (J_r*p_r*swro_x_r/rho_r) / (J_r^3/rho_r^2/swro_x_r); % scaling factor between f and P*J*Z/rho 
    W_p2 = 0;
    W_p4= 1/LP_eff * (pE - p_f(1))*(Q_f(1)*Z)./local_ro_f(1)/rho_r; W_p4 = W_p4*W_r;
    if version(6) == 0 % co-current
        f_2 = ERD_fric * mix_density * (J_E1*swro_b2/rho_ERD2)*((qq*(1-mix_M3)*Q_d(end)*swro_b2)/(local_ro_d(end)*swro_b1*swro_b2))^2;
        f_1 = ERD_fric * mix_density * (qq*Q_exit*swro_b2/rho_exit)*((J_E2*swro_b2)/(rho_E/rho_r*swro_b1*swro_b2))^2; 
        pERD2 = rho_ERD2*(ERD_eff*(p_d(end)*qq*(1-mix_M3)*Q_d(end)*swro_b2*(1-eta_ERD)/local_ro_d(end)-pE*qq*Q_exit*swro_b2/rho_exit - f_2/f_r) + pE*J_E2*swro_b2/rho_E*rho_r + f_1/f_r)/J_E1/swro_b2;
        W_t = T_eff * (p_d(end)-pE)*(mix_M3*Q_d(end)*Z)/local_ro_d(end)/rho_r; W_t = W_t*W_r;
    else % counter-current
        f_2 = ERD_fric * mix_density * (J_E1*swro_b2/rho_ERD2)*((qq*(1-mix_M3)*abs(Q_d(1))*swro_b2)/(local_ro_d(1)*swro_b1*swro_b2))^2;
        f_1 = ERD_fric * mix_density * (qq*Q_exit*swro_b2/rho_exit)*((J_E2*swro_b2)/(rho_E/rho_r*swro_b1*swro_b2))^2; 
        pERD2 = rho_ERD2*(ERD_eff*(p_d(1)*qq*(1-mix_M3)*abs(Q_d(1))*swro_b2*(1-eta_ERD)/local_ro_d(1)-pE*qq*abs(Q_exit)*swro_b2/rho_exit - f_2/f_r) + pE*J_E2*swro_b2/rho_E*rho_r + f_1/f_r)/J_E2/swro_b2;
        W_t = T_eff * (p_d(1)-pE)*(mix_M3*abs(Q_d(1))*Z)/local_ro_d(1)/rho_r; W_t = W_t*W_r;
    end
    norm_f3=min(abs(f_1-f_2),abs(f_2-f_1))/J_r/p_r/swro_x_r*rho_r; % dimensions: J_r*p_r*swro_x_r*/rho_r [W];
    % Calculation for 1st ERD ---------------------------------------------
    f_2 = ERD_fric * mix_density * (J_ERD1*swro_b2/rho_ERD1)*((J_d(end)*swro_b2)/(swro_local_ro_d(end)*swro_b1*swro_b2))^2;
    if version(6) == 0 % co-current
        f_1 = ERD_fric * mix_density * (qq*Q_d(1)*b2/local_ro_d(1))*((J_E1*(1-mix_M1)*swro_b2)/(rho_ERD2*swro_b1*swro_b2))^2; 
        pERD = rho_ERD1*(ERD_eff*(P_d(end)*J_d(end)*swro_b2*(1-eta_ERD)/swro_local_ro_d(end)-p_d(1)*qq*Q_d(1)*b2/local_ro_d(1) - f_2/f_r) + pERD2*J_E1*(1-mix_M1)*swro_b2/rho_ERD2 + f_1/f_r)/J_ERD1/swro_b2;
    else % counter-current
        f_1 = ERD_fric * mix_density * (qq*abs(Q_d(end))*b2/local_ro_d(end))*((J_E1*(1-mix_M1)*swro_b2)/(rho_ERD2*swro_b1*swro_b2))^2; 
        pERD = rho_ERD1*(ERD_eff*(P_d(end)*J_d(end)*swro_b2*(1-eta_ERD)/swro_local_ro_d(end)-p_d(end)*qq*abs(Q_d(end))*b2/local_ro_d(end) - f_2/f_r) + pERD2*J_E1*(1-mix_M1)*swro_b2/rho_ERD2 + f_1/f_r)/J_ERD1/swro_b2;
    end
    norm_f1=min(abs(f_1-f_2),abs(f_2-f_1))/J_r/p_r/swro_x_r*rho_r; % dimensions: J_r*p_r*swro_x_r*/rho_r [W]; 
    W_p3 = 1/HP_eff * (pERD-P_d(1))*(J_ERD1*swro_Z)/rho_ERD1/rho_r; W_p3 = W_p3*swro_W_r;
    W_p1 = 1/HP_eff * (pERD2-P_d(1))*(J_E1*mix_M1*swro_Z)/rho_ERD2/rho_r; W_p1 = W_p1*swro_W_r;
end

%% final output:
W_net= W_p1 + W_p2 + W_p3 + W_p4 + W_t; % in [W]

%% Calculate the final output
SEC_net = W_net*(swro_local_ro_f(end)*rho_r)./(J_f(end)*J_r*swro_Z)/1000/3600; % in [kWh/m^3]  
    if SEC_net > 0 && contains('solfig',obj)==1 && version(6)~=2
        fprintf(2,' \nWaring: SEC_net is positive! \n');
    end
FW = (J_wf(end)*J_r*swro_Z/(swro_local_ro_f(end)*rho_r))*3600; % in [m^3/h] 
    if FW < 0 && contains('solfig',obj)==1 && version(6)~=2
        fprintf(2,' \nWaring: Freshwater production is negative! \n');
    end
Rev= pw*FW + pe*SEC_net*FW; %in [$/h]
SWRO_Recovery =(J_wf(end)./J_d(1))*100;     % in [%]
    if SWRO_Recovery > 100 && contains('solfig',obj)==1 && version(6)~=2
        fprintf(2,' \nWaring: SWRO recovery is greater than 100 %% \n');
    end
    if SWRO_Recovery < 0 && contains('solfig',obj)==1 && version(6)~=2
        fprintf(2,' \nWaring: SWRO recovery is negative! \n');
    end
PRO_Recovery = NaN; % in [%]
    if version(6) > 1
        PRO_Recovery = (1-Q_f(end)/Q_f(1))*100;
        if PRO_Recovery > 100 && contains('solfig',obj)==1
            fprintf(2,' \nWaring: PRO recovery is greater than 100 %% \n');
        end
        if PRO_Recovery < 0 && contains('solfig',obj)==1
            fprintf(2,' \nWaring: PRO recovery is negative! \n');
        end
    end
RO_inflow = swro_Z*J_d(1)*J_r/swro_local_ro_d(1)/rho_r; % in [m^3/s]
Permeate_outflow = swro_Z*J_f(end)*J_r/swro_local_ro_f(end)/rho_r; % in [m^3/s]
Wastewater_inflow = NaN; % in [m^3/s]
    if version(6) > 1
        Wastewater_inflow = Z*Q_f(1)*Q_r/local_ro_f(1)/rho_r; % in [m^3/s]
    end
C_permeate = 10000*C_f(end); % in [ppm]
C_brine = 100*C_d(end); % in [%]
C_dilluted = NaN; % in [%]
    if version(6) > 1
        if version(1)==1
            C_brine=c_d(end);
            C_dilluted =c_d(1);
        else
            C_brine=c_d(1);
            C_dilluted =c_d(end);
        end
    end

if version(6)==0; mix_M1=NaN; norm_f1=NaN; end
if version(6)==2; mix_M1=NaN; norm_f1=NaN; end
if version(6)<4; mix_M3=NaN; norm_f3=NaN; end
if version(6)==2
    SEC_net=NaN;FW=NaN;Rev=NaN;SWRO_Recovery=NaN;C_permeate=NaN;mix_M1=NaN;
    RO_inflow=W_net./(Z*L); 
    if version(6)==0; Permeate_outflow=p_d(end); else; Permeate_outflow=p_d(1); end
end

output1 = [SEC_net, FW, Rev, SWRO_Recovery, PRO_Recovery, ...
           RO_inflow, Permeate_outflow, Wastewater_inflow, C_permeate, C_brine, C_dilluted, ...
           W_net, W_p1, W_p2,  W_p3, W_p4, W_t, ...
           mix_M1, norm_f1, mix_M3, norm_f3]; % length(output1) = 5+6+6+4 = 21

if version(6)==2
%    PD_net = W_net./(Z*L); %in [W/m^2]
%    disp(['PD_net     = ', num2str(PD_net), ' [kWh/m^2]'])
%    if version(1)==0 
%        SE_net = W_net./(Q_f(1)*Q_r*Z/(rho_r*local_ro_f(1)) + Q_d(1)*Q_r*Z/(rho_r*local_ro_d(1)))/1000/3600; %in [kWh/m^3] 
%        SE_f   = W_net./(Q_f(1)*Q_r*Z/(rho_r*local_ro_f(1)))/1000/3600; %in [kWh/m^3]
%    else 
%        SE_net = W_net./(Q_f(1)*Q_r*Z/(rho_r*local_ro_f(1)) + abs(Q_d(end))*Q_r*Z/(rho_r*local_ro_d(end)))/1000/3600; %in [kWh/m^3] 
%        SE_f   = W_net./(Q_f(1)*Q_r*Z/(rho_r*local_ro_f(1)))/1000/3600; %in [kWh/m^3]
%    end 
end


%% output, if BVP-solver fails
catch
    fprintf(2,' \nERROR: bvp5c failed! \n');
    sol.stats=struct('nmeshpoints', NaN, 'maxerr', 1e20, 'nODEevals', NaN, 'nBCevals', NaN);
end
output2 = sol.stats;

%% Output of objective function
switch (obj)
        case 'SEC' % for SEC_net maximization  -->  output1 = -SEC_net
            if sol.stats.maxerr > option_BVP
            output1 = 20;
            else; output1 = -SEC_net;
            end 
        case 'FW' % for FW maximization  -->  output1 = -FW
            if sol.stats.maxerr > option_BVP
            output1 = 20; 
            else; output1 = -FW;
            end 
        case 'Rev' % for Rev maximization  -->  output1 = -Rev
            if sol.stats.maxerr > option_BVP
            output1 = 0;
            else; output1 = -Rev;
            end
        case 'PD_net' % only PRO
            if sol.stats.maxerr > option_BVP
            output1 = -20;
            else; output1 = -RO_inflow;
            end
        case 'REC' % only PRO
            if sol.stats.maxerr > option_BVP
            output1 = -20;
            else; output1 = -PRO_Recovery;
            end
        case 'Pareto' % SEC_net and FW maximization  -->  output1 = [-SEC_net, -FW]
            if sol.stats.maxerr > option_BVP    
            output1 = NaN(1,2); 
            else
                output1 = [-SEC_net, -FW];
                if FW < 0.1
                    output1 = NaN(1,2); 
                end
            end 
        case 'sol' 
            if sol.stats.maxerr > option_BVP
            output1 = NaN(1,21);
            end
    case 'fig' % case ends after the figures
            if sol.stats.maxerr > option_BVP
            output1 = NaN(1,21);
            else

%% figure 1
if sol.stats.maxerr ==1000
fprintf(2,' \nERROR: bvp5c could not satisfy the relative error tolerance ---> no figures could be displayed\n');
else
close all
lw=1.5; %Linewidth for all figures
if fig(1) == 1
f=figure(1); %
f.Position = [277.6667 275 1100 1000];
x2=x(2:end); x3=x(2:end-1);
% total mass flows
subplot(3,2,1);lc='#0072bD';rc='#77AC30';
plot(x*swro_L, J_d*J_r,'Color', lc,'LineWidth',lw);xlabel('x [m]','Fontsize',10); ylabel('[kg/sm]','Fontsize',10); ay=gca; ay.YAxis.Exponent = 0;xlim([0 swro_L]); hold on
yyaxis right
plot(x*swro_L, J_f*J_r,'Color', rc,'LineWidth',lw); ylabel('[kg/sm]','Fontsize',10); legend('J_{d}^{RO}(x)','J_{f}^{RO}(x)','Location','East');xlim([0 swro_L]); ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
% concentrations
subplot(3,2,3);lc='#0072bD';rc='#77AC30';
plot(x*swro_L, 100*C_d*C_r,'Color', lc,'LineWidth',lw);xlabel('x [m]','Fontsize',10); ylabel('[%]','Fontsize',10); ay=gca; ay.YAxis.Exponent = 0;hold on
yyaxis right; C_f=C_f(2:end);
%plot(x2*swro_L, 100*C_f*C_r,'Color', rc,'LineWidth',lw); ylabel('[%]','Fontsize',10); legend('C_d^{RO}(x)','C_f^{RO}(x)','Location','best');xlim([0 swro_L]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
plot(x2*swro_L, 1e4*C_f*C_r,'Color', rc,'LineWidth',lw); ylabel('[ppm]','Fontsize',10); legend('C_d^{RO}(x)','C_f^{RO}(x)','Location','best');xlim([0 swro_L]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
% Pressures
subplot(3,2,5);lc='#0072bD';rc='#77AC30';
plot(x*swro_L, P_d*p_r,'Color', lc,'LineWidth',lw);xlabel('x [m]','Fontsize',10); ylabel('[Pa]','Fontsize',10); ay=gca; ay.YAxis.Exponent = 5; hold on
yyaxis right
plot(x*swro_L, P_f*p_r,'Color', rc,'LineWidth',lw); ylabel('[Pa]','Fontsize',10); legend('P_d^{RO}(x)','P_f^{RO}(x)','Location','best');xlim([0 swro_L]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
% total mass flows
subplot(3,2,2);lc='#0072bD';rc='#77AC30';
plot(x*L,  Q_d*Q_r,'Color', lc,'LineWidth',lw);xlabel('x [m]','Fontsize',10); ylabel('[kg/sm]','Fontsize',10); ay=gca; ay.YAxis.Exponent = 0;xlim([0 swro_L]); hold on
yyaxis right
plot(x*L, Q_f*Q_r,'Color', rc,'LineWidth',lw); ylabel('[kg/sm]','Fontsize',10); legend('J_{d}^{PRO}(x)','J_{f}^{PRO}(x)','Location','East');xlim([0 L]); ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
% concentrations
subplot(3,2,4);lc='#0072bD';rc='#77AC30';
plot(x*L, 100*c_d*C_r,'Color', lc,'LineWidth',lw);xlabel('x [m]','Fontsize',10); ylabel('[%]','Fontsize',10); ay=gca; ay.YAxis(1).Exponent = 0;hold on
yyaxis right; c_f=c_f(2:end);
plot(x2*L, 100*abs(c_f)*C_r,'Color', rc,'LineWidth',lw); ylabel('[%]','Fontsize',10); legend('C_d^{PRO}(x)','C_f^{PRO}(x)','Location','best');xlim([0 L]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;ay.YAxis(2).Exponent = 0;
% Pressures
subplot(3,2,6);lc='#0072bD';rc='#77AC30';
plot(x*L, p_d*p_r,'Color', lc,'LineWidth',lw);xlabel('x [m]','Fontsize',10); ylabel('[Pa]','Fontsize',10);ay=gca; ay.YAxis.Exponent = 5; hold on
yyaxis right
plot(x*L, p_f*p_r,'Color', rc,'LineWidth',lw); ylabel('[Pa]','Fontsize',10); legend('P_d^{PRO}(x)','P_f^{PRO}(x)','Location','best');xlim([0 L]);ay.YAxis(1).Color = lc; ay.YAxis(2).Exponent = 5; ay.YAxis(2).Color = rc;
end
%% figure 2
if fig(2) == 1
f=figure(2); %
f.Position = [1.2977e+03 635.6667 1.0987e+03 639.3333];
% Permeate flows 
subplot(2,2,1); lc='k';rc='#b81414'; J_cross2=J_cross(2:end);J_sin2=J_sin(2:end);
plot(x2*swro_L, J_cross2*p_r*swro_alpha, 'Color', lc, 'LineWidth',lw);xlabel('x [m]','Fontsize',10);  ylabel('[kg/sm^2]','Fontsize',10);ay=gca; hold on
yyaxis right
plot(x2*swro_L, J_sin2*p_r*swro_alpha, 'Color', rc, 'LineWidth',lw); ylabel('[kg/sm^2]','Fontsize',10); legend('J_{w,in}^{RO}(x)','J_{s,in}^{RO}(x)','Location','best');xlim([0 swro_L]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;hold on
% Osmotic/ Hydraulic pressure difference
subplot(2,2,3);lc='k';rc='#b81414'; osm_diff=swro_p_osm_d(2:end)-swro_p_osm_f(2:end); 
a=min([(P_d-P_f).*p_r; osm_diff.*p_r]);b=max([(P_d-P_f).*p_r; osm_diff.*p_r]);
plot(x*swro_L, (P_d-P_f)*p_r,'Color', lc,'LineWidth',lw);xlabel('x [m]','Fontsize',10); ylabel('[Pa]','Fontsize',10); ay=gca; ay.YAxis.Exponent = 6;ylim([0.95*a 1.05*b]);hold on
yyaxis right; 
plot(x2*swro_L, osm_diff*p_r,'Color', rc,'LineWidth',lw);ylim([0.95*a 1.05*b]); ylabel('[Pa]','Fontsize',10); legend('\Delta P^{RO}(x)','\Delta \pi^{RO}(x)','Location','best');xlim([0 swro_L]);ay=gca;ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;hold on
% Permeate flows 
subplot(2,2,2); lc='k';rc='#b81414'; Q_cross2=Q_cross(2:end-1);Q_sin2=Q_sin(2:end-1);
plot(x3*L, -Q_cross2*p_r*alpha, 'Color', lc, 'LineWidth',lw);xlabel('x [m]','Fontsize',10);  ylabel('[kg/sm^2]','Fontsize',10);ay=gca; hold on
yyaxis right
plot(x3*L, Q_sin2*p_r*alpha, 'Color', rc, 'LineWidth',lw); ylabel('[kg/sm^2]','Fontsize',10); legend('J_{w,in}^{PRO}(x)','J_{s,in}^{PRO}(x)','Location','West');xlim([0 L]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
% Osmotic/ Hydraulic pressure difference
subplot(2,2,4);lc='k';rc='#b81414';osm_diff=p_osm_d(2:end)-p_osm_f(2:end);
a=min([(p_d-p_f).*p_r; osm_diff.*p_r]);b=max([(p_d-p_f).*p_r; osm_diff.*p_r]);
plot(x*L, (p_d-p_f)*p_r,'Color', lc,'LineWidth',lw);xlabel('x [m]','Fontsize',10); ylabel('[Pa]','Fontsize',10); ay=gca; ay.YAxis.Exponent = 6;ylim([0 1.05*b]);hold on
yyaxis right; 
plot(x2*L, osm_diff*p_r,'Color', rc,'LineWidth',lw); ylabel('[Pa]','Fontsize',10); legend('\Delta P^{PRO}(x)','\Delta \pi^{PRO}(x)','Location','NorthEast');xlim([0 L]);ylim([0 1.05*b]); ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc; hold on
end

%% figure 3
if fig(3) == 1 && version(6) > 0
f=figure(3);  % ERD quantities
f.Position = [637.6667 210.3333 900 347.3334]; tiledlayout(1,4); 
rx=0; %rotation and shift of texts
if version(6) == 1
vals1 = [pE*p_r; pERD*p_r; 0; P_d(end)*p_r; pE*p_r];
vals2 = 100*[cE; C_ERD; 0; C_d(end); c_exit];
vals3 = [J_M1*J_r; J_ERD*J_r; 0; J_d(end)*J_r; J_exit*J_r];
end
if version(6) == 3
vals1 = [pE*p_r; pERD*p_r; 0; P_d(end)*p_r; max(p_d(1),p_d(end))*p_r];
vals2 = 100*[cE; C_ERD; 0; C_d(end); max(c_d(1),c_d(end))];
vals3 = [J_M1*J_r; J_ERD*J_r; 0; J_d(end)*J_r; min(abs(Q_d(1))*Q_r, abs(Q_d(end))*Q_r)];
end
if version(6) == 4
vals1 = [pERD2*p_r; pERD*p_r; 0; P_d(end)*p_r; max(p_d(1),p_d(end))*p_r];
vals2 = 100*[C_ERD2; C_ERD1; 0; C_d(end); max(c_d(1),c_d(end))];
vals3 = [J_E1*(1-mix_M1)*J_r; J_ERD1*J_r; 0; J_d(end)*J_r; min(abs(Q_d(1))*Q_r, abs(Q_d(end))*Q_r)];
end
nexttile
b = bar(1,vals1, 'FaceColor', 'b'); b(4).FaceColor = "#77AC30"; b(5).FaceColor = [.2 .6 .5]; b(1).FaceColor = '#0072bD';
xticklabels({'[Pa]'}); ylim([0 1.1*max(vals1)]);ay=gca;ay.YAxis(1).Exponent=5;
text(b(1).XEndPoints,b(1).YEndPoints,"P_{s}^{1;in} ",'HorizontalAlignment','center','VerticalAlignment','bottom');
text(b(2).XEndPoints+rx ,b(2).YEndPoints,"P_{s}^{1;out}",'HorizontalAlignment','center','VerticalAlignment','bottom'); 
text(b(4).XEndPoints,b(4).YEndPoints,"P_{b}^{1;in}",'HorizontalAlignment','center','VerticalAlignment','bottom'); 
text(b(5).XEndPoints,b(5).YEndPoints,"   P_{b}^{1;out}",'HorizontalAlignment','center','VerticalAlignment','bottom'); 
nexttile
b = bar(1,vals2, 'FaceColor', 'b'); b(4).FaceColor = "#77AC30"; b(5).FaceColor = [.2 .6 .5]; b(1).FaceColor = '#0072bD';
xticklabels({'[%]'}); ylim([0 1.12*max(vals2)]);ay=gca;ay.YAxis(1).Exponent=0;
text(b(1).XEndPoints,b(1).YEndPoints,"C_{s}^{1;in}  ",'HorizontalAlignment','center','VerticalAlignment','bottom');
text(b(2).XEndPoints+rx ,b(2).YEndPoints,"   C_{s}^{1;out}",'HorizontalAlignment','center','VerticalAlignment','bottom');
text(b(4).XEndPoints,b(4).YEndPoints,"C_b^{1;in}  ",'HorizontalAlignment','center','VerticalAlignment','bottom');
text(b(5).XEndPoints,b(5).YEndPoints,"  C_{b}^{1;out}",'HorizontalAlignment','center','VerticalAlignment','bottom');
nexttile
b = bar(1,vals3, 'FaceColor', 'b'); b(4).FaceColor = "#77AC30"; b(5).FaceColor = [.2 .6 .5]; b(1).FaceColor = '#0072bD';
xticklabels({'[kg/ms]'}); ylim([0 1.1*max(vals3)]);ay=gca;ay.YAxis(1).Exponent=0;
text(b(1).XEndPoints,b(1).YEndPoints,"J_{s}^{1;in} ",'HorizontalAlignment','center','VerticalAlignment','bottom');
text(b(2).XEndPoints+rx ,b(2).YEndPoints," J_{s}^{1;out}",'HorizontalAlignment','center','VerticalAlignment','bottom'); 
text(b(4).XEndPoints ,b(4).YEndPoints,"J_{b}^{1;in} ",'HorizontalAlignment','center','VerticalAlignment','bottom'); 
text(b(5).XEndPoints,b(5).YEndPoints," J_{b}^{1;out}",'HorizontalAlignment','center','VerticalAlignment','bottom');
nexttile
lgd = legend([b(1) b(2) b(5) b(4)],'Seawater inlet','Seawater outlet','Brine inlet','Brine outlet' , 'Location', 'EastOutside');
lgd.Layout.Tile = 4; axis off ; 
end
%% figure 4
if fig(3) == 1 && version(6) ==4
f=figure(4);  % ERD quantities
f.Position = [1497 210.3333 900 347.3334]; tiledlayout(1,4); 
rot=0; rx=0; %rotation and shift of texts
vals1 = [pE*p_r; pERD2*p_r; 0; min(p_d(1),p_d(end))*p_r; p_exit*p_r];
vals2 = 100*[cE; C_ERD2; 0; min(c_d(1),c_d(end)); c_exit];
vals3 = [J_E2*J_r; J_E1*J_r; 0; (1-mix_M3)*max(abs(Q_d(1)),abs(Q_d(end)))*Q_r; abs(Q_exit)*Q_r];
nexttile
b = bar(1,vals1, 'FaceColor', 'b'); b(4).FaceColor = "#77AC30"; b(5).FaceColor = [.2 .6 .5]; b(1).FaceColor = '#0072bD';
xticklabels({'[Pa]'}); ylim([0 1.1*max(vals1)]);ay=gca;ay.YAxis(1).Exponent=5;
t=text(b(1).XEndPoints,b(1).YEndPoints,"P_E",'HorizontalAlignment','center','VerticalAlignment','bottom');
t=text(b(2).XEndPoints+rx ,b(2).YEndPoints,"P_{s}^{2;out}",'HorizontalAlignment','center','VerticalAlignment','bottom'); 
t=text(b(4).XEndPoints,b(4).YEndPoints,"P_b^{2;in}",'HorizontalAlignment','center','VerticalAlignment','bottom');
t=text(b(5).XEndPoints,b(5).YEndPoints,"P_E",'HorizontalAlignment','center','VerticalAlignment','bottom'); 
nexttile
b = bar(1,vals2, 'FaceColor', 'b'); b(4).FaceColor = "#77AC30"; b(5).FaceColor = [.2 .6 .5]; b(1).FaceColor = '#0072bD';
xticklabels({'[%]'}); ylim([0 1.12*max(vals2)]);ay=gca;ay.YAxis(1).Exponent=0;
t= text(b(1).XEndPoints,b(1).YEndPoints,"C_E",'HorizontalAlignment','center','VerticalAlignment','bottom'); t.Rotation = rot;
t=text(b(2).XEndPoints+rx ,b(2).YEndPoints,"  C_{s}^{2;out}",'HorizontalAlignment','center','VerticalAlignment','bottom'); t.Rotation = rot;
t=text(b(4).XEndPoints,b(4).YEndPoints,"C_{b}^{2;in}  ",'HorizontalAlignment','center','VerticalAlignment','bottom'); t.Rotation = rot; 
t=text(b(5).XEndPoints,b(5).YEndPoints,"  C_{b}^{2;out}",'HorizontalAlignment','center','VerticalAlignment','bottom'); t.Rotation = rot; 
nexttile
b = bar(1,vals3, 'FaceColor', 'b'); b(4).FaceColor = "#77AC30"; b(5).FaceColor = [.2 .6 .5]; b(1).FaceColor = '#0072bD';
xticklabels({'[kg/ms]'}); ylim([0 1.1*max(vals3)]);ay=gca;ay.YAxis(1).Exponent=0;
t=text(b(1).XEndPoints,b(1).YEndPoints,"J_{s}^{2;in} ",'HorizontalAlignment','center','VerticalAlignment','bottom'); t.Rotation = rot;
t=text(b(2).XEndPoints+rx ,b(2).YEndPoints,"  J_s^{2;out}",'HorizontalAlignment','center','VerticalAlignment','bottom'); t.Rotation = rot;
t=text(b(4).XEndPoints,b(4).YEndPoints,"J_b^{2;in}  ",'HorizontalAlignment','center','VerticalAlignment','bottom'); t.Rotation = rot;
t=text(b(5).XEndPoints,b(5).YEndPoints," J_b^{2;out}",'HorizontalAlignment','center','VerticalAlignment','bottom'); t.Rotation = rot;
nexttile
lgd = legend([b(1) b(2) b(5) b(4)],'Seawater inlet','Seawater outlet','Diluted brine inlet','Diluted brine outlet' , 'Location', 'EastOutside');
lgd.Layout.Tile = 4; axis off ; 
end
end
end
end