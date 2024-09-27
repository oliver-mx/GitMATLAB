function [output1, output2, output3] = fun_unscaled(input,option_data,obj,option_mesh,option_BVP)
%%  fun_unscaled   Solves the unscaled BVP for the SWRO plant using bvp5c
%
%       fun(input,option_data,obj,option_mesh,option_BVP) will solve the
%       nonlinear BVP of the selected SWRO-PRO system. In "option_data" 
%       all Informations of the system are given with exception of the 
%       additional "input". 
%
%   Input:
%       input         -   Input for the Data function
%       option_data   -   Selects Data function used: data(input)
%       obj           -   'SEC', 'FW', 'Rev', 'Pareto', 'PD_net', 'SE_net', 'SE_f', 'sol';
%       option_mesh   -   NMax of bvp5c
%       option_BVP    -   RelTol of bvp5c
%
%   Output:
%       output1 = [SEC_net, FW, Rev, SWRO_Rec, PRO_Rec, PD_net, SE_net];
%       output2 = [W_net, -W_p1, -W_p2,  -W_p3, -W_p4, W_t];
%       figures may be displayed if obj='sol'
%
%   Example:
%       x0=1;       
%       D=.11;      
%       obj='sol';  
%       mesh=2e4;
%       tol=1e-6;  
%       [a,b]=fun(x0,D,obj,mesh,tol),
%

%% Read data
%
if option_data == 0; DATA = @(x)Test_01_data(input); end
%
if option_data == 1; DATA = @(x)Pareto_1_data(input); end
if option_data == 2; DATA = @(x)Pareto_2_data(input); end
if option_data == 3; DATA = @(x)Pareto_3_data(input); end
%

[H,Z,swro_Z,ro_water,ro_salt,Mw,Ms,Rw,T0,eta,sigma,p_r,rho_r,C_r,swro_L,swro_alpha,swro_R,swro_KK,swro_x_r,swro_b1,swro_b2,J_r,swro_gamma,swro_gamma2,swro_W_r,L,alpha,R,KK,x_r,b1,b2,Q_r,gamma,gamma2,W_r,cE,pE,rho_E,J_sf_0,J_wf_0,Pd_0,Pd_L,Pf_L,Q_sf_0,pd_0,pf_0,pd_L,pf_L,HP_eff,LP_eff,T_eff,V_m,ERD_eff,ERD_fric,A_ERD,eta_ERD,mix_density,pw,pe,swro_beta_fix,beta_fix,mixer_ERD,version,fig,swro_KF,swro_KD, KF, KD] ...
    = DATA(1);

%% Set options
ode_options = bvpset('Stats','off','NMax', option_mesh, 'RelTol', option_BVP);      

n=25;
x = linspace(0,swro_L,n);    
  
%% Create initial guess for the BVP
J_wd_0 = 37.36;          % guess for J_wd(0)
J_wf_0 = 35;             % should be close to 0
Pf_0   = 1.1e5*Pf_L;                % guess for P_f(0)

y_init = [cE; J_wd_0; J_sf_0; J_wf_0; Pd_0; Pf_0];
solinit = bvpinit(x,y_init);

%% Solve the BVP
if version(6) == 0; sol = bvp5c(@(x,J_p)Unscaled_ODEsystem(x, J_p, DATA), @(ya,yb)Unscaled_Boundary2(ya,yb, DATA),solinit,ode_options); end
if version(6) == 1; sol = bvp5c(@(x,J_p)Unscaled_ODEsystem(x, J_p, DATA), @(ya,yb)Unscaled_Boundary1(ya,yb, DATA),solinit,ode_options); end
if version(6) > 1; fprintf(2,' \nERROR: Unscaled model only supports SWRO with/without ERD! \n'); end

%% Evaluate the solution
if any([0,0,0,0] ~= fig)
    n=100; % higher resolution plots
end
x=linspace(0,swro_L,n);
y = deval(sol,x); Y = y'; Y = real(Y);
%display(Y(:,4)) 
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
swro_p_osm_d = ro_water*Rw*T0*log(1 + 2*Mw*Y(:,1)/Ms);
swro_p_osm_f = ro_water*Rw*T0*log(1 + 2*Mw*Y(:,3)/Ms./Y(:,4));
% local density
swro_local_ro_d = (Y(:,1) + 1)./(Y(:,1)./ro_salt + ones(n,1)./ro_water);
swro_local_ro_f = (Y(:,3) + Y(:,4))./(Y(:,3)./ro_salt + Y(:,4)./ro_water);
    % Salt permeability
    if version(4) == 0
        swro_beta = zeros(1,n);
    elseif version(2) == 0
        swro_beta = swro_beta_fix * J_r / swro_x_r.*ones(n,1);
    else
        swro_beta = (1 - swro_R) .* ((J_p(5) - J_p(6)) - sigma .* (swro_p_osm_d - swro_p_osm_f)) ./ swro_R;
    end
    % Water Permeate flux J_win(x)
    if version(4) == 0 % ideal
        J_cross = swro_alpha.*((Y(:,5) - Y(:,6)) - sigma*(swro_p_osm_d - swro_p_osm_f));
    elseif version(7) == 0 % ICP
        J_cross = swro_alpha .* ((swro_p_osm_d - swro_p_osm_f) - (Y(:,5) - Y(:,6)) - (Y(:,5) - Y(:,6)) .* swro_KK .* swro_beta) ./ (1 + swro_KK .* (swro_beta + swro_alpha .* swro_p_osm_f));
    else % ICP+ECP
        J_cross = swro_alpha .* ((Y(:,5) - Y(:,6)) .* (1 + swro_beta .* (1 ./ swro_KD + swro_KK + 1 ./ swro_KF)) - sigma .* (swro_p_osm_d - swro_p_osm_f)) ./ (1 + swro_beta .* (1 ./ swro_KD + swro_KK + 1 ./ swro_KF) - swro_alpha .* sigma .* (swro_p_osm_d ./ swro_KD + swro_p_osm_f .* swro_KK + swro_p_osm_f ./ swro_KF));
    end
    % Salt Permeate J_sin(x)
    if version(7) == 1 % ICP+ECP
        J_sin = swro_beta .* ((Y(:,1) - Y(:,3) ./ Y(:,4)) - Y(:,1) .* J_cross ./ swro_KD - (Y(:,3) ./ Y(:,4)) .* J_cross .* (swro_KK + 1 ./ swro_KF)) ./ (1 + swro_beta .* (1 ./ swro_KD + swro_KK + 1 ./ swro_KF));
    else
        J_sin=swro_beta_fix*J_r/swro_x_r.*(C_d-C_f);
    end
%% Freswater production and recovery rates
Freshwater = J_wf(end); %in [kg/s] 
FW = (Freshwater*swro_Z/(swro_local_ro_f(end)))*3600; %in [m^3/h] 
SWRO_Recovery =(J_f(end)./J_d(1))*100;     % SWRO freshwater recovery [%]

% warning flags:
    if FW < 0 && strcmp(obj,'sol')==1
        fprintf(2,' \nERROR: Waring: Freshwater production is negative! \n');
    end
    if SWRO_Recovery > 100 && strcmp(obj,'sol')==1
        fprintf(2,' \nERROR: Waring: SWRO recovery is greater than 100 %% \n');
    end
    if SWRO_Recovery < 0 && strcmp(obj,'sol')==1
        fprintf(2,' \nERROR: Waring: SWRO recovery is negative! \n');
    end

%% Output of System
J_E = (C_d(1)*J_wd(1)+J_wd(1))/2; 
W_p1 = 1/HP_eff * (P_d(1)-1e5)*(J_E *swro_Z)/rho_E;
W_p2 = 0; 

%% Version(6)=0 -->  only SWRO (no ERD)
W_p3=W_p1; W_p4=0; W_t=0; 
    
%% final output:
W_net= - W_p1 - W_p2 - W_p3 - W_p4 + W_t; % in [W]

SEC_net = W_net*(swro_local_ro_f(end))./(J_f(end)*swro_Z)/1000/3600; %in [kWh/m^3]  
    if SEC_net > 0 && strcmp(obj,'sol')==1
        fprintf(2,' \nERROR: Waring: SEC_net is positive! \n');
    end
Rev= pw*FW + pe*SEC_net*FW; %in [$/h]

output1=[SEC_net, FW, Rev, SWRO_Recovery, NaN];
output2=[W_net, -W_p1, -W_p2,  -W_p3, -W_p4, W_t];
output3=[swro_Z*J_d(1)/swro_local_ro_d(1), swro_Z*J_f(end)/swro_local_ro_f(end), SWRO_Recovery, 10000*C_f(end)];

%% figure 1
if sol.stats.maxerr ==1000
fprintf(2,' \nERROR: bvp5c could not satisfy the relative error tolerance ---> no figures could be displayed\n');
else
close all
lw=1.5; %Linewidth for all figures
if fig(1) == 1
f=figure(1); %
f.Position = [1839 227 1100 1000];
x2=x(2:end);
% total mass flows
subplot(3,2,1);lc='#0072bD';rc='#77AC30';
plot(x, J_d,'Color', lc,'LineWidth',lw);xlabel('x [m]','Fontsize',10); ylabel('[kg/sm]','Fontsize',10); ay=gca; ay.YAxis.Exponent = 0;xlim([0 swro_L]); hold on
yyaxis right
plot(x, J_f,'Color', rc,'LineWidth',lw); ylabel('[kg/sm]','Fontsize',10); legend('J_{d}^{RO}(x)','J_{f}^{RO}(x)','Location','East');xlim([0 swro_L]); ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
% concentrations
subplot(3,2,3);lc='#0072bD';rc='#77AC30';
plot(x, 100*C_d*C_r,'Color', lc,'LineWidth',lw);xlabel('x [m]','Fontsize',10); ylabel('[%]','Fontsize',10); ay=gca; ay.YAxis.Exponent = 0;hold on
yyaxis right; C_f=C_f(2:end);
plot(x2, 100*C_f*C_r,'Color', rc,'LineWidth',lw); ylabel('[%]','Fontsize',10); legend('C_d^{RO}(x)','C_f^{RO}(x)','Location','best');xlim([0 swro_L]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
% Pressures
subplot(3,2,5);lc='#0072bD';rc='#77AC30';
plot(x, P_d,'Color', lc,'LineWidth',lw);xlabel('x [m]','Fontsize',10); ylabel('[Pa]','Fontsize',10); ay=gca; ay.YAxis.Exponent = 5; hold on
yyaxis right
plot(x, P_f,'Color', rc,'LineWidth',lw); ylabel('[Pa]','Fontsize',10); legend('P_d^{RO}(x)','P_f^{RO}(x)','Location','best');xlim([0 swro_L]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
% water flow
subplot(3,2,2);lc='#0072bD';rc='#77AC30';
plot(x, J_wd,'Color', lc,'LineWidth',lw);xlabel('x [m]','Fontsize',10); ylabel('[kg/sm]','Fontsize',10); ay=gca; ay.YAxis.Exponent = 0;xlim([0 swro_L]); hold on
yyaxis right
plot(x, J_wf,'Color', rc,'LineWidth',lw); ylabel('[kg/sm]','Fontsize',10); legend('J_{w,d}^{RO}(x)','J_{w,f}^{RO}(x)','Location','East');xlim([0 swro_L]); ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
% salt flow
subplot(3,2,4);lc='#0072bD';rc='#77AC30';
plot(x, J_sd,'Color', lc,'LineWidth',lw);xlabel('x [m]','Fontsize',10); ylabel('[kg/sm]','Fontsize',10); ay=gca; hold on
yyaxis right; C_f=C_f(2:end);
plot(x(2:end), J_sf(2:end),'Color', rc,'LineWidth',lw); ylabel('[kg/sm]','Fontsize',10); legend('J_{s,d}^{RO}(x)','J_{s,f}^{RO}(x)','Location','best');xlim([0 swro_L]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
% density
subplot(3,2,6);lc='#0072bD';rc='#77AC30';
plot(x, swro_local_ro_d,'Color', lc,'LineWidth',lw);xlabel('x [m]','Fontsize',10); ylabel('[kg/m^3]','Fontsize',10); ay=gca; hold on
yyaxis right
plot(x(2:end), swro_local_ro_f(2:end),'Color', rc,'LineWidth',lw); ylabel('[kg/m^3]','Fontsize',10);xlim([0 swro_L]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;legend('\rho_d^{RO}(x)','\rho_f^{RO}(x)','Location','best');
end
%% figure 2
if fig(2) == 1
f=figure(2); %
f.Position = [1.2977e+03 635.6667 1.0987e+03 639.3333];
% Permeate flows 
subplot(2,2,1); lc='k';rc='#b81414'; J_cross2=J_cross(2:end);J_sin2=J_sin(2:end);
plot(x2, J_cross2, 'Color', lc, 'LineWidth',lw);xlabel('x [m]','Fontsize',10);  ylabel('[kg/sm^2]','Fontsize',10);ay=gca;hold on
yyaxis right
plot(x2, J_sin2, 'Color', rc, 'LineWidth',lw); ylabel('[kg/sm^2]','Fontsize',10); legend('J_{w,in}^{RO}(x)','J_{s,in}^{RO}(x)','Location','best');xlim([0 swro_L]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
% Osmotic/ Hydraulic pressure difference
subplot(2,2,3);lc='k';rc='#b81414'; osm_diff=swro_p_osm_d(2:end)-swro_p_osm_f(2:end); 
plot(x, (P_d-P_f),'Color', lc,'LineWidth',lw);xlabel('x [m]','Fontsize',10); ylabel('[Pa]','Fontsize',10); ay=gca;hold on
yyaxis right; 
plot(x2, osm_diff,'Color', rc,'LineWidth',lw); ylabel('[Pa]','Fontsize',10); legend('\Delta P^{RO}(x)','\Delta \pi^{RO}(x)','Location','best');xlim([0 swro_L]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
figure(1);
end
end
%x1=swro_Z*J_d(1)/swro_local_ro_d(1);
%x2=swro_Z*J_f(end)/swro_local_ro_f(end);
%x3=1000*C_f(end)*ro_salt;
%display(['SWRO feed flow rate:     ',num2str(x1),'m^3/s'])
%display(['Permeate  flow rate:     ',num2str(x2),'m^3/s'])
%
%display(' ')
%display(['Permeate concentration:      ',num2str(x3),' g/m^3'])
end
