function [output1, output2] = fun_unscaled(input,option_data,obj,option_mesh,option_BVP)
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

%% Read data
%
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
x = linspace(0,swro_L,n);    
  
%% Create initial guess for the BVP
J_wd_0 = 37.36;          % guess for J_wd(0)
J_wf_0 = 35;             % should be close to 0
Pf_0   = 1.1e5*Pf_L;                % guess for P_f(0)

y_init = [cE; J_wd_0; J_sf_0; J_wf_0; Pd_0; Pf_0];
solinit = bvpinit(x,y_init);

%% Solve the BVP
try
if version(6) == 0; sol = bvp5c(@(x,J_p)Unscaled_ODEsystem(x, J_p, DATA), @(ya,yb)Unscaled_Boundary1(ya,yb, DATA),solinit,ode_options); end
if version(6) == 1; sol = bvp5c(@(x,J_p)Unscaled_ODEsystem(x, J_p, DATA), @(ya,yb)Unscaled_Boundary2(ya,yb, DATA),solinit,ode_options); end
if version(6) > 1; fprintf(2,' \nERROR: Unscaled model only supports SWRO with/without ERD! \n'); end

%% Evaluate the solution and calculate new mixing ratio
if version(6)>0
    y = deval(sol,x); Y = y'; Y = real(Y);
    mix_M1 = ((Y(end,1).*Y(end,2)+Y(end,2))*(1-eta_ERD))/(Y(1,1).*Y(1,2)+Y(1,2));  
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
    if option_data == 0; DATA = @(x)Senthil_data(input,mix_M1);end
    if option_data == 1; DATA = @(x)Case_1_data(input,mix_M1); end
    if option_data == 2; DATA = @(x)Case_2_data(input,mix_M1); end
    if option_data == 3; DATA = @(x)Case_3_data(input,mix_M1); end
    % solve with bvp
    if version(6) == 0; sol = bvp5c(@(x,J_p)Unscaled_ODEsystem(x, J_p, DATA), @(ya,yb)Unscaled_Boundary1(ya,yb, DATA),solinit,ode_options); end
    if version(6) == 1; sol = bvp5c(@(x,J_p)Unscaled_ODEsystem(x, J_p, DATA), @(ya,yb)Unscaled_Boundary2(ya,yb, DATA),solinit,ode_options); end
end

% Evaluate the solution
if any([0,0,0,0] ~= fig)
    n=100; % higher resolution plots
end
x=linspace(0,swro_L,n);
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
swro_p_osm_d = ro_water*Rw*T0*log(1 + 2*Mw*Y(:,1)/Ms);
swro_p_osm_f = ro_water*Rw*T0*log(1 + 2*Mw*Y(:,3)/Ms./Y(:,4));
% local density
swro_local_ro_d = (Y(:,1) + 1)./(Y(:,1)./ro_salt + ones(n,1)./ro_water);
swro_local_ro_f = (Y(:,3) + Y(:,4))./(Y(:,3)./ro_salt + Y(:,4)./ro_water);
    % Salt permeability
    if version(4) == 0
        swro_beta = zeros(n,1);
    else
        swro_beta = swro_beta_fix.*ones(n,1);
    end
    % Water Permeate flux J_win(x)
    if version(4) == 0 % ideal
        J_cross = swro_alpha.*((Y(:,5) - Y(:,6)) - sigma*(swro_p_osm_d - swro_p_osm_f));
    elseif version(7) == 0 % ICP
        J_cross = swro_alpha .* ((Y(:,5) - Y(:,6)) - sigma.*(swro_p_osm_d - swro_p_osm_f))./(1 + swro_alpha*swro_KK*sigma.*(swro_p_osm_d - swro_p_osm_f));
    else % ICP+ECP
        J_cross = swro_alpha .* ((Y(:,5) - Y(:,6)) .* (1 + swro_beta .* (-1 ./ swro_KD + swro_KK + 1 ./ swro_KF)) - sigma .* (swro_p_osm_d - swro_p_osm_f)) ./ (1 + swro_beta .* (-1 ./ swro_KD + swro_KK + 1 ./ swro_KF) - swro_alpha .* sigma .* (-swro_p_osm_d ./ swro_KD + swro_p_osm_f .* swro_KK + swro_p_osm_f ./ swro_KF));
    end
    % Salt Permeate J_sin(x)
    if version(7) == 1 % ICP+ECP
        J_sin = swro_beta .* ((C_d-C_f) + C_d .* J_cross ./ swro_KD - C_f .* J_cross .* (swro_KK + 1 ./ swro_KF)) ./ (1 + swro_beta .* (-1 ./ swro_KD + swro_KK + 1 ./ swro_KF));
    else
        J_sin = swro_beta_fix.*(C_d-C_f);
    end
   

%% Output of System
% seawater enters the system
J_E = J_d(1);
J_wE = J_E/(cE+1);
J_sE = J_E-J_wE;

%% Version(6)=0 -->  only SWRO (no ERD)
if version(6)==0
    W_p1 = 1/HP_eff * (1e5 - P_d(1))*(J_E *swro_Z)/rho_E;
    W_p2 = 0; W_p3=0; W_p4=0; W_t=0;
end 

%% Version(6)=1 -->  only SWRO (with ERD)
if version(6)==1
    % flow from M1 to ERD1
    J_M1 = J_E * (1*mix_M1);
    J_w_M1 = J_wE* (1*mix_M1);
    J_s_M1 = J_sE* (1*mix_M1);
    % flow from ERD1 to RO
    rho_ERD= V_m*(swro_local_ro_d(end) - rho_E)+rho_E;
    C_ERD= -ro_salt*(rho_ERD-ro_water)/(ro_water*(rho_ERD-ro_salt));
    J_ERD = J_M1;
    J_wERD= J_ERD/(C_ERD+1);
    J_sERD= J_ERD - J_wERD;
    % flow from ERD1 to EXIT
    J_exit = J_d(end)*(1-eta_ERD);
    c_exit= (J_s_M1 + J_sd(end)*(1-eta_ERD) - J_sERD)/(J_w_M1 + J_wd(end)*(1-eta_ERD) - J_wERD);
    J_w_exit= J_exit/(c_exit+1);
    rho_exit = (J_exit)/(((c_exit*J_w_exit)/ro_salt)+(J_w_exit/ro_water)); 
    % energy balance
    f_1 = ERD_fric * rho_r*mix_density * (J_exit*swro_Z/rho_exit)*((J_M1*swro_Z)/(rho_E*A_ERD))^2; 
    f_2 = ERD_fric * rho_r*mix_density * (J_ERD*swro_Z/rho_ERD)*((J_d(end)*swro_Z)/(swro_local_ro_d(end)*A_ERD))^2;
    pERD = rho_ERD*(ERD_eff*(P_d(end)*J_d(end)*(1-eta_ERD)*swro_Z/swro_local_ro_d(end)-1e5*J_exit*swro_Z/rho_exit - f_2) + 1e5*J_M1*swro_Z/rho_E + f_1)/J_ERD/swro_Z;
    % pumps
    W_p1 = 1/HP_eff * (1e5 - P_d(1))*((J_E-J_M1)*swro_Z)/rho_E;
    W_p3 = 1/HP_eff * (pERD - P_d(1))*(J_ERD*swro_Z)/rho_ERD; 
    W_p2 = 0; W_p4=0; W_t=0;
end

%% Version(6)=2 or higher
if version(6) > 1; fprintf(2,' \nERROR: Unscaled model only supports SWRO with/without ERD! \n'); end

%% final output:
W_net= W_p1 + W_p2 + W_p3 + W_p4 + W_t; % in [W]

%% Calculate the final output
SEC_net = W_net*(swro_local_ro_f(end))./(J_f(end)*swro_Z)/1000/3600; % in [kWh/m^3]  
    if SEC_net > 0 && contains('solfig',obj)==1
        fprintf(2,' \nWaring: SEC_net is positive! \n');
    end
FW = (J_wf(end)*swro_Z/(swro_local_ro_f(end)))*3600; % in [m^3/h] 
    if FW < 0 && contains('solfig',obj)==1
        fprintf(2,' \nWaring: Freshwater production is negative! \n');
    end
Rev= pw*FW + pe*SEC_net*FW; %in [$/h]
SWRO_Recovery =(J_wf(end)./J_d(1))*100;     % in [%]
    if SWRO_Recovery > 100 && contains('solfig',obj)==1
        fprintf(2,' \nWaring: SWRO recovery is greater than 100 %% \n');
    end
    if SWRO_Recovery < 0 && contains('solfig',obj)==1
        fprintf(2,' \nWaring: SWRO recovery is negative! \n');
    end
PRO_Recovery = NaN; % in [%]
RO_inflow = swro_Z*J_d(1)/swro_local_ro_d(1); % in [m^3/s]
Permeate_outflow = swro_Z*J_f(end)/swro_local_ro_f(end); % in [m^3/s]
Wastewater_inflow = NaN; % in [m^3/s]
C_permeate = 10000*C_f(end); % in [ppm]
C_brine = 100*C_d(end); % in [%]
C_dilluted = P_d(end)/1e5; % in [%]
if version(6)==0; mix_M1=NaN; end

output1 = [SEC_net, FW, Rev, SWRO_Recovery, PRO_Recovery, ...
           RO_inflow, Permeate_outflow, Wastewater_inflow, C_permeate, C_brine, C_dilluted, mix_M1...
           W_net, W_p1, W_p2,  W_p3, W_p4, W_t]; % length(output1) = 5+7+6 = 18

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
        case 'Pareto' % SEC_net and FW maximization  -->  output1 = [-SEC_net, -FW]
            if sol.stats.maxerr > option_BVP    
            output1 = NaN(1,2); 
            else
                output1 = [-SEC_net, -FW];
                if FW < 0.15
                    output1 = NaN(1,2); 
                end
            end 
        case 'sol' 
            if sol.stats.maxerr > option_BVP
            output1 = NaN(1,18);
            end
    case 'fig' % case ends after the figures
            if sol.stats.maxerr > option_BVP
            output1 = NaN(1,18);
            else
%% figure 1
close all
lw=1.5; %Linewidth for all figures
if fig(1) == 1
    % initialize
    f=figure(5); 
    f.Position = [1839 227 1100 1000];
    x2=x(2:end);
    % total mass flows
    subplot(3,2,1);lc='#0072bD';rc='#77AC30';
    plot(x, J_d,'Color', lc,'LineWidth',lw); xlabel('x [m]','Fontsize',10); ylabel('[kg/sm]','Fontsize',10); ay=gca; ay.YAxis.Exponent = 0;xlim([0 swro_L]); hold on
    yyaxis right
    plot(x, J_f,'Color', rc,'LineWidth',lw); ylabel('[kg/sm]','Fontsize',10); legend('J_{d}^{RO}(x)','J_{f}^{RO}(x)','Location','East');xlim([0 swro_L]); ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
    % concentrations
    subplot(3,2,3);lc='#0072bD';rc='#77AC30';
    plot(x, 100*C_d*C_r,'Color', lc,'LineWidth',lw); xlabel('x [m]','Fontsize',10); ylabel('[%]','Fontsize',10); ay=gca; ay.YAxis.Exponent = 0;hold on
    yyaxis right; C_f=C_f(2:end);
    %plot(x2, 100*C_f*C_r,'Color', rc,'LineWidth',lw); ylabel('[%]','Fontsize',10); legend('C_d^{RO}(x)','C_f^{RO}(x)','Location','best');xlim([0 swro_L]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
    plot(x2, 1e4*C_f*C_r,'Color', rc,'LineWidth',lw); ylabel('[ppm]','Fontsize',10); legend('C_d^{RO}(x)','C_f^{RO}(x)','Location','best');xlim([0 swro_L]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
    % Pressures
    subplot(3,2,5);lc='#0072bD';rc='#77AC30';
    plot(x, P_d,'Color', lc,'LineWidth',lw); xlabel('x [m]','Fontsize',10); ylabel('[Pa]','Fontsize',10); ay=gca; ay.YAxis.Exponent = 5; hold on
    yyaxis right
    plot(x, P_f,'Color', rc,'LineWidth',lw); ylabel('[Pa]','Fontsize',10); legend('P_d^{RO}(x)','P_f^{RO}(x)','Location','best');xlim([0 swro_L]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
    % water flow
    subplot(3,2,2);lc='#0072bD';rc='#77AC30';
    plot(x, J_wd,'Color', lc,'LineWidth',lw); xlabel('x [m]','Fontsize',10); ylabel('[kg/sm]','Fontsize',10); ay=gca; ay.YAxis.Exponent = 0;xlim([0 swro_L]); hold on
    yyaxis right
    plot(x, J_wf,'Color', rc,'LineWidth',lw); ylabel('[kg/sm]','Fontsize',10); legend('J_{w,d}^{RO}(x)','J_{w,f}^{RO}(x)','Location','East');xlim([0 swro_L]); ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
    % salt flow
    subplot(3,2,4);lc='#0072bD';rc='#77AC30';
    plot(x, J_sd,'Color', lc,'LineWidth',lw); xlabel('x [m]','Fontsize',10); ylabel('[kg/sm]','Fontsize',10); ay=gca; hold on
    yyaxis right; C_f=C_f(2:end);
    plot(x(2:end), J_sf(2:end),'Color', rc,'LineWidth',lw); ylabel('[kg/sm]','Fontsize',10); legend('J_{s,d}^{RO}(x)','J_{s,f}^{RO}(x)','Location','best');xlim([0 swro_L]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
    % density
    subplot(3,2,6);lc='#0072bD';rc='#77AC30';
    plot(x, swro_local_ro_d,'Color', lc,'LineWidth',lw); xlabel('x [m]','Fontsize',10); ylabel('[kg/m^3]','Fontsize',10); ay=gca; hold on
    yyaxis right
    plot(x(2:end), swro_local_ro_f(2:end),'Color', rc,'LineWidth',lw); ylabel('[kg/m^3]','Fontsize',10);xlim([0 swro_L]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;legend('\rho_d^{RO}(x)','\rho_f^{RO}(x)','Location','best');
end
%% figure 2
if fig(2) == 1
    f=figure(6); %
    f.Position = [1.2977e+03 635.6667 1.0987e+03 639.3333];
    % Permeate flows 
    subplot(2,2,1); lc='k';rc='#b81414'; J_cross2=J_cross(2:end);J_sin2=J_sin(2:end);
    plot(x2, J_cross2, 'Color', lc, 'LineWidth',lw);xlabel('x [m]','Fontsize',10);  ylabel('[kg/sm^2]','Fontsize',10);ay=gca;hold on
    yyaxis right
    plot(x2, J_sin2, 'Color', rc, 'LineWidth',lw); ylabel('[kg/sm^2]','Fontsize',10); legend('J_{w,in}^{RO}(x)','J_{s,in}^{RO}(x)','Location','best');xlim([0 swro_L]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
    % Osmotic/ Hydraulic pressure difference
    subplot(2,2,3);lc='k';rc='#b81414'; osm_diff=swro_p_osm_d(2:end)-swro_p_osm_f(2:end); 
    a=min([(P_d-P_f); osm_diff]);b=max([(P_d-P_f); osm_diff]);
    plot(x, (P_d-P_f),'Color', lc,'LineWidth',lw);xlabel('x [m]','Fontsize',10); ylabel('[Pa]','Fontsize',10); ay=gca;ylim([0.95*a 1.05*b]);hold on
    yyaxis right; 
    plot(x2, osm_diff,'Color', rc,'LineWidth',lw); ylabel('[Pa]','Fontsize',10); legend('\Delta P^{RO}(x)','\Delta \pi^{RO}(x)','Location','best');xlim([0 swro_L]);ylim([0.95*a 1.05*b]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
end
%% figure 3
if fig(3) == 1 && version(6) > 0
f=figure(3);  % ERD quantities
f.Position = [637.6667 210.3333 900 347.3334]; tiledlayout(1,4); 
rot=0; rx=0; % rotation and shift of texts
vals1 = [1e5; pERD; 0; P_d(end); 1e5];
vals2 = 100*[cE; C_ERD; 0; C_d(end); c_exit];
vals3 = [J_E*mix_M1; J_ERD; 0; J_d(end); J_exit];
nexttile
b = bar(1,vals1, 'FaceColor', 'b'); b(4).FaceColor = "#77AC30"; b(5).FaceColor = [.2 .6 .5]; b(1).FaceColor = '#0072bD';
xticklabels({'[Pa]'}); ylim([0 1.1*max(vals1)]);ay=gca;ay.YAxis(1).Exponent=5;
t=text(b(1).XEndPoints,b(1).YEndPoints,"P_{s}^{1;in} ",'HorizontalAlignment','center','VerticalAlignment','bottom');
t=text(b(2).XEndPoints+rx ,b(2).YEndPoints,"P_{s}^{1;out}",'HorizontalAlignment','center','VerticalAlignment','bottom'); 
t=text(b(4).XEndPoints,b(4).YEndPoints,"P_{b}^{1;in}",'HorizontalAlignment','center','VerticalAlignment','bottom'); 
t=text(b(5).XEndPoints,b(5).YEndPoints,"   P_{b}^{1;out}",'HorizontalAlignment','center','VerticalAlignment','bottom'); 
nexttile
b = bar(1,vals2, 'FaceColor', 'b'); b(4).FaceColor = "#77AC30"; b(5).FaceColor = [.2 .6 .5]; b(1).FaceColor = '#0072bD';
xticklabels({'[%]'}); ylim([0 1.12*max(vals2)]);ay=gca;ay.YAxis(1).Exponent=0;
t= text(b(1).XEndPoints,b(1).YEndPoints,"C_{s}^{1;in}  ",'HorizontalAlignment','center','VerticalAlignment','bottom'); t.Rotation = rot;
t=text(b(2).XEndPoints+rx ,b(2).YEndPoints,"   C_{s}^{1;out}",'HorizontalAlignment','center','VerticalAlignment','bottom'); t.Rotation = rot;
t=text(b(4).XEndPoints,b(4).YEndPoints,"C_b^{1;in}  ",'HorizontalAlignment','center','VerticalAlignment','bottom'); t.Rotation = rot;
t=text(b(5).XEndPoints,b(5).YEndPoints,"  C_{b}^{1;out}",'HorizontalAlignment','center','VerticalAlignment','bottom'); t.Rotation = rot;
nexttile
b = bar(1,vals3, 'FaceColor', 'b'); b(4).FaceColor = "#77AC30"; b(5).FaceColor = [.2 .6 .5]; b(1).FaceColor = '#0072bD';
xticklabels({'[kg/ms]'}); ylim([0 1.1*max(vals3)]);ay=gca;ay.YAxis(1).Exponent=0;
t=text(b(1).XEndPoints,b(1).YEndPoints,"J_{s}^{1;in} ",'HorizontalAlignment','center','VerticalAlignment','bottom'); t.Rotation = rot;
t=text(b(2).XEndPoints+rx ,b(2).YEndPoints," J_{s}^{1;out}",'HorizontalAlignment','center','VerticalAlignment','bottom'); t.Rotation = rot;
t=text(b(4).XEndPoints ,b(4).YEndPoints,"J_{b}^{1;in} ",'HorizontalAlignment','center','VerticalAlignment','bottom'); t.Rotation = rot;
t=text(b(5).XEndPoints,b(5).YEndPoints," J_{b}^{1;out}",'HorizontalAlignment','center','VerticalAlignment','bottom'); t.Rotation = rot ;
nexttile
lgd = legend([b(1) b(2) b(5) b(4)],'Seawater inlet','Seawater outlet','Brine inlet','Brine outlet' , 'Location', 'EastOutside');
lgd.Layout.Tile = 4; axis off ; 
end
end
end
end