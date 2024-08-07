function [output1, output2, output3] = fun(input,option_data,obj,option_mesh,option_BVP)
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
%       obj           -   'SEC', 'FW', 'Rev', 'PD_net' , 'SE_net', 'approx'  or 'sol' ;
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
%       [Sim_01_output1, Sim_01_output2]=fun(x0,D,obj,mesh,tol),
%
tic;
%% Read data
if option_data == 0; DATA = @(x)Sim_0_data(input); end
% Simulations
if option_data == 0.11; DATA = @(x)Sim_11_data(input); end
if option_data == 0.12; DATA = @(x)Sim_12_data(input); end
if option_data == 0.13; DATA = @(x)Sim_13_data(input); end
if option_data == 0.21; DATA = @(x)Sim_21_data(input); end
if option_data == 0.22; DATA = @(x)Sim_22_data(input); end
if option_data == 0.23; DATA = @(x)Sim_23_data(input); end
if option_data == 0.24; DATA = @(x)Sim_24_data(input); end
if option_data == 0.31; DATA = @(x)Sim_31_data(input); end
if option_data == 0.32; DATA = @(x)Sim_32_data(input); end
% PRO Membrane Comparison
if option_data == 0.51; DATA = @(x)Sim_PRO_01_data(input); end
if option_data == 0.52; DATA = @(x)Sim_PRO_02_data(input); end
if option_data == 0.53; DATA = @(x)Sim_PRO_03_data(input); end
% Rough vs. Detailled
if option_data == 0.91; DATA = @(x)Sim_rough_approx_data(input); end
% Optimizations
if option_data==1; DATA = @(x)Opt_1_data(input); end
if option_data==2; DATA = @(x)Opt_2_data(input); end
if option_data==3; DATA = @(x)Opt_3_data(input); end
if option_data==4; DATA = @(x)Opt_4_data(input); end
if option_data==5; DATA = @(x)Opt_L_data(input); end
% for 3d plot
if option_data==11; DATA = @(x)plot3d_data(input); end

[H,Z,swro_Z,ro_water,ro_salt,Mw,Ms,Rw,T0,eta,sigma,p_r,rho_r,C_r,swro_L,swro_alpha,swro_R,swro_KK,swro_x_r,swro_b1,swro_b2,J_r,swro_gamma,swro_gamma2,swro_W_r,L,alpha,R,KK,x_r,b1,b2,Q_r,gamma,gamma2,W_r,cE,pE,rho_E,J_sf_0,J_wf_0,Pd_0,Pd_L,Pf_L,Q_sf_0,pd_0,pf_0,pd_L,pf_L,HP_eff,LP_eff,T_eff,V_m,ERD_eff,ERD_fric,A_ERD,eta_ERD,mix_density,pw,pe,swro_beta_fix,beta_fix,mixer_ERD,version,fig,swro_KF,swro_KD, KF, KD] ...
    = DATA(1);

%% Set options
ode_options = bvpset('Stats','off','NMax', option_mesh, 'RelTol', option_BVP);      

try
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

%% Solve the BVP
if version(6) == 0; sol = bvp5c(@(x,J_p)ERD_model(x, J_p, DATA), @(ya,yb)ERD_boundary2(ya,yb, DATA),solinit,ode_options); end
if version(6) == 1; sol = bvp5c(@(x,J_p)ERD_model(x, J_p, DATA), @(ya,yb)ERD_boundary1(ya,yb, DATA),solinit,ode_options); end
if version(6) == 2; sol = bvp5c(@(x,J_p)ERD_model(x, J_p, DATA), @(ya,yb)ERD_boundary2(ya,yb, DATA),solinit,ode_options); end
if version(6) == 3; sol = bvp5c(@(x,J_p)ERD_model(x, J_p, DATA), @(ya,yb)ERD_boundary3(ya,yb, DATA),solinit,ode_options); end       
if version(6) == 4; sol = bvp5c(@(x,J_p)ERD_model(x, J_p, DATA), @(ya,yb)ERD_boundary4(ya,yb, DATA),solinit,ode_options); end   
            
%% Evaluate the solution
n=100;
x=linspace(0,1,n);
y = deval(sol,x); Y = y'; Y = real(Y);
C_d=zeros(1,n);C_f=zeros(1,n);J_sd=zeros(1,n);J_wd=zeros(1,n);J_d=zeros(1,n);J_sf=zeros(1,n);J_wf=zeros(1,n);J_f=zeros(1,n);P_d=zeros(1,n);P_f=zeros(1,n);swro_p_osm_d=zeros(1,n);swro_p_osm_f=zeros(1,n);swro_beta=zeros(1,n);Y(:,4+6)=abs(Y(:,4+6));swro_del_c=zeros(1,n);J_cross=zeros(1,n);J_sin=zeros(1,n);
c_d=zeros(1,n);c_f=zeros(1,n);Q_sd=zeros(1,n);Q_wd=zeros(1,n);Q_d=zeros(1,n);Q_sf=zeros(1,n);Q_wf=zeros(1,n);Q_f=zeros(1,n);p_d=zeros(1,n);p_f=zeros(1,n);p_osm_d=zeros(1,n);p_osm_f=zeros(1,n);beta=zeros(1,n);del_c=zeros(1,n);Q_cross=zeros(1,n);Q_sin=zeros(1,n); local_ro_d=zeros(1,n);local_ro_f=zeros(1,n);swro_local_ro_d=zeros(1,n);swro_local_ro_f=zeros(1,n);

    for i=1:n
    %% SWRO evaluation
    % Concentrations
    C_d(i)= Y(i,1);                 % concentration of draw side
    C_f(i)= Y(i,3)./Y(i,4);         % concentration of fresh side
    % Flow rates in draw side
    J_sd(i)= Y(i,1)*Y(i,2);         % salt flow rate in draw side
    J_wd(i)= Y(i,2);                % water flow rate in draw side
    J_d(i) = Y(i,1)*Y(i,2)+Y(i,2);  % total flow rate draw side
    % Flow rates in fresh side
    J_sf(i)= Y(i,3);                % salt flow rate in fresh side
    J_wf(i)= Y(i,4);                % water flow rate in fresh side
    J_f(i) = Y(i,3)+Y(i,4);         % total flow rate fresh side
    % Pressures 
    P_d(i)=Y(i,5);                  % Pressure draw side
    P_f(i)=Y(i,6);                  % Pressure fresh side

    % Osmotic pressure
    swro_p_osm_d(i) = ro_water*Rw*T0*log(1 + 2*Mw*Y(i,1)/Ms)/p_r;
    swro_p_osm_f(i) = ro_water*Rw*T0*log(1 + 2*Mw*Y(i,3)/Ms./Y(i,4))/p_r;
    % Concentration difference
    swro_del_c(i) = Y(i,1)./(Y(i,1) + 1) - Y(i,3)./(Y(i,3) + Y(i,4));
    % Local density
    swro_local_ro_d(i) = (Y(i,1) + 1)./(ro_water*Y(i,1)/ro_salt + 1);           % local density of draw side
    swro_local_ro_f(i) = (Y(i,3) + Y(i,4))./(ro_water*Y(i,3)/ro_salt + Y(i,4)); % local density of fresh side
    % B(x)
    swro_beta(i) = (1 - swro_R)*((Y(i,5) - Y(i,6)) - sigma.*(swro_p_osm_d(i) - swro_p_osm_f(i)))/swro_R;
    if version(2) == 0 ;swro_beta(i)= swro_beta_fix; end
    if version(4) == 0 ;swro_beta(i)= 0; end
    % J_w,in and J_s,in
    
    J_cross(i) = ((Y(i,5) - Y(i,6)) - sigma.*(swro_p_osm_d(i) - swro_p_osm_f(i)))./(1 + p_r*swro_alpha*swro_KK*sigma.*(swro_p_osm_d(i) - swro_p_osm_f(i)));
    if version(4) == 0; J_cross(i) = ((Y(i,5) - Y(i,6)) - sigma.*(swro_p_osm_d(i) - swro_p_osm_f(i)));  end 
    J_sin(i)=swro_beta(i)*swro_del_c(i);
 
    % J_w,in and J_s,in with ECP and ICP
    if version(7)==1
    J_cross(i) = ((Y(i,5)-Y(i,6))-sigma.*(swro_p_osm_d(i) - swro_p_osm_f(i)) + swro_beta(i)*(Y(i,5)-Y(i,6))*(swro_KK +1/swro_KD + 1/swro_KF)) ./   (1 + swro_alpha*swro_beta(i)*(1/swro_KD + 1/swro_KF + swro_KK) + p_r*swro_alpha*sigma.*(-swro_p_osm_d(i)/swro_KD -swro_p_osm_f(i)*(swro_KK + 1/swro_KF) ));
    J_sin(i) = swro_beta(i)*( swro_del_c(i)-J_cross(i)*( Y(i,1)./(Y(i,1) + 1)./swro_KD +  Y(i,3)./(Y(i,3) + Y(i,4))*(swro_KK+1/swro_KF)))/(1+swro_alpha*swro_beta(i)*(swro_KK+1/swro_KF+1/swro_KD));
    end
    
%% PRO evaluation
    % Concentrations
    c_d(i)= Y(i,1+6);                    % concentration of draw side
    c_f(i)= Y(i,3+6)./Y(i,4+6);          % concentration of fresh side
    % Flow rates in draw side
    Q_sd(i)= Y(i,1+6)*Y(i,2+6);          % salt flow rate in draw side
    Q_wd(i)= Y(i,2+6);                   % water flow rate in draw side
    Q_d(i) = Y(i,1+6)*Y(i,2+6)+Y(i,2+6); % total flow rate draw side
    % Flow rates in fresh side 
    Q_sf(i)= Y(i,3+6);                   % salt flow rate in fresh side
    Q_wf(i)= Y(i,4+6);                   % water flow rate in fresh side
    Q_f(i) = Y(i,3)+Y(i,4+6);            % total flow rate fresh side
    % Pressures
    p_d(i)=Y(i,5+6);                     % Pressure draw side
    p_f(i)=Y(i,6+6);                     % Pressure fresh side

    % Osmotic pressure
    p_osm_d(i) = ro_water*Rw*T0*log(1 + 2*Mw*Y(i,1+6)/Ms)/p_r;
    p_osm_f(i) = ro_water*Rw*T0*log(1 + 2*Mw*Y(i,3+6)/Ms./Y(i,4+6))/p_r;
    % Concentration difference
    del_c(i) = Y(i,1+6)./(Y(i,1+6) + 1) - Y(i,3+6)./(Y(i,3+6) + Y(i,4+6));
    % local density
    local_ro_d(i) = (Y(i,1+6) + 1)./(ro_water*Y(i,1+6)/ro_salt + 1);            	
    local_ro_f(i) = (Y(i,3+6) + Y(i,4+6))./(ro_water*Y(i,3+6)/ro_salt + Y(i,4+6)); 
    % B(x) 
    beta(i) = (1 - R)*((p_osm_d(i) - p_osm_f(i)) - (Y(i,5+6) - Y(i,6+6)))/R;
    if version(3) == 0 ; beta(i) =  beta_fix; end
    if version(5) == 0 ; beta(i) = 0; end
    % Q_w,in and Q_s,in
    Q_cross(i) = (p_osm_d(i) - p_osm_f(i) - (Y(i,5+6) - Y(i,6+6)) - (Y(i,5+6) - Y(i,6+6))*beta(i)*KK*alpha*p_r)./(1 + KK*alpha*p_r*(beta(i) +p_osm_f(i)));
    if version(5) == 0; Q_cross(i) = ((p_osm_d(i) - p_osm_f(i)) - (Y(i,5+6) - Y(i,6+6))) ; end
    Q_sin(i)=beta(i)*del_c(i);
    
    if version(8)==1 %% version with ICP and ECP
    Q_cross(i)=((p_osm_d(i) - p_osm_f(i)) - (Y(i,5+6) - Y(i,6+6)) - (Y(i,5+6) - Y(i,6+6))*beta(i)*(KK+1/KF+1/KD)*p_r*alpha)./(1 + alpha*beta(i)*(KK+1/KF+1/KD) + alpha*p_r*(p_osm_f(i)*(KK+1/KF) + p_osm_d(i)/KD));
    Q_sin(i)=beta(i)*( del_c(i)-Q_cross(i)*( Y(i,1+6)./(Y(i,1+6) + 1)./KD +  Y(i,3+6)./(Y(i,3+6) + Y(i,4+6))*(KK+1/KF)))/(1+alpha*beta(i)*(KK+1/KF+1/KD));
    end
    
    end
    % In ideal cases take mean of some quantities (one could instead use ylim on y-axis) 
    if version(4)==0 ; J_f=J_wf; J_sf=zeros(1,n); J_sd=mean(J_sd)*ones(1,n); C_f=zeros(1,n); end
    if version(5)==0 ; Q_f=Q_wf; Q_sf=zeros(1,n); Q_sd=mean(Q_sd)*ones(1,n); c_f=zeros(1,n); end
    
%% Analysis of Recovery rates
Freshwater = J_wf(end)*J_r; %in [kg/s] 
FW = (Freshwater*swro_Z/(swro_local_ro_f(end)*rho_r))*3600; %in [m^3/h] 
SWRO_Recovery =(J_f(end)./J_d(1))*100;     % SWRO freshwater recovery [%]
PRO_Recovery = (1- Q_f(end)./Q_f(1))*100;  % PRO recovery rate [%]

%% Output of System
J_E = (C_d(1)*J_wd(1)+J_wd(1))/2; J_ERD = J_E; J_wE = J_E/(cE+1);
W_p1 = 1/HP_eff * (P_d(1)-pE)*(J_E *swro_Z)/rho_E; W_p1 = W_p1*swro_W_r;
W_p2 = 0; W_p2 = W_p2*swro_W_r; 
W_p4 = 1/LP_eff * (p_f(1)-pE)*(Q_f(1)*Z)./local_ro_f(1); W_p4 = W_p4*W_r;

    %% Version(6)=0 -->  only SWRO (no ERD)
    if version(6)==0; W_p3=W_p1; W_p4=0; W_t=0; end

    %% Version(6)=1 -->  only SWRO (with ERD)
    if version(6)==1
        rho_ERD = V_m*(swro_local_ro_d(end) - rho_E)+rho_E;
        C_ERD = -ro_salt*(rho_ERD-1)/(rho_ERD*ro_water-ro_salt);
        J_wERD = J_E/(C_ERD+1);
        Q_exit = J_d(end)*(1-eta_ERD);
        c_exit = (C_d(end)*J_wd(end)*(1 - eta_ERD)+cE*J_wE-C_ERD*J_wERD)/(J_wd(end)*(1 - eta_ERD)+J_wE-J_wERD);
        rho_exit= (ro_salt*(c_exit+1))/(c_exit*ro_water+ro_salt); p_exit=pE;
        f_1= (Q_exit*swro_Z/rho_exit*mix_density*(J_E*swro_Z/rho_E/A_ERD)^2 );
        f_2= (J_ERD*swro_Z/rho_ERD*mix_density*(J_d(end)*swro_Z*(1-eta_ERD)/swro_local_ro_d(end)/A_ERD)^2); 
        f_ERD= ERD_fric*(ERD_eff*f_2 - f_1);
        pERD = rho_ERD*(ERD_eff*( P_d(end)*J_d(end)*swro_Z*(1-eta_ERD)/swro_local_ro_d(end) - p_exit*Q_exit*swro_Z/rho_exit ) +  pE*J_E*swro_Z/rho_E -f_ERD )/J_ERD/swro_Z;
        f_ERD=0; pMax=rho_ERD*(ERD_eff*( P_d(end)*J_d(end)*swro_Z*(1-eta_ERD)/swro_local_ro_d(end) - p_exit*Q_exit*swro_Z/rho_exit ) +  pE*J_E*swro_Z/rho_E -f_ERD )/J_ERD/swro_Z;
        if pERD > pMax; pERD=pMax; end
        if pERD < 0; pERD=pMax; end  
        %display(['pERD:              ',num2str(pERD)]) % pERD use in the model
        %display(['pERD_theoretical:  ',num2str(pMax)]) % theoretical pERD with no leak, friction perm and 100% efficiency
        W_p3 = 1/HP_eff * (P_d(1)-pERD)*(J_ERD*swro_Z)/rho_ERD; W_p3=W_p3*swro_W_r; W_p4=0; W_t=0;
    end
    %% Version(6)=2 --> only PRO
    if version(6)==2
        if version(1) ==0; W_t = T_eff * (p_d(end)-pE)*(Q_d(end)*Z)/local_ro_d(end); end
        if version(1) ==1; W_t = T_eff * (p_d(1)-pE)*(abs(Q_d(1))*Z)/local_ro_d(1); end
        W_p1=0; W_p3=0; W_t=W_t*W_r;
        if version(1) ==0; W_p2= 1/LP_eff * (p_d(1)-pE)*(Q_d(1)*Z)./local_ro_d(1); end 
        if version(1) ==1; W_p2= 1/LP_eff * (p_d(end)-pE)*(abs(Q_d(end))*Z)./local_ro_d(end); end 
        W_p2 = W_p2*W_r;
    end
    
    %% Version(6)=3 --> SWRO-PRO hybrid system (with one ERD)
    if version(6)==3
        rho_ERD = V_m*(swro_local_ro_d(end) - rho_E)+rho_E;
        C_ERD = -ro_salt*(rho_ERD-1)/(rho_ERD*ro_water-ro_salt);
        J_wERD = J_E/(C_ERD+1); qq=Q_r/J_r;
            if version(1) == 0 %co-current
                f_1= qq*Q_d(1)*Z/local_ro_d(1)*mix_density*(J_E*swro_Z/rho_E/A_ERD)^2;  
                f_2= J_ERD*swro_Z/rho_ERD*mix_density*(J_d(end)*swro_Z*(1-eta_ERD)/swro_local_ro_d(end)/A_ERD)^2;  
                f_ERD= ERD_fric*(ERD_eff*f_2 - f_1);
                pERD = rho_ERD*(ERD_eff*( P_d(end)*J_d(end)*swro_Z*(1-eta_ERD)/swro_local_ro_d(end) - p_d(1)*qq*Q_d(1)*Z/local_ro_d(1) ) +  pE*J_E*swro_Z/rho_E -f_ERD )/J_ERD/swro_Z;
                f_ERD=0;pMax=rho_ERD*(ERD_eff*( P_d(end)*J_d(end)*swro_Z*(1-eta_ERD)/swro_local_ro_d(end) - p_d(1)*qq*Q_d(1)*Z/local_ro_d(1) ) +  pE*J_E*swro_Z/rho_E -f_ERD )/J_ERD/swro_Z;
                if pERD > pMax; pERD=pMax; end
                %if pERD < 0; pERD=pMax; end   
                W_p3 = 1/HP_eff * (P_d(1)-pERD)*(J_ERD*swro_Z)/rho_ERD;
                W_t = T_eff * (p_d(end)-pE)*(Q_d(end)*Z)/local_ro_d(end); 
            end
            if version (1) == 1 %counter-current   
                f_1= abs(qq*Q_d(end))*Z/local_ro_d(end)*mix_density*(J_E*swro_Z/rho_E/A_ERD)^2  ;
                f_2= J_ERD*swro_Z/rho_ERD*mix_density*(J_d(end)*swro_Z*(1-eta_ERD)/swro_local_ro_d(end)/A_ERD)^2  ;
                f_ERD= ERD_fric*(ERD_eff*f_2 - f_1);
                pERD = rho_ERD*(ERD_eff*( P_d(end)*J_d(end)*swro_Z*(1-eta_ERD)/swro_local_ro_d(end) - p_d(end)*abs(qq*Q_d(end))*Z/local_ro_d(end) ) +  pE*J_E*swro_Z/rho_E -f_ERD )/J_ERD/swro_Z;     
                f_ERD=0;pMax=rho_ERD*(ERD_eff*( P_d(end)*J_d(end)*swro_Z*(1-eta_ERD)/swro_local_ro_d(end) - p_d(end)*abs(qq*Q_d(end))*Z/local_ro_d(end) ) +  pE*J_E*swro_Z/rho_E -f_ERD )/J_ERD/swro_Z;       
                %display(['pERD:              ',num2str(pERD)]) % pERD use in the model
                %display(['pERD_theoretical:  ',num2str(pMax)]) % theoretical pERD with no leak, friction perm and 100% efficiency
                if pERD > pMax; pERD=pMax; end
                %if pERD < 0; pERD=pMax; end     
                W_p3 = 1/HP_eff * (P_d(1)-pERD)*(J_ERD*swro_Z)/rho_ERD;
                W_t = T_eff * (p_d(1)-pE)*(abs(Q_d(1))*Z)/local_ro_d(1);
            end
        W_p3=W_p3*swro_W_r; W_t=W_t*W_r;
    end

    %% Version(6)=4 --> SWRO-PRO hybrid system (with two ERDs)
    if version(6)==4 
        if version(1) == 0; rho_ERD2= V_m*(local_ro_d(end) - rho_E)+rho_E; end
        if version(1) == 1; rho_ERD2= V_m*(local_ro_d(1) - rho_E)+rho_E; end
        C_ERD2= -ro_salt*(rho_ERD2-1)/(rho_ERD2*ro_water-ro_salt);
        rho_ERD= V_m*(swro_local_ro_d(end) - rho_ERD2)+rho_ERD2;
        C_ERD = -ro_salt*(rho_ERD-1)/(rho_ERD*ro_water-ro_salt);
        J_E = (C_d(1)*J_wd(1)+J_wd(1))/2; J_ERD = J_E; 
        J_wE = J_E/(C_ERD2+1); J_wERD= J_ERD/(C_ERD+1);
        qq=Q_r/J_r;
        J_E_sea = 2*J_E;
        J_wE_sea = J_E_sea/(cE+1);
        if version(1) == 0
            Q_exit = mixer_ERD*Q_d(end)*(1-eta_ERD);
            c_exit = (c_d(end)*mixer_ERD*Q_wd(end)*(1 - eta_ERD)+cE*J_wE_sea-C_ERD2*2*J_E)/(mixer_ERD*Q_wd(end)*(1 - eta_ERD)+J_wE_sea-2*J_E);
        end
        if version(1) == 1
            Q_exit = mixer_ERD*abs(Q_d(1))*(1-eta_ERD);
            c_exit = (c_d(1)*mixer_ERD*abs(Q_wd(1))*(1 - eta_ERD)+cE*J_wE_sea-C_ERD2*2*J_E)/(mixer_ERD*abs(Q_wd(1))*(1 - eta_ERD)+J_wE_sea-2*J_E);
        end
        rho_exit= (ro_salt*(c_exit+1))/(c_exit*ro_water+ro_salt); p_exit=pE;
        %calculation of Pressure P_ERD2
            if version(1) == 0 %co-current
                f_1= qq*Q_exit*Z/rho_exit*mix_density*(J_E_sea*swro_Z/rho_E/A_ERD)^2;  
                f_2= J_E*swro_Z/rho_ERD2*mix_density*(qq*mixer_ERD*Q_d(end)*Z*(1-eta_ERD)/local_ro_d(end)/A_ERD)^2;  
                f_ERD= ERD_fric*(ERD_eff*f_2 - f_1);
                pERD2 = rho_ERD2*(ERD_eff*( p_d(end)*qq*mixer_ERD*Q_d(end)*Z*(1-eta_ERD)/local_ro_d(end) - p_exit*qq*Q_exit*Z/rho_exit ) +  pE*J_E_sea*swro_Z/rho_E -f_ERD )/2/J_E/swro_Z;
                f_ERD=0;pMax = rho_ERD2*(ERD_eff*( p_d(end)*qq*mixer_ERD*Q_d(end)*Z*(1-eta_ERD)/local_ro_d(end) - p_exit*qq*Q_exit*Z/rho_exit ) +  pE*J_E_sea*swro_Z/rho_E -f_ERD )/2/J_E/swro_Z;
                if pERD2 > pMax; pERD2=pMax; end
                %if pERD2 < 0; pERD2=pMax; end  
                
            end
            if version (1) == 1 %counter-current 
                f_1= qq*Q_exit*Z/rho_exit*mix_density*(J_E_sea*swro_Z/rho_E/A_ERD)^2;  
                f_2= J_E*swro_Z/rho_ERD2*mix_density*(qq*mixer_ERD*abs(Q_d(1))*Z*(1-eta_ERD)/local_ro_d(1)/A_ERD)^2;  
                f_ERD= ERD_fric*(ERD_eff*f_2 - f_1);
                pERD2 = rho_ERD2*(ERD_eff*( p_d(1)*qq*mixer_ERD*abs(Q_d(1))*Z*(1-eta_ERD)/local_ro_d(1) - p_exit*qq*Q_exit*Z/rho_exit ) +  pE*J_E_sea*swro_Z/rho_E -f_ERD )/2/J_E/swro_Z;            
                f_ERD=0;pMax = rho_ERD2*(ERD_eff*( p_d(1)*qq*mixer_ERD*abs(Q_d(1))*Z*(1-eta_ERD)/local_ro_d(1) - p_exit*qq*Q_exit*Z/rho_exit ) +  pE*J_E_sea*swro_Z/rho_E -f_ERD )/2/J_E/swro_Z;            
                if pERD2 > pMax; pERD2=pMax; end
                %if pERD2 < 0; pERD2=pMax; end  
            end     
            %calculation of Pressure P_ERD
            if version(1) == 0 %co-current
                f_1= qq*Q_d(1)*Z/local_ro_d(1)*mix_density*(J_E*swro_Z/rho_ERD2/A_ERD)^2;  
                f_2= J_ERD*swro_Z/rho_ERD*mix_density*(J_d(end)*swro_Z*(1-eta_ERD)/swro_local_ro_d(end)/A_ERD)^2;  
                f_ERD= ERD_fric*(ERD_eff*f_2 - f_1);
                pERD = rho_ERD*(ERD_eff*( P_d(end)*J_d(end)*swro_Z*(1-eta_ERD)/swro_local_ro_d(end) - p_d(1)*qq*Q_d(1)*Z/local_ro_d(1) ) +  pERD2*J_E*swro_Z/rho_ERD2 -f_ERD )/J_ERD/swro_Z;
                f_ERD=0; pMax = rho_ERD*(ERD_eff*( P_d(end)*J_d(end)*swro_Z*(1-eta_ERD)/swro_local_ro_d(end) - p_d(1)*qq*Q_d(1)*Z/local_ro_d(1) ) +  pERD2*J_E*swro_Z/rho_ERD2 -f_ERD )/J_ERD/swro_Z;
                if pERD > pMax; pERD=pMax; end
                %if pERD < 0; pERD=pMax; end  
                W_p3 = 1/HP_eff * (P_d(1)-pERD)*(J_ERD*swro_Z)/rho_ERD;
                W_t = T_eff * (p_d(end)-pE)*((1-mixer_ERD)*Q_d(end)*Z)/local_ro_d(end); 
            end
            if version (1) == 1 %counter-current   
                f_1= abs(qq*Q_d(end))*Z/local_ro_d(end)*mix_density*(J_E*swro_Z/rho_ERD2/A_ERD)^2  ;
                f_2= J_ERD*swro_Z/rho_ERD*mix_density*(J_d(end)*swro_Z*(1-eta_ERD)/swro_local_ro_d(end)/A_ERD)^2  ;
                f_ERD= ERD_fric*(ERD_eff*f_2 - f_1);
                pERD = rho_ERD*(ERD_eff*( P_d(end)*J_d(end)*swro_Z*(1-eta_ERD)/swro_local_ro_d(end) - p_d(end)*abs(qq*Q_d(end))*Z/local_ro_d(end) ) +  pERD2*J_E*swro_Z/rho_ERD2 -f_ERD )/J_ERD/swro_Z;            
                f_ERD=0; pMax=rho_ERD*(ERD_eff*( P_d(end)*J_d(end)*swro_Z*(1-eta_ERD)/swro_local_ro_d(end) - p_d(end)*abs(qq*Q_d(end))*Z/local_ro_d(end) ) +  pERD2*J_E*swro_Z/rho_ERD2 -f_ERD )/J_ERD/swro_Z;            
                if pERD > pMax; pERD=pMax; end 
                %if pERD < 0; pERD=pMax; end  
                W_p3 = 1/HP_eff * (P_d(1)-pERD)*(J_ERD*swro_Z)/rho_ERD;
                W_t = T_eff * (p_d(1)-pE)*(abs((1-mixer_ERD)*Q_d(1))*Z)/local_ro_d(1);
            end         
    % Recalculation for Pump 1
    W_p1 = 1/HP_eff * (P_d(1)-pERD2)*(J_E *swro_Z)/rho_ERD2; W_p1 = W_p1*swro_W_r;
    end
    
%% final output:
W_net= - W_p1 - W_p2 - W_p3 - W_p4 + W_t; % in [W]
PD_net = W_net./(Z*L); %in [W/m^2]    
%old SE_net
%if version(1)==0; SE_net = W_net./(Q_f(1)*Q_r*Z/(rho_r*local_ro_f(1)) + Q_d(1)*Q_r*Z/(rho_r*local_ro_d(1)))/1000/3600; end %in [kWh/m^3]
%if version(1)==1; SE_net = W_net./(Q_f(1)*Q_r*Z/(rho_r*local_ro_f(1)) + abs(Q_d(end))*Q_r*Z/(rho_r*local_ro_d(end)))/1000/3600; end %in [kWh/m^3]
%new SE_tot
if version(1)==0; SE_net = W_net./(Q_f(1)*Q_r*Z/(rho_r*local_ro_f(1)))/1000/3600; end %in [kWh/m^3]
if version(1)==1; SE_net = W_net./(Q_f(1)*Q_r*Z/(rho_r*local_ro_f(1)))/1000/3600; end %in [kWh/m^3]
SEC_net = W_net*(swro_local_ro_f(end)*rho_r)./(J_f(end)*J_r*swro_Z)/1000/3600; %in [kWh/m^3]  
Rev= pw*FW + pe*SEC_net*FW; %in [$/h]
tt=toc;
%% if BVP-solver failed
catch
    sol.stats.maxerr =1000;
    SEC_net=1000;FW=0;Rev=-20;SWRO_Recovery=NaN;PRO_Recovery=NaN;PD_net=NaN;SE_net=NaN;
    W_net=NaN;W_p1=NaN;W_p2=NaN;W_p3=NaN;W_p4=NaN;W_t=NaN;
end
%% Output of objective function
switch (obj)
        case 'SEC'
            if sol.stats.maxerr > option_BVP
            output1 = 20;
            else; output1 = -SEC_net;
            end 
        case 'FW' 
            if sol.stats.maxerr > option_BVP
            output1 = 20;
            else; output1 = -FW;
            end 
        case 'Rev' 
            if sol.stats.maxerr > option_BVP
            output1 = 0;
            else; output1 = -Rev;
            end
        case 'Pareto' 
            if sol.stats.maxerr > option_BVP    
            output1 = [NaN, NaN];
            else; output1 = [-SEC_net, -FW];
                %if SWRO_Recovery < 30;output1 = [NaN, NaN];end
                %if SWRO_Recovery > 60;output1 = [NaN, NaN];end
                %if PRO_Recovery < 40;output1 = [NaN, NaN];end
                %if PRO_Recovery > 95;output1 = [NaN, NaN];end
            end 
        case 'PD_net' 
            if sol.stats.maxerr > option_BVP
            output1 = 0;
            else; output1 = -PD_net;
            end
        case 'SE_net' 
            if sol.stats.maxerr > option_BVP
            output1 = 0;
            else; output1 = -SE_net;
            end
        case 'approx' 
            if version(6) == 0; output1=[SEC_net, FW, Rev, SWRO_Recovery, NaN, NaN, NaN]; end
            if version(6) == 1; output1=[SEC_net, FW, Rev, SWRO_Recovery, NaN, NaN, NaN]; end
            if version(6) == 2; output1=[NaN, NaN, NaN, NaN, PRO_Recovery, PD_net, SE_net]; end
            if version(6) == 3; output1=[SEC_net, FW, Rev, SWRO_Recovery, PRO_Recovery, NaN, NaN]; end
            if version(6) == 4; output1=[SEC_net, FW, Rev, SWRO_Recovery, PRO_Recovery, NaN, NaN]; end
            %sol.stats.nmeshpoints
            %sol.stats.maxerr
            %sol.stats.nODEevals
            %sol.stats.nBCevals
            output2=[sol.stats.nmeshpoints sol.stats.maxerr];
    case 'exe' 
            if sol.stats.maxerr > option_BVP
            output1 = [NaN, NaN, NaN, NaN, NaN, NaN, NaN];
            output2=[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN,NaN, NaN, NaN, NaN, NaN, NaN,NaN, NaN, NaN, NaN, NaN, NaN];
            else
            if version(6) == 0; output1=[SEC_net, FW, Rev, SWRO_Recovery, NaN, NaN, NaN]; end
            if version(6) == 1; output1=[SEC_net, FW, Rev, SWRO_Recovery, NaN, NaN, NaN]; end
            if version(6) == 2; output1=[NaN, NaN, NaN, NaN, PRO_Recovery, PD_net, SE_net]; end
            if version(6) == 3; output1=[SEC_net, FW, SWRO_Recovery, PRO_Recovery]; end
            if version(6) == 4; output1=[SEC_net, FW, SWRO_Recovery, PRO_Recovery]; end
            output2=[x*swro_L; 100*C_d; J_r*J_d; p_r*P_d; p_r*P_f; x*L; 100*c_d; Q_r*Q_d; Q_r*Q_f; p_r*p_d; p_r*p_f]; % 11
            output3=[x(2:end)*swro_L; 100*C_f(2:end); J_r*J_f(2:end); J_r/swro_x_r*J_cross(2:end); J_r/swro_x_r*J_sin(2:end); x(2:end)*L; 100*c_f(2:end); Q_r/x_r*Q_cross(2:end); Q_r/x_r*Q_sin(2:end)]; %9
            %
            output1=output1';output2=output2';output3=output3';
            [a,b]=size(output1);
            output1=[1:1:b; output1];
            [a,b]=size(output2);
            output2=[1:1:b; output2];
            [a,b]=size(output3);
            output3=[1:1:b; output3];
            end
    case 'sol' 
            if sol.stats.maxerr > option_BVP
            output1 = [NaN, NaN, NaN, NaN, NaN, NaN, NaN];
            output2=[NaN, NaN, NaN, NaN, NaN, NaN];
            output3=NaN;
            else
            if version(6) == 0; output1=[SEC_net, FW, Rev, SWRO_Recovery, NaN, NaN, NaN]; end
            if version(6) == 1; output1=[SEC_net, FW, Rev, SWRO_Recovery, NaN, NaN, NaN]; end
            if version(6) == 2; output1=[NaN, NaN, NaN, NaN, PRO_Recovery, PD_net, SE_net]; end
            if version(6) == 3; output1=[SEC_net, FW, Rev, SWRO_Recovery, PRO_Recovery, NaN, NaN]; end
            if version(6) == 4; output1=[SEC_net, FW, Rev, SWRO_Recovery, PRO_Recovery, NaN, NaN]; end
            output2=[W_net, -W_p1, -W_p2,  -W_p3, -W_p4, W_t];
            output3=tt;
            end
  
            
%% figure 1
if sol.stats.maxerr ==1000
fprintf(2,' \nERROR: bvp5c could not satisfy the relative error tolerance ---> no figures could be displayed\n');
else
lw=1.5; %Linewidth for all figures
if fig(1) == 1
f=figure(1); %
f.Position = [276.3333 196.3333 1100 1000];
x2=x(2:end);
% total mass flows
subplot(3,2,1);lc='#0072bD';rc='#77AC30';
plot(x*swro_L, J_d*J_r,'Color', lc,'LineWidth',lw);xlabel('x [m]','Fontsize',10); ylabel('[kg/sm]','Fontsize',10); ay=gca; ay.YAxis.Exponent = 0;xlim([0 swro_L]); hold on
yyaxis right
plot(x*swro_L, J_f*J_r,'Color', rc,'LineWidth',lw); ylabel('[kg/sm]','Fontsize',10); legend('J_{d}^{RO}(x)','J_{f}^{RO}(x)','Location','northeast');xlim([0 swro_L]); ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
% concentrations
subplot(3,2,3);lc='#0072bD';rc='#77AC30';
plot(x*swro_L, 100*C_d*C_r,'Color', lc,'LineWidth',lw);xlabel('x [m]','Fontsize',10); ylabel('[%]','Fontsize',10); ay=gca; ay.YAxis.Exponent = 0;hold on
yyaxis right; C_f=C_f(2:end);
plot(x2*swro_L, 100*C_f*C_r,'Color', rc,'LineWidth',lw); ylabel('[%]','Fontsize',10); legend('C_d^{RO}(x)','C_f^{RO}(x)','Location','southeast');xlim([0 swro_L]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
% Pressures
subplot(3,2,5);lc='#0072bD';rc='#77AC30';
plot(x*swro_L, P_d*p_r,'Color', lc,'LineWidth',lw);xlabel('x [m]','Fontsize',10); ylabel('[Pa]','Fontsize',10); ay=gca; ay.YAxis.Exponent = 5; hold on
yyaxis right
plot(x*swro_L, P_f*p_r,'Color', rc,'LineWidth',lw); ylabel('[Pa]','Fontsize',10); legend('P_d^{RO}(x)','P_f^{RO}(x)','Location','northeast');xlim([0 swro_L]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
% total mass flows
subplot(3,2,2);lc='#0072bD';rc='#77AC30';
plot(x*L,  Q_d*Q_r,'Color', lc,'LineWidth',lw);xlabel('x [m]','Fontsize',10); ylabel('[kg/sm]','Fontsize',10); ay=gca; ay.YAxis.Exponent = 0;xlim([0 swro_L]); hold on
yyaxis right
plot(x*L, Q_f*Q_r,'Color', rc,'LineWidth',lw); ylabel('[kg/sm]','Fontsize',10); legend('J_{d}^{PRO}(x)','J_{f}^{PRO}(x)','Location','northeast');xlim([0 L]); ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
% concentrations
subplot(3,2,4);lc='#0072bD';rc='#77AC30';
plot(x*L, 100*c_d*C_r,'Color', lc,'LineWidth',lw);xlabel('x [m]','Fontsize',10); ylabel('[%]','Fontsize',10); ay=gca; ay.YAxis(1).Exponent = 0;hold on
yyaxis right; c_f=c_f(2:end);
plot(x2*L, 100*c_f*C_r,'Color', rc,'LineWidth',lw); ylabel('[%]','Fontsize',10); legend('C_d^{PRO}(x)','C_f^{PRO}(x)','Location','northwest');xlim([0 L]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;ay.YAxis(2).Exponent = 0;
% Pressures
subplot(3,2,6);lc='#0072bD';rc='#77AC30';
plot(x*L, p_d*p_r,'Color', lc,'LineWidth',lw);xlabel('x [m]','Fontsize',10); ylabel('[Pa]','Fontsize',10);ay=gca; ay.YAxis.Exponent = 5; hold on
yyaxis right
plot(x*L, p_f*p_r,'Color', rc,'LineWidth',lw); ylabel('[Pa]','Fontsize',10); legend('P_d^{PRO}(x)','P_f^{PRO}(x)','Location','northeast');xlim([0 L]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
f=figure(2); %
f.Position = [199 232.3333 1.0987e+03 639.3333];
% Permeate flows 
subplot(2,2,1); lc='k';rc='#b81414'; J_cross2=J_cross(2:end);J_sin2=J_sin(2:end);
plot(x2*swro_L, J_cross2*J_r/swro_x_r, 'Color', lc, 'LineWidth',lw);xlabel('x [m]','Fontsize',10);  ylabel('[kg/sm^2]','Fontsize',10);ay=gca; ay.YAxis.Exponent = 0;ylim([min( J_cross2*J_r/swro_x_r) 1.1*max( J_cross2*J_r/swro_x_r)]); hold on
yyaxis right
plot(x2*swro_L, J_sin2*J_r/swro_x_r, 'Color', rc, 'LineWidth',lw); ylabel('[kg/sm^2]','Fontsize',10); legend('J_{w,in}^{RO}(x)','J_{s,in}^{RO}(x)','Location','northeast');xlim([0 swro_L]);ylim([min( J_sin2*J_r/swro_x_r) 1.1*max( J_sin2*J_r/swro_x_r)]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
% Osmotic/ Hydraulic pressure difference
subplot(2,2,3);lc='k';rc='#b81414'; osm_diff=swro_p_osm_d(2:end)-swro_p_osm_f(2:end); a=min([min(osm_diff*p_r) min((P_d-P_f)*p_r)]); b=max( [max(osm_diff*p_r), max((P_d-P_f)*p_r)]);
plot(x*swro_L, (P_d-P_f)*p_r,'Color', lc,'LineWidth',lw);xlabel('x [m]','Fontsize',10); ylabel('[Pa]','Fontsize',10); ylim([a b]);ay=gca; ay.YAxis.Exponent = 6; hold on
yyaxis right; 
plot(x2*swro_L, osm_diff*p_r,'Color', rc,'LineWidth',lw); ylabel('[Pa]','Fontsize',10); legend('\Delta P^{RO}(x)','\Delta \pi^{RO}(x)','Location','northeast');xlim([0 swro_L]);ylim([a b]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
% Permeate flows 
subplot(2,2,2); lc='k';rc='#b81414'; Q_cross2=Q_cross(2:end);Q_sin2=Q_sin(2:end);
plot(x2*L, -Q_cross2*Q_r/x_r, 'Color', lc, 'LineWidth',lw);xlabel('x [m]','Fontsize',10);  ylabel('[kg/sm^2]','Fontsize',10);ay=gca; ay.YAxis.Exponent = 0; hold on
yyaxis right
plot(x2*L, Q_sin2*Q_r/x_r, 'Color', rc, 'LineWidth',lw); ylabel('[kg/sm^2]','Fontsize',10); legend('J_{w,in}^{PRO}(x)','J_{s,in}^{PRO}(x)','Location','southwest');xlim([0 L]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
% Osmotic/ Hydraulic pressure difference
subplot(2,2,4);lc='k';rc='#b81414';osm_diff=p_osm_d(2:end)-p_osm_f(2:end); %a=min([min(osm_diff*p_r) min((p_d-p_f)*p_r)]); b=max( [max(osm_diff*p_r), max((p_d-p_f)*p_r)]);
plot(x*L, (p_d-p_f)*p_r,'Color', lc,'LineWidth',lw);xlabel('x [m]','Fontsize',10); ylabel('[Pa]','Fontsize',10); ay=gca;ay.YAxis.Exponent = 6; hold on; %ylim([a b]); 
yyaxis right; 
plot(x2*L, osm_diff*p_r,'Color', rc,'LineWidth',lw); ylabel('[Pa]','Fontsize',10); legend('\Delta P^{PRO}(x)','\Delta \pi^{PRO}(x)','Location','southeast');xlim([0 L]); ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc; %ylim([a b]);



end
%% figure 2
if fig(2) == 1
figure(2); % pressure+concentration in both PRO channels
subplot(2,2,1);plot(x*L, p_d*p_r,'b','LineWidth',lw);xlabel('x[m]','Fontsize',10); ylabel('[Pa]','Fontsize',10); legend('p_d(x)','Location','northeast'); xlim([0 L]);ay=gca; ay.YAxis.Exponent = 5;
subplot(2,2,2);plot(x*L, p_f*p_r,'g','LineWidth',lw);xlabel('x[m]','Fontsize',10); ylabel('[Pa]','Fontsize',10); legend('p_f(x)','Location','northeast');xlim([0 L]);ay=gca; ay.YAxis.Exponent = 5;
subplot(2,2,3);plot(x*L, 100*c_d*C_r,'b','LineWidth',lw);xlabel('x[m]','Fontsize',10); ylabel('[%]','Fontsize',10); legend('c_d(x)','Location','northeast');xlim([0 L]);ay=gca; ay.YAxis.Exponent = 0;
x2=x(2:end);c_f2=100.*c_f(2:end);
subplot(2,2,4);plot(x2*L, c_f2*C_r,'g','LineWidth',lw);xlabel('x[m]','Fontsize',10); ylabel('[%]','Fontsize',10); legend('c_f(x)','Location','northeast');xlim([0 L]);ay=gca; ay.YAxis.Exponent = 0;
end
%% figure 3
if fig(3) == 1
figure(3); % Q_s,in and Q_w,in
x2=x(2:end);Q_cross2=Q_cross(2:end);Q_sin2=Q_sin(2:end);
subplot(1,2,1);plot(x2*L, Q_cross2*Q_r/x_r, 'b', 'LineWidth',lw);xlabel('x [m]','Fontsize',10); ylabel('[kg/sm^2]','Fontsize',10); legend('Q_{w,in}(x)','Location','northeast');xlim([0 L]);ay=gca; ay.YAxis.Exponent = 0;
subplot(1,2,2);plot(x2*L, Q_sin2*Q_r/x_r, 'g', 'LineWidth',lw);xlabel('x [m]','Fontsize',10); ylabel('[kg/sm^2]','Fontsize',10); legend('Q_{s,in}(x)','Location','northeast');xlim([0 L]);ay=gca; ay.YAxis.Exponent = 0;
end
%% figure 4
lw=1.5; %Linewidth for all figures
if fig(4) == 1
figure(4); % Beta of PRO
    if version(2)==0
    plot(x*swro_L, swro_beta*J_r/swro_x_r,'b','LineWidth',lw);xlabel('x[m]','Fontsize',10); ylabel('[kg/sm^2]','Fontsize',10); legend('B(x)','Location','northeast');xlim([0 swro_L]);ylim([swro_beta(1)*J_r/swro_x_r*.9  swro_beta(1)*J_r/swro_x_r*1.1]); hold on
    else
    plot(x*swro_L, swro_beta*J_r/swro_x_r,'b','LineWidth',lw);xlabel('x[m]','Fontsize',10); ylabel('[kg/sm^2]','Fontsize',10); legend('B(x)','Location','northeast');xlim([0 swro_L]); hold on
    end
end
%% figure 5
if fig(5) == 1
figure(5); % pressure+concentration in both SWRO channels
subplot(2,2,1);plot(x*swro_L, P_d*p_r,'Color', '#0072bD','LineWidth',lw);xlabel('x[m]','Fontsize',10); ylabel('[Pa]','Fontsize',10); legend('P_d(x)','Location','northeast');xlim([0 swro_L]);ay=gca; ay.YAxis.Exponent = 5;
x2=x(2:end);C_f2=100.*C_f(2:end);
subplot(2,2,2);plot(x*swro_L, P_f*p_r,'Color', '#77AC30','LineWidth',lw);xlabel('x[m]','Fontsize',10); ylabel('[Pa]','Fontsize',10); legend('P_f(x)','Location','northeast');xlim([0 swro_L]);ay=gca; ay.YAxis.Exponent = 5;
subplot(2,2,4);plot(x2*swro_L, C_f2*C_r,'Color', '#77AC30','LineWidth',lw);xlabel('x[m]','Fontsize',10); ylabel('[%]','Fontsize',10); legend('C_f(x)','Location','southeast');xlim([0 swro_L]);ay=gca; ay.YAxis.Exponent = 0;
subplot(2,2,3);plot(x*swro_L, 100*C_d*C_r,'Color', '#0072bD','LineWidth',lw);xlabel('x[m]','Fontsize',10); ylabel('[%]','Fontsize',10); legend('C_d(x)','Location','southeast');xlim([0 swro_L]);ay=gca; ay.YAxis.Exponent = 0;
end
%% figure 6
if fig(6) == 1
figure(6); % J_s,in and J_w,in
x2=x(2:end);J_cross2=J_cross(2:end);J_sin2=J_sin(2:end);
subplot(1,2,2);plot(x2*swro_L, J_sin2*J_r/swro_x_r, 'Color', '#77AC30', 'LineWidth',lw);xlabel('x[m]','Fontsize',10); ylabel('[kg/sm^2]','Fontsize',10); legend('J_{s,in}(x)','Location','northeast');xlim([0 swro_L]);ay=gca; ay.YAxis.Exponent = 0;
subplot(1,2,1);plot(x2*swro_L, J_cross2*J_r/swro_x_r, 'Color', '#0072bD', 'LineWidth',lw);xlabel('x[m]','Fontsize',10);  ylabel('[kg/sm^2]','Fontsize',10); legend('J_{w,in}(x)','Location','northeast');xlim([0 swro_L]);ay=gca; ay.YAxis.Exponent = 0;
end
%% figure 7
if fig(7) == 1
f=figure(7); % PRO quantities
f.Position = [100 1 1100 1000];
subplot(4,3,4);plot(x*L, Q_sd*Q_r,'b','LineWidth',lw);xlabel('x[m]','Fontsize',10); ylabel('[kg/sm]','Fontsize',10); legend('J_{s,d}^{PRO}(x)','Location','northeast');xlim([0 L]);ay=gca; ay.YAxis.Exponent = 0;
subplot(4,3,1);plot(x*L, Q_wd*Q_r,'b','LineWidth',lw);xlabel('x[m]','Fontsize',10); ylabel('[kg/sm]','Fontsize',10); legend('J_{w,d}^{PRO}(x)','Location','northeast');xlim([0 L]);ay=gca; ay.YAxis.Exponent = 0;
subplot(4,3,5);plot(x*L, Q_sf*Q_r,'g','LineWidth',lw);xlabel('x[m]','Fontsize',10); ylabel('[kg/sm]','Fontsize',10); legend('J_{s,f}^{PRO}(x)','Location','northeast');xlim([0 L]);ay=gca; ay.YAxis.Exponent = 0;
subplot(4,3,2);plot(x*L, Q_wf*Q_r,'g','LineWidth',lw);xlabel('x[m]','Fontsize',10); ylabel('[kg/sm]','Fontsize',10); legend('J_{w,f}^{PRO}(x)','Location','northeast');xlim([0 L]);ay=gca; ay.YAxis.Exponent = 0; 
subplot(4,3,7);plot(x*L, p_d*p_r,'b','LineWidth',lw);xlabel('x[m]','Fontsize',10); ylabel('[Pa]','Fontsize',10); legend('P_d^{PRO}(x)','Location','northeast'); xlim([0 L]);ay=gca; ay.YAxis.Exponent = 5;
subplot(4,3,8);plot(x*L, p_f*p_r,'g','LineWidth',lw);xlabel('x[m]','Fontsize',10); ylabel('[Pa]','Fontsize',10); legend('P_f^{PRO}(x)','Location','northeast');xlim([0 L]);ay=gca; ay.YAxis.Exponent = 5;
subplot(4,3,10);plot(x*L, 100*c_d*C_r,'b','LineWidth',lw);xlabel('x[m]','Fontsize',10); ylabel('[%]','Fontsize',10); legend('C_d^{PRO}(x)','Location','northeast');xlim([0 L]);ay=gca; ay.YAxis.Exponent = 0;
x2=x(2:end);c_f=100.*c_f(2:end);Q_cross2=Q_cross(2:end);Q_sin2=Q_sin(2:end);
subplot(4,3,11);plot(x2*L, c_f*C_r,'g','LineWidth',lw);xlabel('x[m]','Fontsize',10); ylabel('[%]','Fontsize',10); legend('C_f^{PRO}(x)','Location','northeast');xlim([0 L]);ay=gca; ay.YAxis.Exponent = 0;
subplot(4,3,[3 6]);plot(x2*L, Q_cross2*Q_r/x_r, 'b', 'LineWidth',lw);xlabel('x [m]','Fontsize',10); ylabel('[kg/sm^2]','Fontsize',10); legend('J_{w,in}^{PRO}(x)','Location','northeast');xlim([0 L]);ay=gca; ay.YAxis.Exponent = 0;
subplot(4,3,[9 12]);plot(x2*L, Q_sin2*Q_r/x_r, 'g', 'LineWidth',lw);xlabel('x [m]','Fontsize',10); ylabel('[kg/sm^2]','Fontsize',10); legend('J_{s,in}^{PRO}(x)','Location','northeast');xlim([0 L]);ay=gca; ay.YAxis.Exponent = 0;
end
%% figure 8
if fig(8) == 1
f=figure(8); % SWRO quantities
f.Position = [100 1 1100 1000];
x2=x(2:end);J_wf2=J_wf(2:end);J_sf2=J_sf(2:end);
subplot(4,3,1);plot(x*swro_L, J_wd*J_r,'Color', '#0072bD','LineWidth',lw);xlabel('x[m]','Fontsize',10); ylabel('[kg/sm]','Fontsize',10); legend('J_{w,d}^{RO}(x)','Location','northeast');xlim([0 swro_L]);ay=gca; ay.YAxis.Exponent = 0;
subplot(4,3,4);plot(x*swro_L, J_sd*J_r,'Color', '#0072bD','LineWidth',lw);xlabel('x[m]','Fontsize',10); ylabel('[kg/sm]','Fontsize',10); legend('J_{s,d}^{RO}(x)','Location','northeast');xlim([0 swro_L]);ay=gca; ay.YAxis.Exponent = 0;
subplot(4,3,2);plot(x2*swro_L, J_wf2*J_r,'Color', '#77AC30','LineWidth',lw);xlabel('x[m]','Fontsize',10); ylabel('[kg/sm]','Fontsize',10); legend('J_{w,f}^{RO}(x)','Location','southeast');xlim([0 swro_L]);ay=gca; ay.YAxis.Exponent = 0;
subplot(4,3,5);plot(x2*swro_L, J_sf2*J_r,'Color', '#77AC30','LineWidth',lw);xlabel('x[m]','Fontsize',10); ylabel('[kg/sm]','Fontsize',10); legend('J_{s,f}^{RO}(x)','Location','southeast');xlim([0 swro_L]);ay=gca; ay.YAxis.Exponent = 0;     
subplot(4,3,7);plot(x*swro_L, P_d*p_r,'Color', '#0072bD','LineWidth',lw);xlabel('x[m]','Fontsize',10); ylabel('[Pa]','Fontsize',10); legend('P_d^{RO}(x)','Location','northeast');xlim([0 swro_L]);ay=gca; ay.YAxis.Exponent = 5;
subplot(4,3,8);plot(x*swro_L, P_f*p_r,'Color', '#77AC30','LineWidth',lw);xlabel('x[m]','Fontsize',10); ylabel('[Pa]','Fontsize',10); legend('P_f^{RO}(x)','Location','northeast');xlim([0 swro_L]);ay=gca; ay.YAxis.Exponent = 5;
subplot(4,3,10);plot(x*swro_L, 100*C_d*C_r,'Color', '#0072bD','LineWidth',lw);xlabel('x[m]','Fontsize',10); ylabel('[%]','Fontsize',10); legend('C_d^{RO}(x)','Location','southeast');xlim([0 swro_L]);ay=gca; ay.YAxis.Exponent = 0;
C_f=100.*C_f(2:end);J_cross2=J_cross(2:end);J_sin2=J_sin(2:end);
subplot(4,3,11);plot(x2*swro_L, C_f*C_r,'Color', '#77AC30','LineWidth',lw);xlabel('x[m]','Fontsize',10); ylabel('[%]','Fontsize',10); legend('C_f^{RO}(x)','Location','southeast');xlim([0 swro_L]);ay=gca; ay.YAxis.Exponent = 0;
subplot(4,3,[3 6]);plot(x2*swro_L, J_cross2*J_r/swro_x_r, 'Color', '#0072bD', 'LineWidth',lw);xlabel('x[m]','Fontsize',10);  ylabel('[kg/sm^2]','Fontsize',10); legend('J_{w,in}^{RO}(x)','Location','northeast');xlim([0 swro_L]);ay=gca; ay.YAxis.Exponent = 0;
subplot(4,3,[9 12]);plot(x2*swro_L, J_sin2*J_r/swro_x_r, 'Color', '#77AC30', 'LineWidth',lw);xlabel('x[m]','Fontsize',10); ylabel('[kg/sm^2]','Fontsize',10); legend('J_{s,in}^{RO}(x)','Location','northeast');xlim([0 swro_L]);ay=gca; ay.YAxis.Exponent = 0;
end
%% figure 9
if fig(9) == 1
f=figure(3);  % ERD quantities
f.Position = [100 1 900 500]; tiledlayout(1,4); 
rot=0; rx=0; %rotation and shift of texts
vals1 = [pE*p_r; pERD*p_r; 0; P_d(end)*p_r; p_d(end)*p_r];
vals2 = 100*[cE; C_ERD; 0; C_d(end); c_d(1)];
vals3 = [J_E*J_r; J_ERD*J_r; 0; J_d(end)*J_r; Q_d(1)*Q_r];
%adjust for counter current / SWRO+ERD / 2ERDs
if version(1)==1; vals2(5) = 100*c_d(end); vals3(5) =-Q_d(end)*Q_r;end
if version(6)==1; vals1(5) = p_exit*p_r; vals2(5) = 100*c_exit; vals3(5) = Q_exit*J_r; end
if version(6)==4; vals1(1)= pERD2*p_r; vals2(1)=100*C_ERD2; end
nexttile
b = bar(1,vals1, 'FaceColor', 'b'); b(4).FaceColor = "#77AC30"; b(5).FaceColor = [.2 .6 .5]; b(1).FaceColor = '#0072bD';
xticklabels({'[kg/ms^2]'}); ylim([0 1.1*max(vals1)]);ay=gca;ay.YAxis(1).Exponent=5;
t=text(b(1).XEndPoints,b(1).YEndPoints,"P_{s}^{in}",'HorizontalAlignment','center','VerticalAlignment','bottom');
t=text(b(2).XEndPoints+rx ,b(2).YEndPoints,"P_{s}^{out}",'HorizontalAlignment','center','VerticalAlignment','bottom'); 
t=text(b(4).XEndPoints,b(4).YEndPoints,"P_{b}^{in}",'HorizontalAlignment','center','VerticalAlignment','bottom'); 
t=text(b(5).XEndPoints,b(5).YEndPoints,"P_{b}^{out}",'HorizontalAlignment','center','VerticalAlignment','bottom'); 
nexttile
b = bar(1,vals2, 'FaceColor', 'b'); b(4).FaceColor = "#77AC30"; b(5).FaceColor = [.2 .6 .5]; b(1).FaceColor = '#0072bD';
xticklabels({'[%]'}); ylim([0 1.12*max(vals2)]);ay=gca;ay.YAxis(1).Exponent=0;
t= text(b(1).XEndPoints,b(1).YEndPoints,"C_{s}^{in}",'HorizontalAlignment','center','VerticalAlignment','bottom'); t.Rotation = rot;
t=text(b(2).XEndPoints+rx ,b(2).YEndPoints,"C_{s}^{out}",'HorizontalAlignment','center','VerticalAlignment','bottom'); t.Rotation = rot;
t=text(b(4).XEndPoints,b(4).YEndPoints,"C_b^{in}",'HorizontalAlignment','center','VerticalAlignment','bottom'); t.Rotation = rot;
t=text(b(5).XEndPoints,b(5).YEndPoints,"C_{b}^{out}",'HorizontalAlignment','center','VerticalAlignment','bottom'); t.Rotation = rot;
nexttile
b = bar(1,vals3, 'FaceColor', 'b'); b(4).FaceColor = "#77AC30"; b(5).FaceColor = [.2 .6 .5]; b(1).FaceColor = '#0072bD';
xticklabels({'[kg/ms]'}); ylim([0 1.1*max(vals3)]);ay=gca;ay.YAxis(1).Exponent=0;
t=text(b(1).XEndPoints,b(1).YEndPoints,"J_{s}^{in}",'HorizontalAlignment','center','VerticalAlignment','bottom'); t.Rotation = rot;
t=text(b(2).XEndPoints+rx ,b(2).YEndPoints,"J_{s}^{out}",'HorizontalAlignment','center','VerticalAlignment','bottom'); t.Rotation = rot;
t=text(b(4).XEndPoints ,b(4).YEndPoints,"J_{b}^{in}",'HorizontalAlignment','center','VerticalAlignment','bottom'); t.Rotation = rot;
t=text(b(5).XEndPoints,b(5).YEndPoints,"J_{b}^{out}",'HorizontalAlignment','center','VerticalAlignment','bottom'); t.Rotation = rot ;
nexttile
if version(6)>3
lgd = legend([b(1) b(2) b(5) b(4)],'LP ERD 1 Seawater Inlet','HP ERD 1 Seawater Outlet','HP SWRO Brine Outlet','LP PRO Brine Inlet' , 'Location', 'EastOutside');
else
lgd = legend([b(1) b(2) b(5) b(4)],'Sea water inlet','Sea water outlet','Brine outlet','Brine inlet' , 'Location', 'EastOutside');
end
lgd.Layout.Tile = 4; axis off ; 
end
%% figure 10
if fig(9) == 1 && version(6) ==4
close all;f=figure(3);  % ERD quantities
f.Position = [950 1 900 500]; tiledlayout(1,4); 
rot=0; rx=0; %rotation and shift of texts
vals1 = [pE*p_r; pERD2*p_r; 0; p_d(end)*p_r; p_exit*p_r];
vals2 = 100*[cE; C_ERD2; 0; c_d(end); c_exit];
vals3 = [J_E_sea*J_r; 2*J_E*J_r; 0; mixer_ERD*Q_d(end)*Q_r; Q_exit*Q_r];
%adjust for counter current
if version(1)==1; vals1(4) = p_d(1)*p_r; vals2(4)=c_d(1)*100;vals3(4)=-mixer_ERD*Q_d(1)*Q_r; end
nexttile
b = bar(1,vals1, 'FaceColor', 'b'); b(4).FaceColor = "#77AC30"; b(5).FaceColor = [.2 .6 .5]; b(1).FaceColor = '#0072bD';
xticklabels({'[kg/ms^2]'}); ylim([0 1.1*max(vals1)]);ay=gca;ay.YAxis(1).Exponent=5;
t=text(b(1).XEndPoints,b(1).YEndPoints,"P_E",'HorizontalAlignment','center','VerticalAlignment','bottom');
t=text(b(2).XEndPoints+rx ,b(2).YEndPoints,"P_{ERD 2}",'HorizontalAlignment','center','VerticalAlignment','bottom'); 
if version(1) == 0; t=text(b(4).XEndPoints,b(4).YEndPoints,"p_d^1",'HorizontalAlignment','center','VerticalAlignment','bottom');  end
if version(1) == 1; t=text(b(4).XEndPoints,b(4).YEndPoints,"p_d^0",'HorizontalAlignment','center','VerticalAlignment','bottom');  end
t=text(b(5).XEndPoints,b(5).YEndPoints,"P_E",'HorizontalAlignment','center','VerticalAlignment','bottom'); 
nexttile
b = bar(1,vals2, 'FaceColor', 'b'); b(4).FaceColor = "#77AC30"; b(5).FaceColor = [.2 .6 .5]; b(1).FaceColor = '#0072bD';
xticklabels({'[%]'}); ylim([0 1.12*max(vals2)]);ay=gca;ay.YAxis(1).Exponent=0;
t= text(b(1).XEndPoints,b(1).YEndPoints,"C_E",'HorizontalAlignment','center','VerticalAlignment','bottom'); t.Rotation = rot;
t=text(b(2).XEndPoints+rx ,b(2).YEndPoints," C_{ERD 2}",'HorizontalAlignment','center','VerticalAlignment','bottom'); t.Rotation = rot;
if version(1) == 0; t=text(b(4).XEndPoints,b(4).YEndPoints,"c_d^1",'HorizontalAlignment','center','VerticalAlignment','bottom'); t.Rotation = rot; end
if version(1) == 1; t=text(b(4).XEndPoints,b(4).YEndPoints,"c_d^0",'HorizontalAlignment','center','VerticalAlignment','bottom'); t.Rotation = rot; end
t=text(b(5).XEndPoints,b(5).YEndPoints,"C_{exit}",'HorizontalAlignment','center','VerticalAlignment','bottom'); t.Rotation = rot; 
nexttile
b = bar(1,vals3, 'FaceColor', 'b'); b(4).FaceColor = "#77AC30"; b(5).FaceColor = [.2 .6 .5]; b(1).FaceColor = '#0072bD';
xticklabels({'[kg/ms]'}); ylim([0 1.1*max(vals3)]);ay=gca;ay.YAxis(1).Exponent=0;
t=text(b(1).XEndPoints,b(1).YEndPoints,"J_{sea}",'HorizontalAlignment','center','VerticalAlignment','bottom'); t.Rotation = rot;
t=text(b(2).XEndPoints+rx ,b(2).YEndPoints,"   2J_{E}",'HorizontalAlignment','center','VerticalAlignment','bottom'); t.Rotation = rot;
strg=[num2str(mixer_ERD), '*'];
if version(1) == 0; t=text(b(4).XEndPoints ,b(4).YEndPoints,[strg,"$Q_d^1$"],'Interpreter','latex','HorizontalAlignment','center','VerticalAlignment','bottom'); t.Rotation = rot;end
if version(1) == 1; t=text(b(4).XEndPoints ,b(4).YEndPoints,[strg,"$-Q_d^0$"],'Interpreter','latex','HorizontalAlignment','center','VerticalAlignment','bottom'); t.Rotation = rot;end
t=text(b(5).XEndPoints,b(5).YEndPoints,"Q_{exit}",'HorizontalAlignment','center','VerticalAlignment','bottom'); t.Rotation = rot;
nexttile
lgd = legend([b(1) b(2) b(3) b(4)],'LP ERD 2 Seawater Inlet','HP ERD 2 Seawater Outlet','HP ERD 2 Brine Inlet','LP ERD 2 Brine Outlet' , 'Location', 'EastOutside');
lgd.Layout.Tile = 4; axis off ; 
end
end
%% figures for Paper
if fig(11) == 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=figure(1); % SWRO quantities
f.Position = [276.3333 196.3333 1100 1000];
x2=x(2:end);
% total mass flows
subplot(3,2,1);lc='#0072bD';rc='#77AC30';
plot(x*swro_L, J_d*J_r,'Color', lc,'LineWidth',lw);xlabel('x[m]','Fontsize',10); ylabel('[kg/sm]','Fontsize',10); ay=gca; ay.YAxis.Exponent = 0;xlim([0 swro_L]); hold on
yyaxis right
plot(x*swro_L, J_f*J_r,'Color', rc,'LineWidth',lw); ylabel('[kg/sm]','Fontsize',10); legend('J_{d}^{RO}(x)','J_{f}^{RO}(x)','Location','northeast');xlim([0 swro_L]); ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
% concentrations
subplot(3,2,2);lc='#0072bD';rc='#77AC30';
plot(x*swro_L, 100*C_d*C_r,'Color', lc,'LineWidth',lw);xlabel('x[m]','Fontsize',10); ylabel('[%]','Fontsize',10); ay=gca; ay.YAxis.Exponent = 0;hold on
yyaxis right; C_f=C_f(2:end);
plot(x2*swro_L, 100*C_f*C_r,'Color', rc,'LineWidth',lw); ylabel('[%]','Fontsize',10); legend('C_d^{RO}(x)','C_f^{RO}(x)','Location','southeast');xlim([0 swro_L]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
% Pressures
subplot(3,2,3);lc='#0072bD';rc='#77AC30';
plot(x*swro_L, P_d*p_r,'Color', lc,'LineWidth',lw);xlabel('x[m]','Fontsize',10); ylabel('[Pa]','Fontsize',10); ay=gca; ay.YAxis.Exponent = 5; hold on
yyaxis right
plot(x*swro_L, P_f*p_r,'Color', rc,'LineWidth',lw); ylabel('[Pa]','Fontsize',10); legend('P_d^{RO}(x)','P_f^{RO}(x)','Location','northeast');xlim([0 swro_L]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
% Permeate flows 
subplot(3,2,5); lc='k';rc='#b81414'; J_cross2=J_cross(2:end);J_sin2=J_sin(2:end);
plot(x2*swro_L, J_cross2*J_r/swro_x_r, 'Color', lc, 'LineWidth',lw);xlabel('x[m]','Fontsize',10);  ylabel('[kg/sm^2]','Fontsize',10);ay=gca; ay.YAxis.Exponent = 0;ylim([min( J_cross2*J_r/swro_x_r) 1.1*max( J_cross2*J_r/swro_x_r)]); hold on
yyaxis right
plot(x2*swro_L, J_sin2*J_r/swro_x_r, 'Color', rc, 'LineWidth',lw); ylabel('[kg/sm^2]','Fontsize',10); legend('J_{w,in}^{RO}(x)','J_{s,in}^{RO}(x)','Location','northeast');xlim([0 swro_L]);ylim([min( J_sin2*J_r/swro_x_r) 1.1*max( J_sin2*J_r/swro_x_r)]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
% Osmotic/ Hydraulic pressure difference
subplot(3,2,6);lc='k';rc='#b81414'; osm_diff=swro_p_osm_d(2:end)-swro_p_osm_f(2:end); a=min([min(osm_diff*p_r) min((P_d-P_f)*p_r)]); b=max( [max(osm_diff*p_r), max((P_d-P_f)*p_r)]);
plot(x*swro_L, (P_d-P_f)*p_r,'Color', lc,'LineWidth',lw);xlabel('x[m]','Fontsize',10); ylabel('[Pa]','Fontsize',10); ylim([a b]);ay=gca; ay.YAxis.Exponent = 6; hold on
yyaxis right; 
plot(x2*swro_L, osm_diff*p_r,'Color', rc,'LineWidth',lw); ylabel('[Pa]','Fontsize',10); legend('\Delta P^{RO}(x)','\Delta \pi^{RO}(x)','Location','northeast');xlim([0 swro_L]);ylim([a b]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=figure(2); % PRO quantities
f.Position = [936.3333 266.3333 1.0987e+03 632];
x2=x(2:end);
% total mass flows
subplot(3,2,1);lc='#0072bD';rc='#77AC30';
plot(x*L,  Q_d*Q_r,'Color', lc,'LineWidth',lw);xlabel('x[m]','Fontsize',10); ylabel('[kg/sm]','Fontsize',10); ay=gca; ay.YAxis.Exponent = 0;xlim([0 swro_L]); hold on
yyaxis right
plot(x*L, Q_f*Q_r,'Color', rc,'LineWidth',lw); ylabel('[kg/sm]','Fontsize',10); legend('J_{d}^{PRO}(x)','J_{f}^{PRO}(x)','Location','northeast');xlim([0 L]); ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
% concentrations
subplot(3,2,2);lc='#0072bD';rc='#77AC30';
plot(x*L, 100*c_d*C_r,'Color', lc,'LineWidth',lw);xlabel('x[m]','Fontsize',10); ylabel('[%]','Fontsize',10); ay=gca; ay.YAxis(1).Exponent = 0;hold on
yyaxis right; c_f=c_f(2:end);
plot(x2*L, 100*c_f*C_r,'Color', rc,'LineWidth',lw); ylabel('[%]','Fontsize',10); legend('C_d^{PRO}(x)','C_f^{PRO}(x)','Location','southeast');xlim([0 L]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;ay.YAxis(2).Exponent = 0;
% Pressures
subplot(3,2,3);lc='#0072bD';rc='#77AC30';
plot(x*L, p_d*p_r,'Color', lc,'LineWidth',lw);xlabel('x[m]','Fontsize',10); ylabel('[Pa]','Fontsize',10);ay=gca; ay.YAxis.Exponent = 5; hold on
yyaxis right
plot(x*L, p_f*p_r,'Color', rc,'LineWidth',lw); ylabel('[Pa]','Fontsize',10); legend('P_d^{PRO}(x)','P_f^{PRO}(x)','Location','northeast');xlim([0 L]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
% Permeate flows 
subplot(3,2,5); lc='k';rc='#b81414'; Q_cross2=Q_cross(2:end);Q_sin2=Q_sin(2:end);
plot(x2*L, Q_cross2*Q_r/x_r, 'Color', lc, 'LineWidth',lw);xlabel('x[m]','Fontsize',10);  ylabel('[kg/sm^2]','Fontsize',10);ay=gca; ay.YAxis.Exponent = 0; hold on
yyaxis right
plot(x2*L, Q_sin2*Q_r/x_r, 'Color', rc, 'LineWidth',lw); ylabel('[kg/sm^2]','Fontsize',10); legend('J_{w,in}^{PRO}(x)','J_{s,in}^{PRO}(x)','Location','northeast');xlim([0 L]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
% Osmotic/ Hydraulic pressure difference
subplot(3,2,6);lc='k';rc='#b81414';osm_diff=p_osm_d(2:end)-p_osm_f(2:end); %a=min([min(osm_diff*p_r) min((p_d-p_f)*p_r)]); b=max( [max(osm_diff*p_r), max((p_d-p_f)*p_r)]);
plot(x*L, (p_d-p_f)*p_r,'Color', lc,'LineWidth',lw);xlabel('x[m]','Fontsize',10); ylabel('[Pa]','Fontsize',10); ay=gca;ay.YAxis.Exponent = 6; hold on; %ylim([a b]); 
yyaxis right; 
plot(x2*L, osm_diff*p_r,'Color', rc,'LineWidth',lw); ylabel('[Pa]','Fontsize',10); legend('\Delta P^{PRO}(x)','\Delta \pi^{PRO}(x)','Location','northeast');xlim([0 L]); ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc; %ylim([a b]);
end
%% figures for Paper - all in individual plot
if fig(11) == 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=figure(1);% SWRO quantities
x2=x(2:end);
% total mass flows
lc='#0072bD';rc='#77AC30';
plot(x*swro_L, J_d*J_r,'Color', lc,'LineWidth',lw); ylabel('[kg/sm]','Fontsize',10); ay=gca; ay.YAxis.Exponent = 0; xlim([0 swro_L]); hold on
yyaxis right
plot(x*swro_L, J_f*J_r,'Color', rc,'LineWidth',lw); ylabel('[kg/sm]','Fontsize',10); legend('J_{d}^{RO}(x)','J_{f}^{RO}(x)','Location','northeast');xlim([0 swro_L]); ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
set(gcf,'color','w'); f = gcf; exportgraphics(f,'Simulation_SWROa.png');
% concentrations
close all; f=figure(1);
lc='#0072bD';rc='#77AC30';
plot(x*swro_L, 100*C_d*C_r,'Color', lc,'LineWidth',lw); ylabel('[%]','Fontsize',10); ay=gca; ay.YAxis.Exponent = 0;hold on
yyaxis right; C_f=C_f(2:end);
plot(x2*swro_L, 100*C_f*C_r,'Color', rc,'LineWidth',lw); ylabel('[%]','Fontsize',10); legend('C_d^{RO}(x)','C_f^{RO}(x)','Location','southeast');xlim([0 swro_L]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
set(gcf,'color','w'); f = gcf; exportgraphics(f,'Simulation_SWROb.png');
% Pressures
close all; f=figure(6);
lc='#0072bD';rc='#77AC30';
plot(x*swro_L, P_d*p_r,'Color', lc,'LineWidth',lw);ylabel('[Pa]','Fontsize',10); ay=gca; ay.YAxis.Exponent = 5; hold on
yyaxis right
plot(x*swro_L, P_f*p_r,'Color', rc,'LineWidth',lw); ylabel('[Pa]','Fontsize',10); legend('P_d^{RO}(x)','P_f^{RO}(x)','Location','northeast');xlim([0 swro_L]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
set(gcf,'color','w'); f = gcf; exportgraphics(f,'Simulation_SWROc.png');
% Permeate flows 
close all;f=figure(1);
lc='k';rc='#b81414'; J_cross2=J_cross(2:end);J_sin2=J_sin(2:end);
plot(x2*swro_L, J_cross2*J_r/swro_x_r, 'Color', lc, 'LineWidth',lw); ylabel('[kg/sm^2]','Fontsize',10);ay=gca; ay.YAxis.Exponent = 0;ylim([min( J_cross2*J_r/swro_x_r) 1.1*max( J_cross2*J_r/swro_x_r)]); hold on
yyaxis right
plot(x2*swro_L, J_sin2*J_r/swro_x_r, 'Color', rc, 'LineWidth',lw); ylabel('[kg/sm^2]','Fontsize',10); legend('J_{w,in}^{RO}(x)','J_{s,in}^{RO}(x)','Location','northeast');xlim([0 swro_L]);ylim([min( J_sin2*J_r/swro_x_r) 1.1*max( J_sin2*J_r/swro_x_r)]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
set(gcf,'color','w'); f = gcf; exportgraphics(f,'Simulation_SWROd.png');
% Osmotic/ Hydraulic pressure difference
close all;f=figure(1);
lc='k';rc='#b81414'; osm_diff=swro_p_osm_d(2:end)-swro_p_osm_f(2:end); a=min([min(osm_diff*p_r) min((P_d-P_f)*p_r)]); b=max( [max(osm_diff*p_r), max((P_d-P_f)*p_r)]);
plot(x*swro_L, (P_d-P_f)*p_r,'Color', lc,'LineWidth',lw);ylabel('[Pa]','Fontsize',10); ylim([a b]);ay=gca; ay.YAxis.Exponent = 6; hold on
yyaxis right; 
plot(x2*swro_L, osm_diff*p_r,'Color', rc,'LineWidth',lw); ylabel('[Pa]','Fontsize',10); legend('\Delta P^{RO}(x)','\Delta \pi^{RO}(x)','Location','northeast');xlim([0 swro_L]);ylim([a b]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
set(gcf,'color','w'); f = gcf; exportgraphics(f,'Simulation_SWROe.png');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=figure(2); % PRO quantities
x2=x(2:end);
% total mass flows
lc='#0072bD';rc='#77AC30';
plot(x*L,  Q_d*Q_r,'Color', lc,'LineWidth',lw);ylabel('[kg/sm]','Fontsize',10); ay=gca; ay.YAxis.Exponent = 0;xlim([0 swro_L]); hold on
yyaxis right
plot(x*L, Q_f*Q_r,'Color', rc,'LineWidth',lw); ylabel('[kg/sm]','Fontsize',10); legend('J_{d}^{PRO}(x)','J_{f}^{PRO}(x)','Location','northeast');xlim([0 L]); ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
set(gcf,'color','w'); f = gcf; exportgraphics(f,'Simulation_PROa.png');
% concentrations
close all; f=figure(2);lc='#0072bD';rc='#77AC30';
plot(x*L, 100*c_d*C_r,'Color', lc,'LineWidth',lw); ylabel('[%]','Fontsize',10); ay=gca; ay.YAxis(1).Exponent = 0;hold on
yyaxis right; c_f=c_f(2:end);
plot(x2*L, 100*c_f*C_r,'Color', rc,'LineWidth',lw); ylabel('[%]','Fontsize',10); legend('C_d^{PRO}(x)','C_f^{PRO}(x)','Location','southeast');xlim([0 L]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;ay.YAxis(2).Exponent = 0;
set(gcf,'color','w'); f = gcf; exportgraphics(f,'Simulation_PROb.png');
% Pressures
close all; f=figure(2);lc='#0072bD';rc='#77AC30';
plot(x*L, p_d*p_r,'Color', lc,'LineWidth',lw);ylabel('[Pa]','Fontsize',10);ay=gca; ay.YAxis.Exponent = 5; hold on
yyaxis right
plot(x*L, p_f*p_r,'Color', rc,'LineWidth',lw); ylabel('[Pa]','Fontsize',10); legend('P_d^{PRO}(x)','P_f^{PRO}(x)','Location','northeast');xlim([0 L]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
set(gcf,'color','w'); f = gcf; exportgraphics(f,'Simulation_PROc.png');
% Permeate flows 
close all; f=figure(2);lc='k';rc='#b81414'; Q_cross2=Q_cross(2:end);Q_sin2=Q_sin(2:end);
plot(x2*L, Q_cross2*Q_r/x_r, 'Color', lc, 'LineWidth',lw); ylabel('[kg/sm^2]','Fontsize',10);ay=gca; ay.YAxis.Exponent = 0; hold on
yyaxis right
plot(x2*L, Q_sin2*Q_r/x_r, 'Color', rc, 'LineWidth',lw); ylabel('[kg/sm^2]','Fontsize',10); legend('J_{w,in}^{PRO}(x)','J_{s,in}^{PRO}(x)','Location','northeast');xlim([0 L]);ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc;
set(gcf,'color','w'); f = gcf; exportgraphics(f,'Simulation_PROd.png');
% Osmotic/ Hydraulic pressure difference
close all; f=figure(2);lc='k';rc='#b81414';osm_diff=p_osm_d(2:end)-p_osm_f(2:end); %a=min([min(osm_diff*p_r) min((p_d-p_f)*p_r)]); b=max( [max(osm_diff*p_r), max((p_d-p_f)*p_r)]);
plot(x*L, (p_d-p_f)*p_r,'Color', lc,'LineWidth',lw);ylabel('[Pa]','Fontsize',10); ay=gca;ay.YAxis.Exponent = 6; hold on; %ylim([a b]); 
yyaxis right; 
plot(x2*L, osm_diff*p_r,'Color', rc,'LineWidth',lw); ylabel('[Pa]','Fontsize',10); legend('\Delta P^{PRO}(x)','\Delta \pi^{PRO}(x)','Location','northeast');xlim([0 L]); ay.YAxis(1).Color = lc; ay.YAxis(2).Color = rc; %ylim([a b]);
set(gcf,'color','w'); f = gcf; exportgraphics(f,'Simulation_PROe.png');
end


%ax = axes; yyaxis('left'); yyaxis('right');
%ax.YAxis(1).Color = [0 0 0]; ax.YAxis(2).Color = [0 0 0];


























end
end