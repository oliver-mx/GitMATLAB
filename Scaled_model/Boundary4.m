function [ res ] = Boundary4(ya, yb, DATA)
% Residual function for boundary conditions
% version(6) == 4
% --> SWRO-PRO hybrid system (with two ERD)
% important: mixer_ERD=0; is not allowed! (use hybrid system I instead)

%% Read data
[H, Z, swro_Z, ro_water, ro_salt, Mw, Ms, Rw, T0, eta, sigma, p_r, rho_r, C_r, swro_L, swro_alpha, swro_KK, swro_x_r, swro_b1, swro_b2, J_r, swro_gamma, swro_gamma2, swro_W_r, L, alpha, KK, x_r, b1, b2, Q_r, gamma, gamma2, W_r, cE, pE, rho_E, J_sf_0, J_wf_0, Pd_0, Pd_L, Pf_L, Q_sf_0, pd_0, pf_0, pd_L, pf_L, HP_eff, LP_eff, T_eff, V_m, ERD_eff, ERD_fric, A_ERD, eta_ERD, mix_density, pw, pe, swro_beta_fix, beta_fix, mix_M1, mix_M3, version, fig, swro_KF, swro_KD, KF, KD]...
    = DATA(1);

%% no Mixing i.e. V_m = 0
if V_m ==0 
swroC_in = cE;
proC_in = yb(1);
PRO_brine = J_r/Q_r*(yb(1)*yb(2)+yb(2))*(1-eta_ERD);
else
%% with Mixing i.e. V_m > 0
% second ERD (i.e. index 2)
if version(1)==0
    rho_d2= real((yb(6+1)+1)./(yb(6+1)/ro_salt + 1/ro_water));
else
    rho_d2= real((ya(6+1)+1)./(ya(6+1)/ro_salt + 1/ro_water));
end
rho_ERD2 = V_m*(rho_d2 - rho_E)+rho_E;
C_ERD2= -ro_salt*(rho_ERD2-ro_water)/(ro_water*(rho_ERD2-ro_salt));
% first ERD (i.e. index 1)
rho_d1= real((yb(1)+1)./(yb(1)/ro_salt + 1/ro_water));
rho_ERD1= V_m*(rho_d1 - rho_ERD2) + rho_ERD2;
C_ERD1= -ro_salt*(rho_ERD1-ro_water)/(ro_water*(rho_ERD1-ro_salt));
% seawater after passing ERD 2
J_E1=(ya(1)*ya(2)+ya(2));
J_ERD1 = J_E1*(1-mix_M1);
J_wE1 = J_E1/(C_ERD2+1);
J_wERD1= J_ERD1/(C_ERD1+1);
% final calculations
swroC_in = (C_ERD2*J_wE1*mix_M1 + C_ERD1*J_wERD1)/(J_wE1*mix_M1+J_wERD1);
proC_in = (C_ERD2*J_wE1*(1-mix_M1) + yb(1)*yb(2)*(1-eta_ERD) - C_ERD1*J_wERD1 )/(J_wE1*(1-mix_M1) + yb(2)*(1-eta_ERD) - J_wERD1);
PRO_brine = J_r*( (yb(1)*yb(2)+yb(2))*(1-eta_ERD) )/Q_r;
end

%% define vector for residual error
res =   [ % SWRO part:
          ya(1)- swroC_in       % MIXTURE OF ERD AND SEA WATER
          ya(3)- J_sf_0         % salt flux in fresh side at 0
          ya(4)- J_wf_0         % water flux in fresh side at L 
          ya(5)- Pd_0           % pressure draw side at 0
          yb(5)- Pd_L           % pressure draw side at L
          yb(6)- Pf_L           % pressure fresh side at L

          % PRO part:
          ya(7)-  proC_in       % INCOMMING CONZENTRATION FROM ERD TO PRO
          yb(8) - PRO_brine
          ya(9)-  Q_sf_0       	% salt flux in fresh side at 0
          ya(12)- pf_0      	% pressure of fresh side at 0
          yb(11)- pd_L          % pressure draw side at L 
          yb(12)- pf_L];    	% pressure of fresh side at L
	
%% counter current (just that PRO model wont break)  
if version(1) ==1
    res(7)=yb(7) - proC_in; 
    res(8)=yb(6+1)*yb(6+2)+yb(6+2) + PRO_brine;
end 
end