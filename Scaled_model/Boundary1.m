function [ res ] = ERD_boundary1(ya, yb, DATA)
% Residual function for boundary conditions
% version(6) = 1 
% --> only SWRO (with ERD)

%% Read data
[H, Z, swro_Z, ro_water, ro_salt, Mw, Ms, Rw, T0, eta, sigma, p_r, rho_r, C_r, swro_L, swro_alpha, swro_KK, swro_x_r, swro_b1, swro_b2, J_r, swro_gamma, swro_gamma2, swro_W_r, L, alpha, KK, x_r, b1, b2, Q_r, gamma, gamma2, W_r, cE, pE, rho_E, J_sf_0, J_wf_0, Pd_0, Pd_L, Pf_L, Q_sf_0, pd_0, pf_0, pd_L, pf_L, HP_eff, LP_eff, T_eff, V_m, ERD_eff, ERD_fric, A_ERD, eta_ERD, mix_density, pw, pe, swro_beta_fix, beta_fix, mix_M1, mix_M3, version, fig, swro_KF, swro_KD, KF, KD]...
    = DATA(1);
 
%% Mixing and Leak
rho_d1= (yb(1) + 1)./(ro_water*yb(1)/ro_salt + 1);
rho_ERD= V_m*(rho_d1 - rho_E)+rho_E;
C_ERD= -ro_salt*(rho_ERD-1)/(rho_ERD*ro_water-ro_salt);
J_E=(ya(1)*ya(2)+ya(2))/2;
J_ERD = J_E;
J_wE = J_E/(cE+1);
J_wERD= J_ERD/(C_ERD+1);
swroC_in = (cE*J_wE + C_ERD*J_wERD)/(J_wE+J_wERD);
 
%% define vector for residual error
res =   [ % SWRO part:
          ya(1)- swroC_in       % seawater concentration
          ya(3)- J_sf_0         % salt flux in fresh side at 0
          ya(4)- J_wf_0         % water flux in fresh side at L 
          ya(5)- Pd_0           % pressure draw side at 0
          yb(5)- Pd_L           % pressure draw side at L
          yb(6)- Pf_L           % pressure fresh side at L

          % PRO part:
          ya(7)-  cE            % seawater concentration
          ya(9)-  Q_sf_0       	% salt flux in fresh side at 0
          ya(11)- pd_0      	% pressure draw side at 0
          ya(12)- pf_0      	% pressure of fresh side at 0
          yb(11)- pd_L          % pressure draw side at L 
          yb(12)- pf_L];    	% pressure of fresh side at L

%% (for comparison with senthil; set target production rate)
if Pd_L == 0
    res(5) = swro_Z*2*J_E*J_r/rho_d1/rho_r - 0.291;
end

%% counter current (just that PRO model wont break)  
if version(1) ==1; res(7)=yb(7) - cE; end 

end