function [ res ] = ERD_boundary3(ya, yb, DATA)
% Residual function for boundary conditions
% version(6) == 3
% --> SWRO-PRO hybrid system (with one ERD)

%% Read data
[H,Z,swro_Z,ro_water,ro_salt,Mw,Ms,Rw,T0,eta,sigma,p_r,rho_r,C_r,swro_L,swro_alpha,swro_R,swro_KK,swro_x_r,swro_b1,swro_b2,J_r,swro_gamma,swro_gamma2,swro_W_r,L,alpha,R,KK,x_r,b1,b2,Q_r,gamma,gamma2,W_r,cE,pE,rho_E,J_sf_0,J_wf_0,Pd_0,Pd_L,Pf_L,Q_sf_0,pd_0,pf_0,pd_L,pf_L,HP_eff,LP_eff,T_eff,V_m,ERD_eff,ERD_fric,A_ERD,eta_ERD,mix_density,pw,pe,swro_beta_fix,beta_fix,mixer_ERD,version,fig]...
     = DATA(1);

%% no Mixing i.e. V_m = 0
if V_m ==0 
swroC_in = cE;
c_in = yb(1);
if version(1)==0; J_wd_out = Z/swro_Z*Q_r/J_r*ya(6+2)*(c_in+1)/((1-eta_ERD)*(c_in+1)); end
if version(1)==1; J_wd_out = Z/swro_Z*abs(Q_r/J_r*yb(6+2))*(c_in+1)/((1-eta_ERD)*(c_in+1)); end
else
%% with Mixing i.e. V_m > 0
rho_d1= (yb(1) + 1)./(ro_water*yb(1)/ro_salt + 1);
rho_ERD= V_m*(rho_d1 - rho_E)+rho_E;
C_ERD= -ro_salt*(rho_ERD-1)/(rho_ERD*ro_water-ro_salt);
J_E=(ya(1)*ya(2)+ya(2))/2; J_ERD = J_E;
J_wE = J_E/(cE+1); J_wERD = J_ERD/(C_ERD+1);
swroC_in = (cE*J_wE + C_ERD*J_wERD)/(J_wE+J_wERD);
if version(1) ==0 % co-current
c_in = -(C_ERD * J_wERD * swro_Z * J_r + J_r * J_wE * yb(1) * swro_Z - cE * J_wE * swro_Z * J_r - J_r * J_wERD * yb(1) * swro_Z - ya(8) * yb(1) * Q_r * Z) / Z / Q_r / ya(8);
J_wd_out = (swro_Z*J_wERD+Z*Q_r/J_r*ya(2+6)-swro_Z*J_wE)/((1-eta_ERD)*swro_Z);
end
if version(1) ==1 % counter-current
c_in = -(C_ERD * J_wERD * swro_Z * J_r + J_r * J_wE * yb(1) * swro_Z - cE * J_wE * swro_Z * J_r - J_r * J_wERD * yb(1) * swro_Z - abs(yb(8)) * yb(1) * Q_r * Z) / Z / Q_r / abs(yb(8));
J_wd_out = (swro_Z*J_wERD+Z*Q_r/J_r*abs(yb(2+6))-swro_Z*J_wE)/((1-eta_ERD)*swro_Z);
end 
end

%% define vector for residual error
res =   [ % SWRO part:
          ya(1)- swroC_in       % MIXTURE OF ERD AND SEA WATER
          ya(3)- J_sf_0         % salt flux in fresh side at 0
          ya(4)- J_wf_0         % water flux in fresh side at L 
          ya(5)- Pd_0           % pressure draw side at 0
          %yb(5)- Pd_L          % pressure draw side at 0
          yb(2)- J_wd_out       % CORRECT MASS FLOWS BETWEEN SWRO AND PRO
          yb(6)- Pf_L           % pressure fresh side at L

          % PRO part:
          ya(7)-  c_in          % INCOMMING CONZENTRATION FROM ERD TO PRO
          ya(9)-  Q_sf_0       	% salt flux in fresh side at 0
          ya(11)- pd_0      	% pressure draw side at 0
          ya(12)- pf_0      	% pressure of fresh side at 0
          yb(11)- pd_L          % pressure draw side at L 
          yb(12)- pf_L];    	% pressure of fresh side at L
	
%% counter current    
if version(1) ==1; res(7)=yb(7) - c_in; end 

end