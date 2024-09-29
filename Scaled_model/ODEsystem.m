function [dJ_dp] = ODEsystem(x, J_p, DATA)
    %%  ODEsystem(x, J_p, DATA)
    %
    %   Scaled ODE system
    %
    %   contains 12 differential equations
    %
    %   Input:
    %       input1        -   [x, J_p, DATA]
    %       input2        -   x

    %% Read data
[H,Z,swro_Z,ro_water,ro_salt,Mw,Ms,Rw,T0,eta,sigma,p_r,rho_r,C_r,swro_L,swro_alpha,swro_R,swro_KK,swro_x_r,swro_b1,swro_b2,J_r,swro_gamma,swro_gamma2,swro_W_r,L,alpha,R,KK,x_r,b1,b2,Q_r,gamma,gamma2,W_r,cE,pE,rho_E,J_sf_0,J_wf_0,Pd_0,Pd_L,Pf_L,Q_sf_0,pd_0,pf_0,pd_L,pf_L,HP_eff,LP_eff,T_eff,V_m,ERD_eff,ERD_fric,A_ERD,eta_ERD,mix_density,pw,pe,swro_beta_fix,beta_fix,mixer_ERD,version,fig,swro_KF,swro_KD,KF,KD] ...
    = DATA(1);

%% Define Vector for unknowns of the model
dJ_dp = zeros(12,1);

%% SWRO
% Reynolds number
swro_ReH_d = J_r*((J_p(2) .* J_p(1)) + J_p(2)) / 2 / eta;   % Reynolds number for draw side
swro_ReH_f = J_r*(J_p(3) + J_p(4)) / 2 / eta;               % Reynolds number for fresh side% Function F_mix(x)
% Function F_mix(x)
swro_fmix_d = 96/swro_ReH_d*(4.86 + 0.65*sqrt(swro_ReH_d));        % F_mix,d(x)
swro_fmix_f = 96/swro_ReH_f*(4.86 + 0.65*sqrt(swro_ReH_f));        % F_mix,f(x)
% Local densitys
swro_local_ro_d = (J_p(1) + 1) ./ (J_p(1) / ro_salt + 1 / ro_water)/rho_r;           % local density of draw side
swro_local_ro_f = (J_p(3) + J_p(4)) ./ (J_p(3) / ro_salt + J_p(4) / ro_water)/rho_r; % local density of fresh side
% Osmotic pressure
swro_p_osm_d = ro_water*Rw*T0*log(1 + 2*Mw*J_p(1)/Ms)/p_r;         % osmotic pressure in draw side
swro_p_osm_f = ro_water*Rw*T0*log(1 + 2*Mw*J_p(3)/Ms./J_p(4))/p_r; % osmotic pressure in feed side
% Salt permeability
    if version(4) == 0
        swro_beta = 0;
    elseif version(2) == 0
        swro_beta = swro_beta_fix / p_r/ swro_alpha;
    else
        swro_beta = (1 - swro_R)*((J_p(5) - J_p(6)) - sigma.*(swro_p_osm_d - swro_p_osm_f) )./swro_R;
    end
    % Water Permeate flux J_win(x)
    if version(4) == 0 % ideal
        J_cross =((J_p(5) - J_p(6)) - sigma.*(swro_p_osm_d - swro_p_osm_f));
    elseif version(7) == 0 % ICP
        J_cross = ((J_p(5) - J_p(6)) - sigma.*(swro_p_osm_d - swro_p_osm_f)) ./   (1 + p_r*swro_alpha*swro_KK*sigma.*(swro_p_osm_d - swro_p_osm_f));
    else % ICP+ECP
        J_cross = ((J_p(5)-J_p(6))-sigma.*(swro_p_osm_d - swro_p_osm_f) + swro_beta*(J_p(5)-J_p(6))*(swro_KK +1/swro_KD + 1/swro_KF)) ./   (1 + swro_alpha*swro_beta*(1/swro_KD + 1/swro_KF + swro_KK) + p_r*swro_alpha*sigma.*(-swro_p_osm_d/swro_KD -swro_p_osm_f*(swro_KK + 1/swro_KF) ));
    end
    % Salt Permeate J_sin(x)
    if version(7) == 1 % ICP+ECP
        J_sin = swro_beta*( (J_p(1) - J_p(3) / J_p(4))-J_cross*( J_p(1)./(J_p(1) + 1)./swro_KD +  J_p(3)./(J_p(3) + J_p(4))*(swro_KK+1/swro_KF)))/(1+swro_alpha*swro_beta*(swro_KK+1/swro_KF+1/swro_KD));
    else
        J_sin = swro_beta*(J_p(1) - J_p(3) / J_p(4));
    end
    % Hydraulic diameter
    swro_DH_rect = 2*(swro_b2*swro_b1)/(swro_b2+swro_b1);

%% PRO
% Reynolds number
if version(1)==0
ReH_d = (Q_r) *((J_p(2+6).*J_p(1+6)) + J_p(2+6))/2/eta;   % Reynolds number for draw side
ReH_f = (Q_r) *(J_p(3+6) + J_p(4+6))/2/eta ;              % Reynolds number for fresh side
else
ReH_d =(Q_r)*abs((J_p(2+6).*J_p(1+6)) + J_p(2+6))/2/eta;    % for draw side
ReH_f =(Q_r)*abs(J_p(3+6) + J_p(4+6))/2/eta;              % for fresh side
end
% Function f_mix(x)
fmix_d = 96/ReH_d*(4.86 + 0.65*sqrt(ReH_d));              % f_mix,d(x)
fmix_f = 96/ReH_f*(4.86 + 0.65*sqrt(ReH_f));              % f_mix,f(x)
% Local densitys
local_ro_d = (J_p(1+6) + 1)./(ro_water*J_p(1+6)/ro_salt + 1);            	% local density of draw side
local_ro_f = (J_p(3+6) + J_p(4+6))./(ro_water*J_p(3+6)/ro_salt + J_p(4+6)); % local density of fresh side
% Osmotic pressure
p_osm_d = ro_water*Rw*T0*log(1 + 2*Mw*J_p(1+6)/Ms)/p_r;             % osmotic pressure in draw side
p_osm_f = ro_water*Rw*T0*log(1 + 2*Mw*J_p(3+6)/Ms./J_p(4+6))/p_r;   % osmotic pressure in fresh side
% B(from permeate salt flux)
beta = (1 - R)*((p_osm_d - p_osm_f) - (J_p(5+6) - J_p(6+6)))./R;
if version(3)==0; beta = beta_fix; end
if version(5)==0; beta = 0; end
% Difference in salt concentrations
del_c = J_p(1+6)./(J_p(1+6) + 1) - J_p(3+6)./(J_p(3+6) + J_p(4+6)); % delta c_salt
% Permeate flux Q_win(x)
Q_cross = ((p_osm_d - p_osm_f) - (J_p(5+6) - J_p(6+6)) - (J_p(5+6) - J_p(6+6))*beta*KK*alpha*p_r)./(1 + KK*alpha*p_r*(beta +p_osm_f));
if version(5)==0; Q_cross = ((p_osm_d - p_osm_f) - (J_p(5+6) - J_p(6+6))) ; end 
Q_sin = beta*del_c;
if version(8)==1 && version(4)>0
    Q_cross=((p_osm_d - p_osm_f) - (J_p(5+6) - J_p(6+6)) - (J_p(5+6) - J_p(6+6))*beta*(KK+1/KF+1/KD)*p_r*alpha)./(1 +alpha* beta*(KK+1/KF+1/KD) + alpha*p_r*(p_osm_f*(KK+1/KF) + p_osm_d/KD));
    Q_sin=beta*( del_c-Q_cross*( J_p(1+6)./(J_p(1+6) + 1)./KD +  J_p(3+6)./(J_p(3+6) + J_p(4+6))*(KK+1/KF)))/(1+alpha*beta*(KK+1/KF+1/KD));
end
% Hydraulic diameter
DH_rect = 2*(b2*b1)/(b2+b1);

%% derivatives of the SWRO flow model:
dJ_dp(1) =  +swro_gamma*(-J_sin+J_p(1)*J_cross)/J_p(2);   %C_d' 
dJ_dp(2) =  -swro_gamma*(J_cross);                        %J_wd'
dJ_dp(3) =  +swro_gamma*(J_sin);                          %J_sf'
dJ_dp(4) =  +swro_gamma*(J_cross);                        %J_wf'
dJ_dp(5) =  +swro_gamma2*(-swro_fmix_d*(J_p(2).*(J_p(1) + 1))*abs(J_p(2).*(J_p(1) + 1))/(2*swro_local_ro_d*swro_b1^2*swro_DH_rect)- swro_b1^-2*( dJ_dp(2).*(J_p(1).*J_p(2).*(1/ro_salt + 1/ro_water) + 2*J_p(2)/ro_water)));                                                                                                                                                                                                                        
dJ_dp(6) =  +swro_gamma2*(-swro_fmix_f*(J_p(3) + J_p(4))*abs(J_p(3) + J_p(4))/(2*swro_local_ro_f*swro_b1^2*swro_DH_rect) - swro_b1^-2*( dJ_dp(4).*(J_p(3)*(1/ro_salt + 1/ro_water) + 2*J_p(4)/ro_water)));

%% derivatives of the PRO flow model:
dJ_dp(1+6) = -gamma*(Q_sin+J_p(1+6)*Q_cross)/J_p(2+6); %c_d'
dJ_dp(2+6) = +gamma*(Q_cross);                         %Q_wd'
dJ_dp(3+6) = +gamma*(Q_sin);                           %Q_sf'
dJ_dp(4+6) = -gamma*(Q_cross);                         %Q_wf'
dJ_dp(5+6) = +gamma2*(-fmix_d*(J_p(2+6).*(J_p(1+6) + 1))*abs(J_p(2+6).*(J_p(1+6) + 1))/(2*local_ro_d*b1^2*DH_rect) - b1^-2*( dJ_dp(2+6).*(J_p(1+6).*J_p(2+6).*(1/ro_salt + 1/ro_water) + 2*J_p(2+6)/ro_water) ));
dJ_dp(6+6) = +gamma2*(-fmix_f*(J_p(3+6) + J_p(4+6))*abs(J_p(3+6) + J_p(4+6))/(2*local_ro_f*b1^2*DH_rect)- b1^-2*( dJ_dp(4+6).*(J_p(3+6)*(1/ro_salt + 1/ro_water) + 2*J_p(4+6)/ro_water)) );

end