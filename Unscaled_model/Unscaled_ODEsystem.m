function [dJ_dp] = Unscaled_ODEsystem(x, J_p, DATA)
    %%  Unscaled_ODEsystem(x, J_p, DATA)
    %
    %   Unscaled ODE system
    %
    %   contains ONLY 6 differential equations
    %
    %   Input:
    %       input1        -   [x, J_p, DATA]
    %       input2        -   x

    %% Read data
    [H, Z, swro_Z, ro_water, ro_salt, Mw, Ms, Rw, T0, eta, sigma, p_r, rho_r, C_r, swro_L, swro_alpha, swro_R, swro_KK, swro_x_r, swro_b1, swro_b2, J_r, swro_gamma, swro_gamma2, swro_W_r, L, alpha, R, KK, x_r, b1, b2, Q_r, gamma, gamma2, W_r, cE, pE, rho_E, J_sf_0, J_wf_0, Pd_0, Pd_L, Pf_L, Q_sf_0, pd_0, pf_0, pd_L, pf_L, HP_eff, LP_eff, T_eff, V_m, ERD_eff, ERD_fric, A_ERD, eta_ERD, mix_density, pw, pe, swro_beta_fix, beta_fix, mixer_ERD, version, fig, swro_KF, swro_KD, KF, KD] ...
        = DATA(1);

    %% Define Vector for unknowns of the model
    dJ_dp = zeros(6, 1);

    %% SWRO
    % Reynolds number
    swro_ReH_d = ((J_p(2) .* J_p(1)) + J_p(2)) / 2 / eta;   % Reynolds number for draw side
    swro_ReH_f = (J_p(3) + J_p(4)) / 2 / eta;               % Reynolds number for fresh side
    % Function F_mix(x)
    swro_fmix_d = 96 / swro_ReH_d * (4.86 + 0.65 * sqrt(swro_ReH_d));        % F_mix,d(x)
    swro_fmix_f = 96 / swro_ReH_f * (4.86 + 0.65 * sqrt(swro_ReH_f));        % F_mix,f(x)
    % Local densitys
    swro_local_ro_d = (J_p(1) + 1) ./ (J_p(1) / ro_salt + 1 / ro_water);           % local density of draw side
    swro_local_ro_f = (J_p(3) + J_p(4)) ./ (J_p(3) / ro_salt + J_p(4) / ro_water); % local density of fresh side
    % Osmotic pressure
    swro_p_osm_d = ro_water * Rw * T0 * log(1 + 2 * Mw * J_p(1) / Ms);           % osmotic pressure in draw side
    swro_p_osm_f = ro_water * Rw * T0 * log(1 + 2 * Mw * J_p(3) / Ms ./ J_p(4)); % osmotic pressure in feed side
    % Salt permeability
    if version(4) == 0
        swro_beta = 0;
    elseif version(2) == 0
        swro_beta = swro_beta_fix;
    else
        swro_beta = (1 - swro_R) * ((J_p(5) - J_p(6)) - sigma .* (swro_p_osm_d - swro_p_osm_f)) ./ swro_R;
    end
    % Water Permeate flux J_win(x)
    if version(4) == 0 % ideal
        J_cross = swro_alpha * ((J_p(5) - J_p(6)) - sigma .* (swro_p_osm_d - swro_p_osm_f));
    elseif version(7) == 0 % ICP
        J_cross = swro_alpha * ((J_p(5) - J_p(6)) - sigma.*(swro_p_osm_d - swro_p_osm_f)) ./   (1 + swro_alpha*swro_KK*sigma.*(swro_p_osm_d - swro_p_osm_f));
    else % ICP+ECP
        J_cross = swro_alpha * ((J_p(5) - J_p(6)) * (1 + swro_beta * (1 / swro_KD + swro_KK + 1 / swro_KF)) - sigma * (swro_p_osm_d - swro_p_osm_f)) / (1 + swro_beta * (1 / swro_KD + swro_KK + 1 / swro_KF) - swro_alpha * sigma * (swro_p_osm_d / swro_KD + swro_p_osm_f * swro_KK + swro_p_osm_f / swro_KF));
    end
    % Salt Permeate J_sin(x)
    if version(7) == 1 % ICP+ECP
        J_sin = swro_beta * ((J_p(1) - J_p(3) / J_p(4)) - J_p(1) * J_cross / swro_KD - (J_p(3) / J_p(4)) * J_cross * (swro_KK + 1 / swro_KF)) / (1 + swro_beta * (1 / swro_KD + swro_KK + 1 / swro_KF));
    else
        J_sin = swro_beta * (J_p(1) - J_p(3) / J_p(4));
    end
    % Hydraulic diameter
    swro_DH_rect = 2 * (H * swro_Z) / (H + swro_Z);
    
    %% derivatives of the SWRO flow model:
    dJ_dp(1) = (-J_sin + J_p(1) * J_cross) / J_p(2);  % C_d' 
    dJ_dp(2) = -(J_cross);                            % J_wd'
    dJ_dp(3) = +(J_sin);                              % J_sf'
    dJ_dp(4) = +(J_cross);                            % J_wf'
    dJ_dp(5) = -swro_fmix_d * (J_p(2) .* (J_p(1) + 1)) * abs(J_p(2) .* (J_p(1) + 1)) / (2 * swro_local_ro_d * H^2 * swro_DH_rect) - H^-2 * ( ...                                                                                                                                                                                  
        2 * J_p(2) * (J_p(1) + 1) * (J_p(1) / ro_salt + 1 / ro_water) * dJ_dp(2) + J_p(2)^2 * dJ_dp(1) * (J_p(1) / ro_salt + 1 / ro_water) + J_p(2)^2 * (J_p(1) + 1) * dJ_dp(1) / ro_salt);
    dJ_dp(6) = -swro_fmix_f * (J_p(3) + J_p(4)) * abs(J_p(3) + J_p(4)) / (2 * swro_local_ro_f * H^2 * swro_DH_rect) - H^-2 * (...
        (dJ_dp(3) + dJ_dp(4)) * (J_p(3) / ro_salt + J_p(4) / ro_water) + (J_p(3) + J_p(4)) * (dJ_dp(3) / ro_salt + dJ_dp(4) / ro_water));

end