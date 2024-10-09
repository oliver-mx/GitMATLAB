function [H, Z, swro_Z, ro_water, ro_salt, Mw, Ms, Rw, T0, eta, sigma, p_r, rho_r, C_r, swro_L, swro_alpha, swro_KK, swro_x_r, swro_b1, swro_b2, J_r, swro_gamma, swro_gamma2, swro_W_r, L, alpha, KK, x_r, b1, b2, Q_r, gamma, gamma2, W_r, cE, pE, rho_E, J_sf_0, J_wf_0, Pd_0, Pd_L, Pf_L, Q_sf_0, pd_0, pf_0, pd_L, pf_L, HP_eff, LP_eff, T_eff, V_m, ERD_eff, ERD_fric, A_ERD, eta_ERD, mix_density, pw, pe, swro_beta_fix, beta_fix, mix_M1, mix_M3, version, fig, swro_KF, swro_KD, KF, KD] = data6(input1, input2)
    %%  data3(input)
    %
    %   Data for Test zwischen dem skalierten und unskalierten ODE system
    %
    %   only supports singleSWRO and SWRO+ERD
    %
    %   option_data = -.3
    %
    %   Input:
    %       input1        -   [FlowRate, Pd_0, Pd_L] (if FlowRate==0 --> Pd_L can be set)
    %       input2        -    mix_M1

    if nargin == 1
        input2=.5;
    end

    %% model versions
    version = zeros(1, 10);
    % version(1) = 0 if co-current PRO, 1 otherwise
    version(2) = input1(1);  % <-- unscaled RO feed massflow rate 
    version(3) = 0;  % <-- unscaled PRO feed massflow rate
    version(4) = 1;  % 0 = ideal SWRO
    version(5) = 0;  % 0 = ideal PRO
    % configuration:
    version(6) = 1;
    % version(6) = 0 --> only SWRO (no ERD)
    % version(6) = 1 --> only SWRO (with ERD)
    % version(6) = 2 --> only PRO
    % version(6) = 3 --> SWRO-PRO hybrid system (with one ERD)
    % version(6) = 4 --> SWRO-PRO hybrid system (with two ERDs) 
    version(7) = 1; % 1 = ICP and ECP for SWRO enabled
    version(8) = 0; % 1 = ICP and ECP for PRO enabled

    %% Membrane unit properties
    H = 1e-3;           % height of the membrane [m]
    ro_water = 1000;    % mass density of water  [kg/m^3] 
    ro_salt = 2165;     % mass density of salt [kg/m^3] 
    Mw = 18;            % molecular weight of water [kg/kmol]
    Ms = 58.44;         % molecular weight of salt  [kg/kmol]
    Rw = 462;           % gas constant of water  [J/(kg K)] 
    T0 = 298;           % temperature [K] (i.e. 24.85°C) 
    eta = 1.3e-3;       % seawater viscosity [kg/(m s)]
    p_r = 1e5;          % pressure [Pa]=[kg/ms^2]
    rho_r = 1e3;        % density [kg/m^3]
    C_r = 1;            % salt concentration [%]

    %% SWRO
    swro_Z = 7.7210;        % width of SWRO membrane [m]
    swro_L = 7 * 0.9626;    % length of SWRO membrane [m]
    swro_alpha = 5.0815e-9; % SWRO water permeability coefficient [s/m]
    swro_KK = 1e2;          % 
    swro_KD = 1/swro_KK;    % SWRO ECP draw side mass transfer coefficient
    swro_KF = 1/swro_KK;    % SWRO ECP fresh side mass transfer coefficient
    swro_x_r = swro_L;      % x_r = swro_L^2 since x = linspace(0,1,n) (if x = linspace(0,swro_L,n) then x_r=swro_L;)     
    swro_b1 = H / swro_x_r;                           % H/swro_L ratio
    swro_b2 = swro_Z / swro_x_r;                      % Z/swro_L ratio
    J_r = sqrt(H^3 / swro_x_r * p_r * rho_r);         % flux [kg/s^2]
    swro_gamma = swro_x_r * p_r * swro_alpha / J_r;   % SWRO scaling factor - mass balance
    swro_gamma2 = J_r^2 / (swro_x_r^2 * p_r * rho_r); % SWRO scaling factor - momentum balance
    swro_W_r = J_r * p_r;                             % net work [W/m^2]
    sigma = 0.999 ;                                   % Reflection coefficient
    swro_beta_fix = 4.43e-4;                          % value for fixed SWRO beta [kg/sm^2]

    %% PRO
    Z = 7.7210;         % width of the PRO membrane [m]
    L = 7 * 0.9626;     % length of the PRO membrane [m]      
    alpha = 5.47e-9;    % water permeability coefficient [s/m]
    KK = 7.13e2;        % mass transfer coefficient [sm^2/kg]
    KD = -1 / KK;       % PRO ECP draw side mass transfer coefficient [sm^2/kg]
    KF = -1 / KK;       % PRO ECP fresh side mass transfer coefficient [sm^2/kg]
    x_r = L;            % x_r = L^2 since x = linspace(0,1,n) (if x = linspace(0,L,n) then x_r=L;)              
    b1 = H / x_r;       % H/L ratio
    b2 = Z / x_r;       % Z/L ratio
    Q_r = sqrt(H^3 / x_r * p_r * rho_r);       % flux [kg/s^2]
    gamma = x_r * p_r * alpha / Q_r ;          % PRO scaling factor - mass balance
    gamma2 = Q_r^2 / (x_r^2 * p_r * rho_r);    % PRO scaling factor - momentum balance
    W_r = Q_r * p_r;                           % net work [W/m^2]
    beta_fix = 1.71e-4;                        % value for fixed PRO beta [kg/sm^2]

    %% Sea Water
    cE =(32/(1-32/2165)/1000) / C_r;                  % salt concentration in seawater (32/(1-32/2165)/1000)
    pE = 1e5 / p_r;                                   % external pressure
    rho_E = (cE + 1) / (cE / ro_salt + 1 / ro_water); % density of incoming seawater

    %% SWRO operating conditions
    J_sf_0 = 0;              % salt flux in fresh side at 0 
    J_wf_0 = 0;              % water flux in fresh side at L               
    Pf_L = pE;               % pressure fresh side at L

    Pd_0 =  input1(2); % pressure draw side at 0
    Pd_L =  input1(3); % pressure draw side at L (not needed in the hybrid system)

    %% PRO operation conditions 
    Q_sf_0 = 0;              % salt flux in fresh side at 0
    pf_L = pE;               % pressure of fresh side at L

    pd_0 = 14;       % pressure draw side at 0
    pd_L = 13.8;     % pressure of fresh side at 0
    pf_0 = 1.01;     % pressure draw side at L 

    %% ERD/Turbine/Pump parameters
    T_eff  = .95;               % turbine efficiency
    HP_eff = .9;            	% high pressure pump efficiency
    LP_eff = .95;               % low pressure pump efficiency
    V_m = 0.052;                % Volumetric mixing
    ERD_eff = .96;              % ERD unit pressure efficiency
    ERD_fric = 5e-04;           % ERD friction coefficient
    A_ERD = H * swro_Z;         % cross-sectional area of ERD inflows/outflows
    eta_ERD = 0.01;             % leak of high-pressure brine
    mix_density = 997 / rho_r;  % density of mixture in ERD
    pw = 2.0;                   % water price [$/m^3]
    pe = 0.3;                   % electricity price [$/kWh]
    mix_M1 = input2;   % spliting rate at M1 (if 0 --> all to ERD1)  
    mix_M3 = 1;        % spliting rate at M3 (if 0 --> all to ERD2)

    %% display figures
    fig = [0, 0, 0, 0]; % f(i) = 1 --> figure i will be displayed

    %% automatic changes:
    % co-current
    if pd_0 > pd_L
        version(1) = 0;
    else
        version(1) = 1;
    end            
    % no SWRO needed (--> use trivial data)
    if version(6) == 2
        version(4) = 0;
        Pd_0 = 40;
        Pd_L = 35;
    end
    
end