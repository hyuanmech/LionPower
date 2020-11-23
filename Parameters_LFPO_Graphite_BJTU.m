%% Battery parameters from [Jiang et al., PS, 2013]

%% Geometric Parameters
% Cell area
p.Area = 12/23;  % Area of cell [m^2] (estimated based on [Smith and Wang, PS, 2006])

% SOC limits of positive and negative electrodes
p.theta_max_p = 0.1;     % at 100% cell SOC, estimated based on [Jiang et al., PS, 2013]
p.theta_max_n = 0.85510;  % at 100% cell SOC from LIONSIMBA
p.theta_min_p = 0.9;     % at 0% cell SOC, estimated based on [Jiang et al., PS, 2013]
p.theta_min_n = 0.01429;  % at 0% cell SOC from LIONSIMBA

% Thickness of each layer
p.L_n = 59e-6;     % Thickness of negative electrode [m]
p.L_s = 20e-6;     % Thickness of separator [m]
p.L_p = 92e-6;     % Thickness of positive electrode [m]
p.L_ccn = 9e-6;   % Thickness of negative current collector [m]
p.L_ccp = 16e-6;   % Thickness of positive current collector [m]

% Particle Radius
p.R_s_n = 14.75e-6;        % Radius of solid particles in negative electrode [m]
p.R_s_p = 1.15e-6;         % Radius of solid particles in positive electrode [m]

% Volume fractions 
p.epsilon_s_n = 0.56;      % Volume fraction in solid for neg. electrode
p.epsilon_s_p = 0.435;     % Volume fraction in solid for pos. electrode

p.epsilon_e_s = 0.4;       % Volume fraction in electrolyte for separator
p.epsilon_e_p = 0.28;      % Volume fraction in electrolyte for pos. electrode
p.epsilon_e_n = 0.3;       % Volume fraction in electrolyte for neg. electrode

% Specific interfacial surface area
p.a_s_n = 3*p.epsilon_s_n / p.R_s_n;  % Negative electrode [1/m]
p.a_s_p = 3*p.epsilon_s_p / p.R_s_p;  % Positive electrode [1/m]

% % Mass densities
% p.rho_sn = 2223;    % Solid phase in negative electrode [kg/m^3]
% p.rho_sp = 1500;    % Solid phase in positive electrode [kg/m^3]
% p.rho_e =  1210;    % Electrolyte [kg/m^3]
% p.rho_ccn = 8954;   % Current collector in negative electrode
% p.rho_ccp = 2707;   % Current collector in positive electrode
% 
% % Compute cell mass [kg/m^2]
% p.m_n = p.L_n * (p.rho_e*p.epsilon_e_n + p.rho_sn*p.epsilon_s_n);
% p.m_s = p.L_s * (p.rho_e*p.epsilon_e_n);
% p.m_p = p.L_p * (p.rho_e*p.epsilon_e_p + p.rho_sp*p.epsilon_s_p);
% p.m_cc = p.rho_ccn*p.L_ccn + p.rho_ccp*p.L_ccp;
% 
% % Lumped density [kg/m^3]
% p.density = (p.m_n+p.m_s+p.m_p+p.m_cc)/(p.L_n+p.L_s+p.L_p+p.L_ccn+p.L_ccp);

%% Concentrations
p.c_s_n_max = 31370;       % Max concentration in anode, [mol/m^3]
p.c_s_p_max = 22806;       % Max concentration in cathode, [mol/m^3]
p.c_e = 1000;              % Initial Li concentratioin in electrolyte, [mol/m^3]


%% Kinetic Parameters
p.R = 8.314472;         % Gas constant, [J/mol-K]
p.alph = 0.5;           % Charge transfer coefficients
p.brug = 1.5;           % Bruggeman porosity
p.R_f_n = 0;         % Resistivity of SEI layer, [Ohms*m^2]
p.R_f_p = 0;         % Resistivity of SEI layer, [Ohms*m^2]
p.Faraday = 96487;      % Faraday's constant, [Coulumbs/mol]
p.t_plus = 0.363;

% Activation energy for reactions
p.E_k_n = 20e3;         % Activiation energy for negative electrode, [J/mol]
p.E_k_p = 30e3;         % Activiation energy for positive electrode, [J/mol]

%% Transport Parameters
% Activation energy for diffusion
p.E_D_n = 35e3;        % Activiation energy for negative electrode, [J/mol]
p.E_D_p = 35e3;        % Activiation energy for positive electrode, [J/mol]

%% Cut-off Voltages
p.volt_max = 3.65;
p.volt_min = 2.5;

%% Electrical Conductivity of Electrodes
p.sig_n = 100;      % Conductivity of solid in neg. electrode, [1/Ohms*m]
p.sig_p = 0.5;     % Conductivity of solid in pos. electrode, [1/Ohms*m]

%% Density [kg / m^3 ]
% Aluminium current collector
p.rho_al = 2700; % At room temp. From material datasheets and standard scientific tables
% Positive electrode
p.rho_p  = 2500;
% Separator
p.rho_s  = 1100; % We do not know exactly the separator material used. Common separator materials have a density of 1070 - 1200 kg/m^3  (1070, 1100, 1200 etc.). Using a median value
% Negative electrode
p.rho_n  = 2500;
% Copper current collector
p.rho_cu = 8940; % At room temp. From material datasheets and standard scientific tables
% LiPF6 electrolyte
p.rho_LiPF6 = 1290; % Table 1 from 'Characterization of Lithium-Ion Battery Thermal Abuse Behavior Using Experimental and Computational Analysis',doi: 10.1149/2.0751510jes, J. Electrochem. Soc. 2015 volume 162, issue 10, A2163-A2173, also supported by 'Thermal analysis of a cylindrical lithium-ion battery', Xiongwen Zhang, Electrochimica Acta,  56 (2011) 1246–1255, Table 3

%% Specific heat capacities [ J / (kg K) ]
p.Cpal   = 897; % Aluminium current collector
p.Cpp    = 700; % Positive Electrode
p.Cps    = 700; % Separator
p.Cpn    = 700; % Negative Electrode
p.Cpcu   = 385; % Copper current collector
p.CpLiPF6 = 134.1; % Electrolyte

%% Thermal conductivities [ W / (m K) ]
% Aluminium current collector
p.k_al = 237;  % From material datasheets and standard scientific tables
% Positive electrode
p.k_p  = 2.1;
% Separator
p.k_s  = 0.16;
% Negative Electrode
p.k_n  = 1.7;
% Copper current collector
p.k_cu = 401;  % From material datasheets and standard scientific tables

%% Current collector conductivities [S/m]
p.sig_al = 3.55e7;
p.sig_cu = 5.96e7;

%% Thermodynamic Parameters
p.C_p = 1108;   % Heat capacity, [J/kg-K]
p.h = 20;       % Heat transfer coefficient, [W/K-m^2]
p.T_amb = 10+273.15;  % Ambient temperature, [K]
p.k = 3.9;      % Thermal conductivity, [W/m/K] from Rao et al 2016 Applied Energy