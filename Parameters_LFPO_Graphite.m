%% Parameters from [Ye et al., JPS, 2012] and [Jiang et al., JPS, 2013]

%% Geometric Parameters
% Cell area
p.Area = 12/20.5;  % estimated

% SOC limits of positive and negative electrodes
p.theta_max_p = 0.06;     
p.theta_max_n = 0.9656;   
p.theta_min_p = 0.92;     
p.theta_min_n = 0.012;    

% Thickness of each layer
p.L_n = 59e-6;     % Thickness of negative electrode [m]
p.L_s = 20e-6;     % Thickness of separator [m]
p.L_p = 92e-6;     % Thickness of positive electrode [m]
p.L_ccn = 9e-6;    % Thickness of negative current collector [m]
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

%% Concentrations
p.c_s_n_max = 31370;       % Max concentration in anode, [mol/m^3]
p.c_s_p_max = 22806;       % Max concentration in cathode, [mol/m^3]
p.c_e = 1000;              % Initial Li concentratioin in electrolyte, [mol/m^3]


%% Kinetic Parameters
p.R = 8.314472;         % Gas constant, [J/mol-K]
p.alph = 0.5;           % Charge transfer coefficients
p.brug = 1.5;           % Bruggeman porosity
p.R_f_n = 0;            % Resistivity of SEI layer, [Ohms*m^2]
p.R_f_p = 0;            % Resistivity of SEI layer, [Ohms*m^2]
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
p.sig_p = 0.5;      % Conductivity of solid in pos. electrode, [1/Ohms*m]

%% Density 
p.rho_al = 2700;       % Aluminium current collector, [kg/m^3]
p.rho_p  = 1500;       % Positive electrode, [kg/m^3]
p.rho_s  = 1200;       % Separator, [kg/m^3]
p.rho_n  = 2500;       % Negative electrode, [kg/m^3]
p.rho_cu = 8900;       % Copper current collector, [kg/m^3]
p.rho_LiPF6 = 1290;    % LiPF6 electrolyte, [kg/m^3]

%% Specific heat capacities 
p.Cpal   = 897;     % Aluminium current collector, [J/(kg-K)]
p.Cpp    = 800;     % Positive Electrode, [J/(kg-K)]
p.Cps    = 800;     % Separator, [J/(kg-K)]
p.Cpn    = 800;     % Negative Electrode, [J/(kg-K)]
p.Cpcu   = 385;     % Copper current collector, [J/(kg-K)]
p.CpLiPF6 = 134.1;  % Electrolyte, [J/(kg-K)]

%% Thermal conductivities
p.k_al = 237;       % Aluminium current collector, [W/(m-K)]
p.k_p  = 1.48;      % Positive electrode, [W/(m-K)]
p.k_s  = 1;         % Separator, [W/(m-K)]
p.k_n  = 1.04;      % Negative Electrode, [W/(m-K)]
p.k_cu = 401;       % Copper current collector, [W/(m-K)]

%% Current collector conductivities 
p.sig_al = 3.55e7;   % [S/m]
p.sig_cu = 5.96e7;   % [S/m]

%% Thermodynamic Parameters
p.h = 20;       % Heat transfer coefficient, [W/K-m^2]
p.T_amb = 10+273.15;  % Ambient temperature, [K]
