function [dZ_by_dt, J0_grow, J0_upt_gluc, J0_leak_yX_X, J0_upt_yX_Y, J_grow_yA_gluc, J_grow_yA_B] = ode1strain(T, Z, theta_biomass, theta_yA, omega)

%% === VARIABLES ==========================================================
Z(Z < 0) = 0;
Gluc = Z(1); % glucose concentration in the culture vessel
yA   = Z(2); % A+/B- total population density of strain A
yAa  = Z(3); % A+/B- total active population density of strain A
A    = Z(4); % A concentration in the culture vessel
B    = Z(5); % B concentration in the culture vessel

%% === PARAMETERS =========================================================
% --- Biomass parameters --------------------------------------------------
gamma_gluc = theta_biomass(1); % Biomass yield of S. cerevisiae grown on glucose
gamma_B    = theta_biomass(2); % Biomass yield of S. cerevisiae grown on B

% --- Strain A parameters -------------------------------------------------
Vmax_gluc_yA = theta_yA(1); % Maximum glucose uptake rate
Km_gluc_yA   = theta_yA(2); % Michaelis constant for glucose uptake
eta_yA       = theta_yA(3); % Mortality rate constant of strain A
Vmax_yA_B    = theta_yA(4); % Maximum B uptake rate by the strain A
Km_yA_B      = theta_yA(5); % Michaelis constant for B uptake
delta_yA_A   = theta_yA(6); % Number of A molecules produced per molecule of glucose is consumed
phi_yA_A     = omega( 1);   % Fraction of glucose uptake leaked in form of A production by strain A

%% === CALCULATE RATES ====================================================
% --- Glucose uptake rate -------------------------------------------------
J_upt_yA_gluc = (Vmax_gluc_yA*Gluc)/(Km_gluc_yA + Gluc);

% --- Leakage flux --------------------------------------------------------
J_leak_yA_A = J_upt_yA_gluc*delta_yA_A*phi_yA_A; % A production from strain A

% --- Amino acid uptake rate ----------------------------------------------
J_upt_yA_B = (Vmax_yA_B*B)/(Km_yA_B + B);

% --- Specific growth rate ------------------------------------------------
% J_grow_yA_gluc = gamma_gluc*(1 - phi_yA_A)*J_upt_yA_gluc;
% J_grow_yA_B = gamma_B*J_upt_yA_B;
J_grow_yA = min([gamma_B*J_upt_yA_B, gamma_gluc*(1 - phi_yA_A)*J_upt_yA_gluc]);

%% === DIFFERENTIAL EQUATIONS =============================================
% --- Glucose dynamics ----------------------------------------------------
dGluc_by_dt = - J_upt_yA_gluc*yAa;

% --- Dynamics of the total population density of strain A ----------------
dyA_by_dt = J_grow_yA*yA;

% --- Dynamics of the active population density of strain A ---------------
dyAa_by_dt = (J_grow_yA - eta_yA)*yAa;

% --- Amino acid A dynamics -----------------------------------------------
dA_by_dt = J_leak_yA_A*yAa;

% --- Amino acid B dynamics -----------------------------------------------
dB_by_dt = - J_upt_yA_B*yAa;

%% === RETURN =============================================================
dZ_by_dt = [dGluc_by_dt; dyA_by_dt; dyAa_by_dt; dA_by_dt; dB_by_dt];

% --- Return fluxes -------------------------------------------------------
J0_grow      = [J_grow_yA        ];
J0_upt_gluc  = [J_upt_yA_gluc*yAa];
J0_leak_yX_X = [J_leak_yA_A*yAa  ];
J0_upt_yX_Y  = [J_upt_yA_B*yAa   ];    

end