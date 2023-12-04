function [dZ_by_dt, J0_grow, J0_upt_gluc, J0_leak_yX_X, J0_upt_yX_Y] = ode2REUPstrains(T, Z, theta_biomass, theta_yA, theta_yB, omega)

%% === VARIABLES ==========================================================
Z(Z < 0) = 0;
Gluc = Z(1); % glucose concentration in the culture vessel
yA   = Z(2); % A+/B- total population density of strain A
yB   = Z(3); % A-/B+ total population density of strain B
yAa  = Z(4); % A+/B- total active population density of strain A
yBa  = Z(5); % A-/B+ total active population density of strain B
A    = Z(6); % A conconcentration in the culture vessel
B    = Z(7); % B concentration in the culture vessel

%% === PARAMETERS =========================================================
% --- Biomass parameters --------------------------------------------------
gamma_gluc = theta_biomass(1); % Biomass yield of S. cerevisiae grown on glucose
gamma_A    = theta_biomass(2); % Biomass yield of S. cerevisiae grown on A
gamma_B    = theta_biomass(3); % Biomass yield of S. cerevisiae grown on B

% --- Strain A parameters -------------------------------------------------
Vmax_gluc_yA = theta_yA(1); % Maximum glucose uptake rate
Km_gluc_yA   = theta_yA(2); % Michaelis constant for glucose uptake
eta_yA       = theta_yA(3); % Mortality rate constant of strain A
Vmax_yA_A    = theta_yA(4); % Maximum A uptake rate by the strain A
Km_yA_A      = theta_yA(5); % Michaelis constant for A uptake
Vmax_yA_B    = theta_yA(6); % Maximum B uptake rate by the strain A
Km_yA_B      = theta_yA(7); % Michaelis constant for B uptake
delta_yA_A   = theta_yA(8); % Number of A molecules produced per molecule of glucose is consumed
phi_yA_A     = omega( 1);   % Fraction of glucose uptake leaked in form of A production by strain A

% --- Strain B parameters -------------------------------------------------
Vmax_gluc_yB = theta_yB(1); % Maximum glucose uptake rate
Km_gluc_yB   = theta_yB(2); % Michaelis constant for glucose uptake
eta_yB       = theta_yB(3); % Mortality rate constant of strain B
Vmax_yB_A    = theta_yB(4); % Maximum A uptake rate by the strain B
Km_yB_A      = theta_yB(5); % Michaelis constant for A uptake
Vmax_yB_B    = theta_yB(6); % Maximum A uptake rate by the strain B
Km_yB_B      = theta_yB(7); % Michaelis constant for A uptake
delta_yB_B   = theta_yB(8); % Number of B molecules produced per molecule of glucose is consumed
phi_yB_B     = omega( 2);   % Fraction of glucose uptake leaked in form of B production by the strain B

%% === CALCULATE RATES ====================================================
% --- Glucose uptake rate -------------------------------------------------
J_upt_yA_gluc = (Vmax_gluc_yA*Gluc)/(Km_gluc_yA + Gluc);
J_upt_yB_gluc = (Vmax_gluc_yB*Gluc)/(Km_gluc_yB + Gluc);

% --- Leakage flux --------------------------------------------------------
J_leak_yA_A = J_upt_yA_gluc*delta_yA_A*phi_yA_A; % A production from strain A
J_leak_yB_B = J_upt_yB_gluc*delta_yB_B*phi_yB_B; % B production from strain B

% --- Amino acid uptake rate ----------------------------------------------
J_upt_yA_A = (Vmax_yA_A*A)/(Km_yA_A + A);
J_upt_yA_B = (Vmax_yA_B*B)/(Km_yA_B + B);
J_upt_yB_A = (Vmax_yB_A*A)/(Km_yB_A + A);
J_upt_yB_B = (Vmax_yB_B*B)/(Km_yB_B + B);

% --- Specific growth rate ------------------------------------------------
J_grow_yA = min([gamma_A*J_upt_yA_A, gamma_B*J_upt_yA_B, gamma_gluc*(1 - phi_yA_A)*J_upt_yA_gluc]);
J_grow_yB = min([gamma_A*J_upt_yB_A, gamma_B*J_upt_yB_B, gamma_gluc*(1 - phi_yB_B)*J_upt_yB_gluc]);

%% === DIFFERENTIAL EQUATIONS =============================================
% --- Glucose dynamics ----------------------------------------------------
dGluc_by_dt = - J_upt_yA_gluc*yAa - J_upt_yB_gluc*yBa;

% --- Dynamics of the total population density of strain A ----------------
dyA_by_dt = J_grow_yA*yA;

% --- Dynamics of the total population density of strain B ----------------
dyB_by_dt = J_grow_yB*yB;

% --- Dynamics of the active population density of strain A ---------------
dyAa_by_dt = (J_grow_yA - eta_yA)*yAa;

% --- Dynamics of the active population density of strain B ---------------
dyBa_by_dt = (J_grow_yB - eta_yB)*yBa;

% --- Amino acid A dynamics -----------------------------------------------
dA = J_leak_yA_A*yAa - J_upt_yB_A*yBa - J_upt_yA_A*yAa;

% --- Amino acid B dynamics -----------------------------------------------
dB = J_leak_yB_B*yBa - J_upt_yA_B*yAa - J_upt_yB_B*yBa;

%% === RETURN =============================================================
dZ_by_dt = [dGluc_by_dt; dyA_by_dt; dyB_by_dt; dyAa_by_dt; dyBa_by_dt; dA; dB];

% --- Return fluxes -------------------------------------------------------
J0_grow      = [J_grow_yA*yAa                  , J_grow_yB*yBa];
J0_upt_gluc  = [J_upt_yA_gluc*yAa              , J_upt_yB_gluc*yBa];
J0_leak_yX_X = [J_leak_yA_A*yAa                , J_leak_yB_B*yBa];
J0_upt_yX_Y  = [J_upt_yB_A*yBa + J_upt_yA_A*yAa, J_upt_yA_B*yAa + J_upt_yB_B*yBa];    

end