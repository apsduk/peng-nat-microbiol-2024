function [dZ_by_dt, J0_grow, J0_upt_gluc, J0_leak_Y, J0_upt_Y] = ode3v2strains(T, Z, theta_biomass, theta_yA, theta_yB, theta_yC, omega)

%% ===== VARIABLES ========================================================
Z(Z < 0) = 0;
Gluc = Z( 1); % glucose concentration in the culture vessel
yA   = Z( 2); % A+/B-/C- total population density of strain A
yB   = Z( 3); % A-/B+/C- total population density of strain B
yC   = Z( 4); % A-/B-/C+ total population density of strain C
yAa  = Z( 5); % A+/B-/C- total active population density of strain A
yBa  = Z( 6); % A-/B+/C- total active population density of strain B
yCa  = Z( 7); % A-/B-/C+ total active population density of strain C
A    = Z( 8); % A concentration in the culture vessel
B    = Z( 9); % B concentration in the culture vessel
C    = Z(10); % C concentration in the culture vessel

%% ===== PARAMETERS =======================================================
% ----- Biomass parameters ------------------------------------------------
gamma_gluc = theta_biomass(1); % Biomass yield of S. cerevisiae grown on glucose
gamma_A    = theta_biomass(2); % Biomass yield of S. cerevisiae grown on A
gamma_B    = theta_biomass(3); % Biomass yield of S. cerevisiae grown on B
gamma_C    = theta_biomass(4); % Biomass yield of S. cerevisiae grown on C

% ----- Strain A parameters (A+/B-/C-) ------------------------------------
% ----- (e.g. lys production, his uptake, ade uptake) ---------------------
Vmax_gluc_yA = theta_yA(1); % Maximum glucose uptake rate
Km_gluc_yA   = theta_yA(2); % Michaelis constant for glucose uptake
eta_yA       = theta_yA(3); % Mortality rate constant of strain A
Vmax_yA_B    = theta_yA(4); % Maximum B uptake rate by strain A
Km_yA_B      = theta_yA(5); % Michaelis constant for B uptake by strain A
Vmax_yA_C    = theta_yA(6); % Maximum C uptake rate by strain A
Km_yA_C      = theta_yA(7); % Michaelis constant for C uptake by strain A
delta_yA_A   = theta_yA(8); % Number of A molecules produced per molecule of glucose is consumed
phi_yA_A     = omega( 1);   % Fraction of glucose uptake leaked in form of A production by strain A

% ----- Strain B parameters (A-/B+/C-) ------------------------------------
% ----- (e.g. his production, lys uptake, ade uptake) ---------------------
Vmax_gluc_yB = theta_yB(1); % Maximum glucose uptake rate
Km_gluc_yB   = theta_yB(2); % Michaelis constant for glucose uptake
eta_yB       = theta_yB(3); % Mortality rate constant of strain B
Vmax_yB_A    = theta_yB(4); % Maximum A uptake rate by strain B
Km_yB_A      = theta_yB(5); % Michaelis constant for A uptake by strain B
Vmax_yB_C    = theta_yB(6); % Maximum C uptake rate by strain B
Km_yB_C      = theta_yB(7); % Michaelis constant for C uptake by strain B
delta_yB_B   = theta_yB(8); % Number of B molecules produced per molecule of glucose is consumed
phi_yB_B     = omega( 2);   % Fraction of glucose uptake leaked in form of B production by the strain B

% ----- Strain C parameters (A-/B-/C+) ------------------------------------
% ----- (e.g. ade production, lys uptake, his uptake) ---------------------
Vmax_gluc_yC = theta_yC(1); % Maximum glucose uptake rate
Km_gluc_yC   = theta_yC(2); % Michaelis constant for glucose uptake
eta_yC       = theta_yC(3); % Mortality rate constant of strain C
Vmax_yC_A    = theta_yC(4); % Maximum A uptake rate by strain C
Km_yC_A      = theta_yC(5); % Michaelis constant for A uptake by strain C
Vmax_yC_B    = theta_yC(6); % Maximum B uptake rate by strain C
Km_yC_B      = theta_yC(7); % Michaelis constant for B uptake by strain C
delta_yC_C   = theta_yC(8); % Number of C molecules produced per molecule of glucose is consumed
phi_yC_C     = omega( 3);   % Fraction of glucose uptake leaked in form of C production by the strain C

%% ===== CALCULATE RATES ==================================================
% ----- Glucose uptake rate -----------------------------------------------
J_upt_yA_gluc = (Vmax_gluc_yA*Gluc)/(Km_gluc_yA + Gluc);
J_upt_yB_gluc = (Vmax_gluc_yB*Gluc)/(Km_gluc_yB + Gluc);
J_upt_yC_gluc = (Vmax_gluc_yC*Gluc)/(Km_gluc_yC + Gluc);

% ----- Leakage flux ------------------------------------------------------
J_leak_yA_A = J_upt_yA_gluc*delta_yA_A*phi_yA_A;
J_leak_yB_B = J_upt_yB_gluc*delta_yB_B*phi_yB_B;
J_leak_yC_C = J_upt_yC_gluc*delta_yC_C*phi_yC_C;

% ----- Amino acid uptake rate --------------------------------------------
J_upt_yA_B = (Vmax_yA_B*B)/(Km_yA_B + B);
J_upt_yA_C = (Vmax_yA_C*C)/(Km_yA_C + C);
J_upt_yB_A = (Vmax_yB_A*A)/(Km_yB_A + A);
J_upt_yB_C = (Vmax_yB_C*C)/(Km_yB_C + C);
J_upt_yC_A = (Vmax_yC_A*A)/(Km_yC_A + A);
J_upt_yC_B = (Vmax_yC_B*B)/(Km_yC_B + B);

% ----- Specific growth rate ----------------------------------------------
J_grow_yA_G = gamma_gluc*(1 - phi_yA_A)*J_upt_yA_gluc;
J_grow_yA_B = gamma_B*J_upt_yA_B;
J_grow_yA_C = gamma_C*J_upt_yA_C;

J_grow_yB_A = gamma_A*J_upt_yB_A;
J_grow_yB_G = gamma_gluc*(1 - phi_yB_B)*J_upt_yB_gluc;
J_grow_yB_C = gamma_C*J_upt_yB_C;

J_grow_yC_A = gamma_A*J_upt_yC_A;
J_grow_yC_B = gamma_B*J_upt_yC_B;
J_grow_yC_G = gamma_gluc*(1 - phi_yC_C)*J_upt_yC_gluc;

% ----- Calculate specific growth rates -----------------------------------
J_grow_yA = min([J_grow_yA_G, J_grow_yA_B, J_grow_yA_C]);
J_grow_yB = min([J_grow_yB_A, J_grow_yB_G, J_grow_yB_C]);
J_grow_yC = min([J_grow_yC_A, J_grow_yC_B, J_grow_yC_G]);

%% === DIFFERENTIAL EQUATIONS =============================================
% --- Glucose dynamics ----------------------------------------------------
dGluc_by_dt = - J_upt_yA_gluc*yAa - J_upt_yB_gluc*yBa - J_upt_yC_gluc*yCa;

% --- Dynamics of the total population density of strain A ----------------
dyA_by_dt = J_grow_yA*yA;

% --- Dynamics of the total population density of strain B ----------------
dyB_by_dt = J_grow_yB*yB;

% --- Dynamics of the total population density of strain C ----------------
dyC_by_dt = J_grow_yC*yC;

% --- Dynamics of the active population density of strain A ---------------
dyAa_by_dt = (J_grow_yA - eta_yA)*yAa;

% --- Dynamics of the active population density of strain B ---------------
dyBa_by_dt = (J_grow_yB - eta_yB)*yBa;

% --- Dynamics of the active population density of strain C ---------------
dyCa_by_dt = (J_grow_yC - eta_yC)*yCa;

% --- Amino acid A dynamics -----------------------------------------------
dA_by_dt = J_leak_yA_A*yAa - J_upt_yB_A*yBa - J_upt_yC_A*yCa;

% --- Amino acid B dynamics -----------------------------------------------
dB_by_dt = J_leak_yB_B*yBa - J_upt_yA_B*yAa - J_upt_yC_B*yCa;

% --- Amino acid C dynamics -----------------------------------------------
dC_by_dt = J_leak_yC_C*yCa - J_upt_yA_C*yAa - J_upt_yB_C*yBa;

%% === RETURN =============================================================
dZ_by_dt = [dGluc_by_dt; dyA_by_dt; dyB_by_dt; dyC_by_dt; dyAa_by_dt; dyBa_by_dt; dyCa_by_dt; dA_by_dt; dB_by_dt; dC_by_dt];

% --- Return fluxes -------------------------------------------------------
J0_grow      = [J_grow_yA*yAa                  , J_grow_yB*yBa                  , J_grow_yC*yCa];
J0_upt_gluc  = [J_upt_yA_gluc*yAa              , J_upt_yB_gluc*yBa              , J_upt_yC_gluc*yCa];
J0_leak_Y    = [J_leak_yA_A*yAa                , J_leak_yB_B*yBa                , J_leak_yC_C*yCa];
J0_upt_Y     = [J_upt_yB_A*yBa + J_upt_yC_A*yCa, J_upt_yA_B*yAa + J_upt_yC_B*yCa, J_upt_yA_C*yAa + J_upt_yB_C*yBa];

end