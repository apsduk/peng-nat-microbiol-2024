function [performancevector, performancenames, T, yNa, Z] = efast2MEv1strains(X, kvector, kidx, tmax, odeoptions)

%% ===== SUBSTITUTE PARAMETERS DEFINED BY kidx ============================
kvector(kidx) = X;

%% ===== GET SPECIFIC VECTORS OUT OF INPUT ================================
% ----- Initial conditions ------------------------------------------------
dN0 = kvector(1); % total initial population
r0A = kvector(2); % A ratio
r0B = kvector(3); % B ratio
G0  = kvector(4); % glucose conc.
aA0 = kvector(5); % amino acid A
aB0 = kvector(6); % amino acid B

% ------ Other parameters -------------------------------------------------
theta_biomass = kvector( 7: 9);
theta_yA      = kvector(10:15);
theta_yB      = kvector(16:21);
omega         = kvector(22:23);
theta_meta    = kvector(24:27);

% ----- Initial conditions ------------------------------------------------
Y0 = [G0 r0A*dN0 r0B*dN0 r0A*dN0 r0B*dN0 aA0 aB0 0 0];

%% ===== SIMULATE MODEL ===================================================
[T, Z] = ode15s( @(T, Z) ode2MEv1strains(T, Z, theta_biomass, theta_yA, theta_yB, theta_meta, omega), [0, tmax], Y0, odeoptions);

% ----- Get species -------------------------------------------------------
Gluc = Z(:,1); % glucose concentration in the culture vessel
yA   = Z(:,2); % A+/B- total population density of strain A
yB   = Z(:,3); % A-/B+ total population density of strain B
yAa  = Z(:,4); % A+/B- total active population density of strain A
yBa  = Z(:,5); % A-/B+ total active population density of strain B
A    = Z(:,6); % A concentration in the culture vessel
B    = Z(:,7); % B concentration in the culture vessel
X    = Z(:,8); % X concentration in the culture vessel
Y    = Z(:,9); % Y concenttation in the cutlure vessel
yNa  = yAa + yBa;

% ----- Iterate over T ----------------------------------------------------
dY_by_dt     = zeros(length(T),length(Y0));
J0_grow      = zeros(length(T),2);
J0_upt_gluc  = zeros(length(T),2);
J0_leak_yX_X = zeros(length(T),2);
J0_upt_yX_Y  = zeros(length(T),2);
J0_meta      = zeros(length(T),2);
for t = 1:length(T)
    [dY_by_dt(t,:), J0_grow(t,:), J0_upt_gluc(t,:), J0_leak_yX_X(t,:), J0_upt_yX_Y(t,:), J0_meta(t,:)] = ode2MEv1strains(T(t), Z(t,:), theta_biomass, theta_yA, theta_yB, theta_meta, omega);
end

% ----- Find where gluc < 1 -----------------------------------------------
tdx = sum(Gluc > 0.001*Gluc(1));

% --- Calculate volumetric productivity -----------------------------------
vProd = (Y(end) - Y(1))/(T(tdx) - T(1));

% --- Calcualte yield -----------------------------------------------------
pYield = Y(end)/Gluc(1);

%% ===== RETURN PERFORMANCE ===============================================
performancevector( 1) = yNa(end);               performancenames{ 1} = 'N(end)';
performancevector( 2) = T(tdx);                 performancenames{ 2} = 'T(tdx)';
performancevector( 3) = vProd;                  performancenames{ 3} = 'vProd';
performancevector( 4) = pYield;                 performancenames{ 4} = 'pYield';
performancevector( 5) = yAa(end)./yNa(end);     performancenames{ 5} = 'yAa(end)_by_yNa(end)';
performancevector( 6) = yBa(end)./yNa(end);     performancenames{ 6} = 'yBa(end)_by_yNa(end)';
performancevector( 7) = max(J0_grow(:,1));      performancenames{ 7} = 'max(J_grow_yA)';
performancevector( 8) = max(J0_grow(:,2));      performancenames{ 8} = 'max(J_grow_yB)';
performancevector( 9) = max(J0_leak_yX_X(:,1)); performancenames{ 9} = 'max(J_leak_yA_A)';
performancevector(10) = max(J0_leak_yX_X(:,2)); performancenames{10} = 'max(J_leak_yB_B)';
performancevector(11) = max(J0_upt_yX_Y(:,1));  performancenames{11} = 'max(J_upt_yB_A)';
performancevector(12) = max(J0_upt_yX_Y(:,2));  performancenames{12} = 'max(J_upt_yA_B)';
performancevector(13) = max(J0_meta(:,1));      performancenames{13} = 'max(J_conv_yA_GX)';
performancevector(14) = max(J0_meta(:,2));      performancenames{14} = 'max(J_conv_yB_XY)';


end