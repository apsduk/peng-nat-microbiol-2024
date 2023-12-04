function [performancevector, performancenames, T, yNa, Y] = efast1strainV2(X, kvector, kidx, tmax, odeoptions)

%% ===== SUBSTITUTE PARAMETERS DEFINED BY kidx ============================
kvector(kidx) = X;

%% ===== GET SPECIFIC VECTORS OUT OF INPUT ================================
% ----- Initial conditions ------------------------------------------------
dN0 = kvector(1); % total initial population
G0  = kvector(2); % glucose conc.
aA0 = kvector(3); % amino acid A
aB0 = kvector(4); % amino acid B

% ------ Other parameters -------------------------------------------------
theta_biomass = kvector( 5: 6);
theta_yA      = kvector( 7:12);
omega         = kvector(13);

% ----- Initial conditions ------------------------------------------------
Y0 = [G0 dN0 dN0 aA0 aB0];

%% ===== SIMULATE MODEL ===================================================
[T, Y] = ode15s( @(T, Z) ode1strain(T, Z, theta_biomass, theta_yA, omega), [0, tmax], Y0, odeoptions);

% ----- Get species -------------------------------------------------------
Gluc = Y(:,1); % glucose concentration in the culture vessel
yA   = Y(:,2); % A+/B- total population density of strain A
yAa  = Y(:,3); % A+/B- total active population density of strain A
A    = Y(:,4); % A conconcentration in the culture vessel
B    = Y(:,5); % B concentration in the culture vessel
yNa  = yA;

% ----- Iterate over T ----------------------------------------------------
dY_by_dt     = zeros(length(T),length(Y0));
J0_grow      = zeros(length(T),1);
J0_upt_gluc  = zeros(length(T),1);
J0_leak_yX_X = zeros(length(T),1);
J0_upt_yX_Y  = zeros(length(T),1);
for t = 1:length(T)
    [dY_by_dt(t,:), J0_grow(t), J0_upt_gluc(t), J0_leak_yX_X(t), J0_upt_yX_Y(t)] = ode1strain(T(t), Y(t,:), theta_biomass, theta_yA, omega);
end

% ----- Find where gluc < 1 -----------------------------------------------
tdx = sum(Gluc > 0.001*Gluc(1));

%% ===== RETURN PERFORMANCE ===============================================
performancevector( 1) = yNa(end);                performancenames{ 1} = 'N(end)';
performancevector( 2) = T(tdx);                  performancenames{ 2} = 'T(tdx)';
performancevector( 3) = yAa(end)./yNa(end);      performancenames{ 3} = 'yAa(end)_by_yNa(end)';
performancevector( 4) = max(J0_grow(:,1));       performancenames{ 4} = 'max(J_grow_yA)';
performancevector( 5) = max(J0_leak_yX_X(:,1));  performancenames{ 5} = 'max(J_leak_yA_A)';
performancevector( 6) = max(J0_upt_yX_Y(:,1));   performancenames{ 6} = 'max(J_upt_yA_B)';

end