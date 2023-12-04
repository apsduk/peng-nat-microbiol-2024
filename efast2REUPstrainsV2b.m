function [performancevector, performancenames, T, yNa, Y] = efast2REUPstrainsV2b(X, kvector, kidx, tmax, odeoptions)

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
theta_yA      = kvector(10:17);
theta_yB      = kvector(18:25);
omega         = kvector(26:27);

% ----- Initial conditions ------------------------------------------------
Y0 = [G0 r0A*dN0 r0B*dN0 r0A*dN0 r0B*dN0 aA0 aB0];

%% ===== SIMULATE MODEL ===================================================
[T, Y] = ode15s( @(T, Z) ode2REUPstrains(T, Z, theta_biomass, theta_yA, theta_yB, omega), [0, tmax], Y0, odeoptions);

% ----- Get species -------------------------------------------------------
Gluc = Y(:,1); % glucose concentration in the culture vessel
yA   = Y(:,2); % A+/B- total population density of strain A
yB   = Y(:,3); % A-/B+ total population density of strain B
yAa  = Y(:,4); % A+/B- total active population density of strain A
yBa  = Y(:,5); % A-/B+ total active population density of strain B
A    = Y(:,6); % A conconcentration in the culture vessel
B    = Y(:,7); % B concentration in the culture vessel
yNa  = yAa + yBa;

% ----- Iterate over T ----------------------------------------------------
dY_by_dt     = zeros(length(T),length(Y0));
J0_grow      = zeros(length(T),2);
J0_upt_gluc  = zeros(length(T),2);
J0_leak_yX_X = zeros(length(T),2);
J0_upt_yX_Y  = zeros(length(T),2);
for t = 1:length(T)
    [dY_by_dt(t,:), J0_grow(t,:), J0_upt_gluc(t,:), J0_leak_yX_X(t,:), J0_upt_yX_Y(t,:)] = ode2REUPstrains(T(t), Y(t,:), theta_biomass, theta_yA, theta_yB, omega);
end

% ----- Find where gluc < 1 -----------------------------------------------
tdx = sum(Gluc > 0.001*Gluc(1));

%% ===== RETURN PERFORMANCE ===============================================
performancevector( 1) = yNa(end);                performancenames{ 1} = 'N(end)';
performancevector( 2) = T(tdx);                  performancenames{ 2} = 'T(tdx)';
performancevector( 3) = yAa(end)./yNa(end);      performancenames{ 3} = 'yAa(end)_by_yNa(end)';
performancevector( 4) = yBa(end)./yNa(end);      performancenames{ 4} = 'yBa(end)_by_yNa(end)';
performancevector( 5) = max(J0_grow(:,1));       performancenames{ 5} = 'max(J_grow_yA)';
performancevector( 6) = max(J0_grow(:,2));       performancenames{ 6} = 'max(J_grow_yB)';
performancevector( 7) = max(J0_leak_yX_X(:,1));  performancenames{ 7} = 'max(J_leak_yA_A)';
performancevector( 8) = max(J0_leak_yX_X(:,2));  performancenames{ 8} = 'max(J_leak_yB_B)';
performancevector( 9) = max(J0_upt_yX_Y(:,1));   performancenames{ 9} = 'max(J_upt_A)';
performancevector(10) = max(J0_upt_yX_Y(:,2));   performancenames{10} = 'max(J_upt_B)';

end