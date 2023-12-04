function [performancevector, performancenames, T, yNa, Y] = efast3v1bstrainsV2(X, kvector, kidx, tmax, odeoptions)

%% ===== SUBSTITUTE PARAMETERS DEFINED BY kidx ============================
kvector(kidx) = X;

%% ===== GET SPECIFIC VECTORS OUT OF INPUT ================================
% ----- Initial conditions ------------------------------------------------
N0  = kvector(1); % total initial population
r0A = kvector(2); % A ratio
r0B = kvector(3); % B ratio
r0C = kvector(4); % C ratio
G0  = kvector(5); % glucose conc.
aA0 = kvector(6); % amino acid A
aB0 = kvector(7); % amino acid B
aC0 = kvector(8); % amino acid B

% ------ Other parameters -------------------------------------------------
theta_biomass = kvector( 9:12);
theta_yA      = kvector(13:19);
theta_yB      = kvector(20:25);
theta_yC      = kvector(26:31);
omega         = kvector(32:35);

% ----- Initial conditions ------------------------------------------------
yA0 = r0A*N0;
yB0 = r0B*N0;
yC0 = r0C*N0;

Y0 = [G0 yA0 yB0 yC0 yA0 yB0 yC0 aA0 aB0 aC0];

%% ===== SIMULATE MODEL ===================================================
[T, Y] = ode15s( @(T, Z) ode3v1bstrains(T, Z, theta_biomass, theta_yA, theta_yB, theta_yC, omega), [0, tmax], Y0, odeoptions);

% ----- Get species -------------------------------------------------------
Gluc = Y(:, 1); % glucose concentration in the culture vessel
yA   = Y(:, 2); % A+/B- total population density of strain A
yB   = Y(:, 3); % B+/C- total population density of strain B
yC   = Y(:, 4); % C+/A- total population density of strain C
yAa  = Y(:, 5); % A+/B- total active population density of strain A
yBa  = Y(:, 6); % B+/C- total active population density of strain B
yCa  = Y(:, 7); % C+/A- total active population density of strain C
A    = Y(:, 8); % A concentration in the culture vessel
B    = Y(:, 9); % B concentration in the culture vessel
C    = Y(:,10); % C concentration in the culture vessel
yNa  = yAa + yBa + yCa;

% ----- Iterate over T ----------------------------------------------------
dY_by_dt     = zeros(length(T),length(Y0));
J0_grow      = zeros(length(T),3);
J0_upt_gluc  = zeros(length(T),3);
J0_leak_yX_X = zeros(length(T),3);
J0_upt_yX_Y  = zeros(length(T),3);
for t = 1:length(T)
    [dY_by_dt(t,:), J0_grow(t,:), J0_upt_gluc(t,:), J0_leak_yX_X(t,:), J0_upt_yX_Y(t,:)] = ode3v1bstrains(T(t), Y(t,:), theta_biomass, theta_yA, theta_yB, theta_yC, omega);
end

% ----- Find where gluc < 1 -----------------------------------------------
tdx = sum(Gluc > 0.001*Gluc(1));

%% ===== RETURN PERFORMANCE ===============================================
performancevector( 1) = yNa(end);                performancenames{ 1} = 'N(end)';
performancevector( 2) = T(tdx);                  performancenames{ 2} = 'T(tdx)';
performancevector( 3) = yAa(end)./yNa(end);      performancenames{ 3} = 'yAa(end)_by_yNa(end)';
performancevector( 4) = yBa(end)./yNa(end);      performancenames{ 4} = 'yBa(end)_by_yNa(end)';
performancevector( 5) = yCa(end)./yNa(end);      performancenames{ 5} = 'yCa(end)_by_yNa(end)';
performancevector( 6) = max(J0_grow(:,1));       performancenames{ 6} = 'max(J_grow_yA)';
performancevector( 7) = max(J0_grow(:,2));       performancenames{ 7} = 'max(J_grow_yB)';
performancevector( 8) = max(J0_grow(:,3));       performancenames{ 8} = 'max(J_grow_yC)';
performancevector( 9) = max(J0_leak_yX_X(:,1));  performancenames{ 9} = 'max(J_leak_yA_A)';
performancevector(10) = max(J0_leak_yX_X(:,2));  performancenames{10} = 'max(J_leak_yB_B)';
performancevector(11) = max(J0_leak_yX_X(:,3));  performancenames{11} = 'max(J_leak_yC_C)';
performancevector(12) = max(J0_upt_yX_Y(:,1));   performancenames{12} = 'max(J_upt_yX_A)';
performancevector(13) = max(J0_upt_yX_Y(:,2));   performancenames{13} = 'max(J_upt_yX_B)';
performancevector(14) = max(J0_upt_yX_Y(:,3));   performancenames{14} = 'max(J_upt_yX_C)';

end