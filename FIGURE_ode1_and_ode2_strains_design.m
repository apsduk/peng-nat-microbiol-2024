%% %%%%% SET UP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear('all'); close('all')

% ----- Time span ---------------------------------------------------------
tmax = 72;

% ----- Initial conditions ------------------------------------------------
G0 = 20; % g per L
yA0 = 0.03;
aA0 = 0; % mg per L
aB0 = 75; % mg per L

% ----- M-M parameters ----------------------------------------------------
kmG = 5; % g per L
kmB = 50; % mg per L

% --- Biomass parameters --------------------------------------------------
gamma_gluc = 0.05; % 0.025; % Biomass yield of S. cerevisiae grown on glucose
gamma_B    = 0.5; % 0.025; % Biomass yield of S. cerevisiae grown on B
% --- Strain A parameters (Genotype A+/B-) --------------------------------
Vmax_gluc_yA = 7.2; % Maximum glucose uptake rate
Km_gluc_yA   = kmG; % Michaelis constant for glucose uptake
eta_yA       = 1e-3;   % Mortality rate constant of strain A
Vmax_yA_B    = 30;  % Maximum B uptake rate by strain A
Km_yA_B      = kmB; % Michaelis constant for B uptake
delta_yA_A   = 1;   % Number of B molecules produced per molecule of glucose is consumed

% ----- Build vectors -----------------------------------------------------
theta_biomass = [gamma_gluc gamma_B];
theta_yA      = [Vmax_gluc_yA Km_gluc_yA eta_yA Vmax_yA_B Km_yA_B delta_yA_A];

% ----- Set up options ----------------------------------------------------
odeoptions = odeset('NonNegative', 1:7);

% ----- Leaky vectors -----------------------------------------------------
omega(1) = 0;

%% ===== SIMULATE CO-CULTURE ==============================================
v_phi = linspace(0, 1, 21);

cmap = parula(length(v_phi));

fplot = figure('Units', 'centimeters', 'Position', [5 5 40 20]);

% ----- Iterate over phi --------------------------------------------------
for i_a = 1:length(v_phi)

    % ----- Update omega ----------------------------------------------
    omega = v_phi(i_a);

    % ----- Update initial conditions ---------------------------------
    Z0 = [G0 yA0 yA0 aA0 aB0];

    % ----- Simulate model --------------------------------------------
    [T, Z] = ode15s( @(T, Z) ode1strain(T, Z, theta_biomass, theta_yA, omega), [0, tmax], Z0);

    % ----- Get species -----------------------------------------------
    Gluc = Z(:,1); % glucose concentration in the culture vessel
    yA   = Z(:,2); % A+/B- total population density of strain A
    yAa  = Z(:,3); % A+/B- total active population density of strain A
    A    = Z(:,4); % A conconcentration in the culture vessel
    B    = Z(:,5); % B concentration in the culture vessel

    % ----- Iterate over T ----------------------------------------------------
    dZ_by_dt     = zeros(length(T),length(Z0));
    J0_grow      = zeros(length(T),1);
    J0_upt_gluc  = zeros(length(T),1);
    J0_leak_yX_X = zeros(length(T),1);
    J0_upt_yX_Y  = zeros(length(T),1);
    for t = 1:length(T)
        [dZ_by_dt(t,:), J0_grow(t,:), J0_upt_gluc(t,:), J0_leak_yX_X(t,:), J0_upt_yX_Y(t,:)] = ode1strain(T(t), Z(t,:), theta_biomass, theta_yA, omega);
    end

    % ----- Plot results ----------------------------------------------
    figure(fplot.Number);
     %subplot(1, 3, 1); hold('on');
    % plot(T, yAa, '-', 'LineWidth', 4, 'Color', cmap(i_a,:));
    % xlabel('Time(h)');
    % ylabel('y_N = y_1 + y_2');
    % set(gca, 'PlotBoxAspectRatio', [1 1 1], 'Box', 'on', 'LineWidth', 2, 'FontSize', 16, 'FontWeight', 'bold');
    % xticks([0:24:72]); xlim([0 72]);
    subplot(1, 3, 2); hold('on');
    plot(T, J0_leak_yX_X, '-', 'LineWidth', 4, 'Color', cmap(i_a,:));
    ylabel('J^{y_1}_{leak,1}');
    xlabel('Time (h)');
    set(gca, 'PlotBoxAspectRatio', [1 1 1], 'Box', 'on', 'LineWidth', 2, 'FontSize', 16, 'FontWeight', 'bold');
    xticks([0:24:72]); xlim([0 72]);

    output(i_a,1) = max(J0_grow)
    output(i_a,2) = max(J0_leak_yX_X);

end

subplot(1, 3, 3);
hold('on');
plot(output(:,1), output(:,2), '-k', 'LineWidth', 2);
for i_a = 1:length(v_phi)
    plot(output(i_a,1), output(i_a,2), 'o', 'LineWidth', 4, 'Color', cmap(i_a,:), 'MarkerFaceColor', cmap(i_a,:), 'MarkerSize', 8);
end
set(gca, 'PlotBoxAspectRatio', [1 1 1], 'Box', 'on', 'LineWidth', 2, 'FontSize', 16, 'FontWeight', 'bold');
ylabel('J^{y_1}_{leak,1}');
xlabel('J_{grow}');
saveas(fplot.Number, 'FIG_1_B', 'svg');

%% ===== SIMULATE CO-CULTURE ==============================================
phi_yA_A = linspace(0.01, 0.99, 50);
phi_yB_B = linspace(0.01, 0.99, 50);

fid = 0;
omega = [];
theta_biomass = [];
theta_yA = [];
theta_yB = [];
yB0 = yA0;
aA0 = aB0;
Z0 = [];
for ia = 1:length(phi_yA_A)
    for ib = 1:length(phi_yB_B)
        
        % --- Update fid ID -----------------------------------------------
        fid = fid + 1;
        
        % --- Update omega ------------------------------------------------
        omega(fid,:) = [phi_yA_A(ia) phi_yB_B(ib)];
        
        % --- Update biomass parameters -----------------------------------
        theta_biomass(fid,:) = [gamma_gluc gamma_B gamma_B];
        
        % --- Update yA strain parameters ---------------------------------
        theta_yA(fid,:) = [Vmax_gluc_yA Km_gluc_yA eta_yA Vmax_yA_B Km_yA_B delta_yA_A];
        
        % --- Update yB strain parameters ---------------------------------
        theta_yB(fid,:) = [Vmax_gluc_yA Km_gluc_yA eta_yA Vmax_yA_B Km_yA_B delta_yA_A];
        
        % --- Update initial conditions -----------------------------------
        if sum(omega(fid,:)) == 0; Z0(fid,:) = [G0   0   0   0   0   0   0];
        elseif omega(fid,2)  == 0; Z0(fid,:) = [G0 yA0   0 yA0   0   0 aB0];
        elseif omega(fid,1)  == 0; Z0(fid,:) = [G0   0 yB0   0 yB0 aA0   0];
        else;                      Z0(fid,:) = [G0 yA0 yB0 yA0 yB0   0   0];
        end
        
    end
end

imax = length(omega(:,1));

%% ===== SIMULATE MODEL ===================================================
for i = 1:imax
    
    % ----- Simulate model --------------------------------------------
    [T, Z] = ode15s( @(T, Z) ode2strains(T, Z, theta_biomass(i,:), theta_yA(i,:), theta_yB(i,:), omega(i,:)), [0, tmax], Z0(i,:));
    
    % ----- Get species -----------------------------------------------
    Gluc = Z(:,1); % glucose concentration in the culture vessel
    yA   = Z(:,2); % A+/B- total population density of strain A
    yB   = Z(:,3); % A-/B+ total population density of strain B
    yAa  = Z(:,4); % A+/B- total active population density of strain A
    yBa  = Z(:,5); % A-/B+ total active population density of strain B
    % A    = Z(:,6); % A conconcentration in the culture vessel
    % B    = Z(:,7); % B concentration in the culture vessel
    
    % ----- Plot results ----------------------------------------------
    %figure(fplot.Number); fid = fid + 1;
    %subplot(length(v_phi), length(v_phi), fid); pbaspect([1 1 1]);
    %yyaxis('left');  plot(T, yA + yB, '-k'); set(gca, 'YColor', 'k');
    %yyaxis('right'); plot(T, yA./(yA + yB), '-r', T, yB./(yA + yB), '-b'); set(gca, 'YColor', 'k');
    %title(['phi_yA_B = ',num2str(v_phi(ia)), ' | phi_yB_A = ',num2str(v_phi(ib))], 'Interpreter', 'none', 'FontSize', 8);
    
    % --- Return output population ----------------------------------------
    yNa(i) = yAa(end) + yBa(end);
    yAa_by_yNa(i) = yAa(end)./(yAa(end) + yBa(end));
    A(i) = Z(end,6);
    B(i) = Z(end,7);
    
    
end

omega_yA = reshape(omega(:,1), [length(phi_yB_B), length(phi_yA_A)]);
omega_yB = reshape(omega(:,2), [length(phi_yB_B), length(phi_yA_A)]);
yNa = reshape(yNa, [length(phi_yB_B), length(phi_yA_A)]);
yAa_by_yNa = reshape(yAa_by_yNa, [length(phi_yB_B), length(phi_yA_A)]);
A = reshape(A, [length(phi_yB_B), length(phi_yA_A)]);
B = reshape(B, [length(phi_yB_B), length(phi_yA_A)]);
B(B<0) = 0;

%% ===== PLOT =============================================================
fplot = figure('Units', 'centimeters', 'Position', [5 5 40 20]);
subplot(1, 2, 1);
surf(omega_yA, omega_yB, yNa, 'LineStyle', 'none');
set(gca, 'View', [0 90]);
xlabel('\phi_{1}'); ylabel('\phi_{2}');
set(gca, 'PlotBoxAspectRatio', [1 1 1], 'Box', 'on', 'LineWidth', 2, 'FontSize', 16, 'FontWeight', 'bold');
cbar = colorbar;
ylabel(cbar,'y_N = y_1 + y_2','FontSize',16,'Rotation',90);
subplot(1, 2, 2);
surf(omega_yA, omega_yB, yAa_by_yNa, 'LineStyle', 'none');
set(gca, 'View', [0 90]);
xlabel('\phi_{1}'); ylabel('\phi_{2}');
set(gca, 'PlotBoxAspectRatio', [1 1 1], 'Box', 'on', 'LineWidth', 2, 'FontSize', 16, 'FontWeight', 'bold');
cbar = colorbar;
ylabel(cbar,'y_1/(y_1 + y_2)','FontSize',16,'Rotation',90);
saveas(fplot.Number, 'FIG_1_C_D', 'svg');
